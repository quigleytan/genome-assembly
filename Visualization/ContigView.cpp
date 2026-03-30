/*
 * ContigView.cpp
 */

#include "ContigView.h"

#include <algorithm>
#include <cmath>
#include <string>


// STATIC MEMBER DEFINITIONS

constexpr float       ContigView::SPEED_PRESETS[];
constexpr const char* ContigView::SPEED_LABELS[];

// CONSTRUCTION

ContigView::ContigView(const VisSession& session)
    : session_(session) {
    buildDisplayData();
}

void ContigView::buildDisplayData() {
    const size_t n = session_.contigs.size();
    if (n == 0) return;

    // Find max contig length for bar scaling
    maxContigLength_ = 1;
    for (const auto& c : session_.contigs)
        maxContigLength_ = std::max(maxContigLength_, c.sequence.length());

    // Assign per-contig display state
    displayStates_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        const VisContig& c = session_.contigs[i];
        displayStates_[i].color = assignColor(c.scaffoldIndex, c.isCircular);
        // visible and fillFraction start at 0 — animation reveals them
    }

    // Build scaffold rows sorted by length descending
    scaffoldRows_.clear();
    scaffoldRows_.reserve(session_.scaffolds.size());

    for (size_t si = 0; si < session_.scaffolds.size(); ++si) {
        const VisScaffold& vs = session_.scaffolds[si];

        ScaffoldRow row;
        row.scaffoldIndex  = si;
        row.isCircular     = vs.isCircular;
        row.sortedContigIndices = vs.contigIndices; // copy, then sort

        std::sort(row.sortedContigIndices.begin(),
                  row.sortedContigIndices.end(),
                  [&](size_t a, size_t b) {
                      return session_.contigs[a].sequence.length() >
                             session_.contigs[b].sequence.length();
                  });

        scaffoldRows_.push_back(std::move(row));
    }

    // Collect unscaffolded contigs
    unscaffoldedContigs_.clear();
    for (size_t i = 0; i < n; ++i) {
        if (session_.contigs[i].scaffoldIndex < 0)
            unscaffoldedContigs_.push_back(i);
    }
    std::sort(unscaffoldedContigs_.begin(), unscaffoldedContigs_.end(),
              [&](size_t a, size_t b) {
                  return session_.contigs[a].sequence.length() >
                         session_.contigs[b].sequence.length();
              });

    // ── If no animation steps, make everything visible immediately ────────
    if (session_.contigSteps.empty()) {
        for (auto& ds : displayStates_) {
            ds.visible      = true;
            ds.fillFraction = 1.0f;
        }
    }
}

// COLOR HELPERS


uint32_t ContigView::hsvToImCol32(float h, float s, float v, float a) {
    // Standard HSV → RGB conversion
    float r, g, b;
    if (s <= 0.0f) {
        r = g = b = v;
    } else {
        h = fmod(h, 1.0f) * 6.0f;
        int   i  = static_cast<int>(h);
        float f  = h - static_cast<float>(i);
        float p  = v * (1.0f - s);
        float q  = v * (1.0f - s * f);
        float t  = v * (1.0f - s * (1.0f - f));
        switch (i) {
            case 0: r = v; g = t; b = p; break;
            case 1: r = q; g = v; b = p; break;
            case 2: r = p; g = v; b = t; break;
            case 3: r = p; g = q; b = v; break;
            case 4: r = t; g = p; b = v; break;
            default:r = v; g = p; b = q; break;
        }
    }
    return IM_COL32(
        static_cast<int>(r * 255),
        static_cast<int>(g * 255),
        static_cast<int>(b * 255),
        static_cast<int>(a * 255));
}

uint32_t ContigView::assignColor(int scaffoldIndex, bool isCircular) const {
    if (scaffoldIndex < 0)
        return IM_COL32(130, 130, 130, 255); // unscaffolded — gray

    // Golden ratio spacing gives visually distinct hues across scaffold indices
    float hue = fmod(static_cast<float>(scaffoldIndex) * 0.618033988f, 1.0f);
    float sat = isCircular ? 0.40f : 0.80f;
    float val = 0.88f;
    return hsvToImCol32(hue, sat, val);
}

// BAR WIDTH

float ContigView::barWidth(size_t contigIndex) const {
    float ratio = static_cast<float>(session_.contigs[contigIndex].sequence.length())
                / static_cast<float>(maxContigLength_);
    return ratio * BAR_MAX_WIDTH;
}

// ANIMATION


void ContigView::applyStep(const TraversalStep& step) {
    const size_t idx = step.contigIndex;
    if (idx >= displayStates_.size()) return;

    ContigDisplayState& ds = displayStates_[idx];

    switch (step.type) {

        case TraversalStep::Type::ContigStarted:
            ds.visible       = true;
            ds.fillFraction  = 0.0f;
            ds.basesAppended = 0;
            break;

        case TraversalStep::Type::BaseAppended: {
            ++ds.basesAppended;
            size_t totalBases = session_.contigs[idx].sequence.length();
            if (totalBases > 0)
                ds.fillFraction = static_cast<float>(ds.basesAppended)
                                / static_cast<float>(totalBases);
            break;
        }

        case TraversalStep::Type::ContigFinished:
            ds.fillFraction  = 1.0f;
            ds.basesAppended = session_.contigs[idx].sequence.length();
            break;

        default:
            break; // EdgeConsumed / NodeCommitted not used in ContigView
    }
}

void ContigView::resetAnimation() {
    currentStep_ = 0;
    accumulator_ = 0.0f;
    playing_     = false;
    for (auto& ds : displayStates_) {
        ds.visible       = false;
        ds.fillFraction  = 0.0f;
        ds.basesAppended = 0;
    }
    // If no steps, show everything immediately
    if (session_.contigSteps.empty()) {
        for (auto& ds : displayStates_) {
            ds.visible      = true;
            ds.fillFraction = 1.0f;
        }
    }
}

void ContigView::seekToStep(size_t targetStep) {
    // Full replay from scratch — only called on scrub, not every frame
    resetAnimation();
    targetStep = std::min(targetStep, session_.contigSteps.size());
    for (size_t i = 0; i < targetStep; ++i)
        applyStep(session_.contigSteps[i]);
    currentStep_ = targetStep;
}

void ContigView::update(float deltaTime) {
    if (!playing_ || session_.contigSteps.empty()) return;
    if (currentStep_ >= session_.contigSteps.size()) {
        playing_ = false;
        return;
    }

    float stepsPerSecond = SPEED_PRESETS[speedPresetIndex_];
    accumulator_ += deltaTime * stepsPerSecond;

    while (accumulator_ >= 1.0f && currentStep_ < session_.contigSteps.size()) {
        applyStep(session_.contigSteps[currentStep_]);
        ++currentStep_;
        accumulator_ -= 1.0f;
    }

    if (currentStep_ >= session_.contigSteps.size())
        playing_ = false;
}

// RENDERING / SUB-RENDERERS

void ContigView::renderContigBar(ImDrawList* drawList, ImVec2 origin, size_t contigIdx) {
    const ContigDisplayState& ds = displayStates_[contigIdx];
    if (!ds.visible) return;

    const float fullW  = barWidth(contigIdx);
    const float filledW = fullW * ds.fillFraction;

    // Background track (unfilled portion)
    uint32_t bgColor = IM_COL32(50, 50, 50, 180);
    drawList->AddRectFilled(
        origin,
        ImVec2(origin.x + fullW, origin.y + BAR_HEIGHT),
        bgColor,
        3.0f); // corner rounding

    // Filled portion
    if (filledW > 0.0f) {
        drawList->AddRectFilled(
            origin,
            ImVec2(origin.x + filledW, origin.y + BAR_HEIGHT),
            ds.color,
            3.0f);
    }

    // Selection highlight — outline when selected
    if (static_cast<int>(contigIdx) == selectedContig_) {
        drawList->AddRect(
            ImVec2(origin.x - 1, origin.y - 1),
            ImVec2(origin.x + fullW + 1, origin.y + BAR_HEIGHT + 1),
            IM_COL32(255, 255, 255, 220),
            3.0f,
            0,
            2.0f); // line thickness
    }

    ImGui::SetCursorScreenPos(origin);
    ImGui::InvisibleButton(
        ("##contig" + std::to_string(contigIdx)).c_str(),
        ImVec2(fullW, BAR_HEIGHT));

    if (ImGui::IsItemClicked())
        selectedContig_ = static_cast<int>(contigIdx);

    if (ImGui::IsItemHovered()) {
        ImGui::BeginTooltip();
        ImGui::Text("Contig %zu", contigIdx);
        ImGui::Text("Length: %zu bases", session_.contigs[contigIdx].sequence.length());
        ImGui::Text("Score:  %.3f", session_.contigs[contigIdx].score);
        ImGui::EndTooltip();
    }
}

void ContigView::renderBarPanel(float availWidth, float availHeight) {
    ImGui::BeginChild("##barview",
                      ImVec2(availWidth, availHeight),
                      false,
                      ImGuiWindowFlags_HorizontalScrollbar);

    ImDrawList* drawList = ImGui::GetWindowDrawList();
    ImVec2 cursor        = ImGui::GetCursorScreenPos();
    float  x0            = cursor.x + 8.0f; // left margin

    auto advanceCursor = [&](float height) {
        ImGui::Dummy(ImVec2(BAR_MAX_WIDTH + 16.0f, height));
        cursor = ImGui::GetCursorScreenPos();
        cursor.x = x0;
    };

    // SCAFFOLDED CONTIGS
    for (const ScaffoldRow& row : scaffoldRows_) {
        const VisScaffold& vs = session_.scaffolds[row.scaffoldIndex];

        // Scaffold label
        ImGui::SetCursorScreenPos(cursor);
        ImGui::TextColored(
            ImVec4(0.7f, 0.7f, 0.7f, 1.0f),
            "Scaffold %zu  (%zu contigs%s)",
            row.scaffoldIndex,
            row.sortedContigIndices.size(),
            row.isCircular ? ", circular" : "");
        advanceCursor(ImGui::GetTextLineHeight() + 4.0f);

        // Bars within this scaffold
        for (size_t ci : row.sortedContigIndices) {
            ImVec2 barOrigin = { x0, cursor.y };
            renderContigBar(drawList, barOrigin, ci);
            advanceCursor(BAR_HEIGHT + BAR_SPACING);

            // Gap indicator after this contig if applicable
            // Find gap value from scaffold entries
            for (size_t ei = 0; ei < vs.contigIndices.size(); ++ei) {
                if (vs.contigIndices[ei] == ci && ei + 1 < vs.gaps.size()) {
                    int gap = vs.gaps[ei];
                    if (gap == -1) {
                        // Unknown gap — draw a small hatched gray block
                        drawList->AddRectFilled(
                            { x0, cursor.y },
                            { x0 + GAP_BAR_WIDTH, cursor.y + BAR_HEIGHT * 0.5f },
                            IM_COL32(100, 100, 100, 180));
                        advanceCursor(BAR_HEIGHT * 0.5f + BAR_SPACING);
                    }
                    break;
                }
            }
        }

        // Extra spacing between scaffold groups
        advanceCursor(SCAFFOLD_GAP);
    }

    // Unscaffolded contigs
    if (!unscaffoldedContigs_.empty()) {
        ImGui::SetCursorScreenPos(cursor);
        ImGui::TextColored(
            ImVec4(0.5f, 0.5f, 0.5f, 1.0f),
            "Unscaffolded  (%zu contigs)",
            unscaffoldedContigs_.size());
        advanceCursor(ImGui::GetTextLineHeight() + 4.0f);

        for (size_t ci : unscaffoldedContigs_) {
            ImVec2 barOrigin = { x0, cursor.y };
            renderContigBar(drawList, barOrigin, ci);
            advanceCursor(BAR_HEIGHT + BAR_SPACING);
        }
    }

    ImGui::EndChild();
}

void ContigView::renderDetailPanel(float panelWidth, float availHeight) {
    ImGui::BeginChild("##detail",
                      ImVec2(panelWidth, availHeight),
                      true); // show border

    if (selectedContig_ < 0 ||
        selectedContig_ >= static_cast<int>(session_.contigs.size())) {
        ImGui::TextDisabled("Click a contig bar\nto see details.");
        ImGui::EndChild();
        return;
    }

    const VisContig& c = session_.contigs[selectedContig_];

    ImGui::TextColored(ImVec4(1.0f, 0.85f, 0.4f, 1.0f),
                       "Contig %d", selectedContig_);
    ImGui::Separator();

    ImGui::Text("Length:    %zu bases", c.sequence.length());
    ImGui::Text("Scaffold:  %s",
        c.scaffoldIndex >= 0
            ? std::to_string(c.scaffoldIndex).c_str()
            : "none");
    ImGui::Text("Score:     %.4f", c.score);
    ImGui::Text("Circular:  %s", c.isCircular ? "yes" : "no");

    ImGui::Spacing();
    ImGui::TextDisabled("Start node:");
    ImGui::TextWrapped("%s", c.startLabel.c_str());
    ImGui::TextDisabled("End node:");
    ImGui::TextWrapped("%s", c.endLabel.c_str());

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::TextDisabled("Sequence preview:");

    constexpr size_t PREVIEW_LEN = 120;
    if (c.sequence.length() <= PREVIEW_LEN) {
        ImGui::TextWrapped("%s", c.sequence.c_str());
    } else {
        std::string preview = c.sequence.substr(0, PREVIEW_LEN) + "...";
        ImGui::TextWrapped("%s", preview.c_str());
    }

    ImGui::Spacing();
    // Copy full sequence to clipboard
    if (ImGui::Button("Copy sequence")) {
        ImGui::SetClipboardText(c.sequence.c_str());
    }

    ImGui::EndChild();
}

void ContigView::renderTimeline() {
    ImGui::Separator();

    const size_t totalSteps = session_.contigSteps.size();

    // Play / Pause
    if (ImGui::Button(playing_ ? "  Pause  " : "  Play   ")) {
        if (currentStep_ >= totalSteps)
            resetAnimation(); // restart from beginning if at end
        playing_ = !playing_;
    }

    ImGui::SameLine();

    // Restart
    if (ImGui::Button("Restart")) {
        resetAnimation();
    }

    ImGui::SameLine();

    // Step counter
    ImGui::Text("Step %zu / %zu", currentStep_, totalSteps);

    if (totalSteps == 0) {
        ImGui::SameLine();
        ImGui::TextDisabled("(no animation data)");
        return;
    }

    // ── Scrub slider ──────────────────────────────────────────────────────
    int sliderStep = static_cast<int>(currentStep_);
    int sliderMax  = static_cast<int>(totalSteps);

    ImGui::SetNextItemWidth(-180.0f); // leave room for speed buttons on right
    if (ImGui::SliderInt("##scrub", &sliderStep, 0, sliderMax)) {
        playing_ = false; // stop playback when user scrubs
        seekToStep(static_cast<size_t>(sliderStep));
    }

    ImGui::SameLine();

    // ── Speed presets ─────────────────────────────────────────────────────
    ImGui::TextDisabled("Speed:");
    for (int i = 0; i < NUM_PRESETS; ++i) {
        ImGui::SameLine();
        bool isSelected = (i == speedPresetIndex_);
        if (isSelected) {
            ImGui::PushStyleColor(ImGuiCol_Button,
                                  ImVec4(0.3f, 0.6f, 0.9f, 1.0f));
        }
        if (ImGui::SmallButton(SPEED_LABELS[i]))
            speedPresetIndex_ = i;
        if (isSelected)
            ImGui::PopStyleColor();
    }
}

// RENDER — TOP LEVEL

void ContigView::render() {
    // Full-window ImGui window — no title bar, no padding, fills the OS window
    ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(io.DisplaySize);
    ImGui::Begin("ContigView",
                 nullptr,
                 ImGuiWindowFlags_NoTitleBar  |
                 ImGuiWindowFlags_NoResize    |
                 ImGuiWindowFlags_NoMove      |
                 ImGuiWindowFlags_NoScrollbar |
                 ImGuiWindowFlags_NoBringToFrontOnFocus);

    // Reserve space at the bottom for the timeline bar
    const float timelineHeight = ImGui::GetTextLineHeightWithSpacing() * 2.5f;
    const float availHeight    = ImGui::GetContentRegionAvail().y - timelineHeight;
    const float totalWidth     = ImGui::GetContentRegionAvail().x;
    const float detailWidth    = totalWidth * DETAIL_WIDTH_FRACTION;
    const float barPanelWidth  = totalWidth - detailWidth - ImGui::GetStyle().ItemSpacing.x;

    // Left: bar view
    renderBarPanel(barPanelWidth, availHeight);

    ImGui::SameLine();

    // Right: detail panel
    renderDetailPanel(detailWidth, availHeight);

    // Bottom: timeline
    renderTimeline();

    ImGui::End();
}