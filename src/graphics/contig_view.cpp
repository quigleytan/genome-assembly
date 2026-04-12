/*
 * contig_view.cpp
 */

#include "assembler/graphics/contig_view.h"

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
    }

    // Build scaffold rows sorted by length descending
    scaffoldRows_.clear();
    scaffoldRows_.reserve(session_.scaffolds.size());

    for (size_t si = 0; si < session_.scaffolds.size(); ++si) {
        const VisScaffold& vs = session_.scaffolds[si];

        ScaffoldRow row;
        row.scaffoldIndex       = si;
        row.isCircular          = vs.isCircular;
        row.sortedContigIndices = vs.contigIndices;

        std::sort(row.sortedContigIndices.begin(),
                  row.sortedContigIndices.end(),
                  [&](size_t a, size_t b) {
                      return session_.contigs[a].sequence.length() >
                             session_.contigs[b].sequence.length();
                  });

        scaffoldRows_.push_back(std::move(row));
    }

    // Collect unscaffolded contigs sorted longest-first
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

    // If no animation steps, make everything visible immediately
    if (session_.contigSteps.empty()) {
        for (auto& ds : displayStates_) {
            ds.visible      = true;
            ds.fillFraction = 1.0f;
        }
    }

    // Build genome map segment list
    buildGenomeSegments();
}

// GENOME MAP

void ContigView::buildGenomeSegments() {
    genomeSegments_.clear();

    const std::string& genome = session_.genomeSequence;
    if (genome.empty()) return;

    size_t scaffoldCount = 0;
    size_t i = 0;

    while (i < genome.size()) {
        if (genome[i] == 'N') {
            // Gap - consume run of Ns
            size_t start = i;
            while (i < genome.size() && genome[i] == 'N')
                ++i;

            GenomeSegment seg;
            seg.type     = GenomeSegment::Type::Gap;
            seg.index    = 0;
            seg.startPos = start;
            seg.length   = i - start;
            genomeSegments_.push_back(seg);

        } else {
            // Scaffold - consume run of non-N chars
            size_t start = i;
            while (i < genome.size() && genome[i] != 'N')
                ++i;

            GenomeSegment seg;
            seg.type     = GenomeSegment::Type::Scaffold;
            seg.index    = scaffoldCount;
            seg.startPos = start;
            seg.length   = i - start;
            genomeSegments_.push_back(seg);
            ++scaffoldCount;
        }
    }
}

// COLOR HELPERS

uint32_t ContigView::hsvToImCol32(float h, float s, float v, float a) {
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
        return IM_COL32(130, 130, 130, 255);

    float hue = fmod(static_cast<float>(scaffoldIndex) * 0.618033988f, 1.0f);
    float sat = isCircular ? 0.40f : 0.80f;
    float val = 0.88f;
    return hsvToImCol32(hue, sat, val);
}

uint32_t ContigView::scaffoldColor(size_t scaffoldIndex) const {
    if (scaffoldIndex >= session_.scaffolds.size())
        return IM_COL32(130, 130, 130, 255);
    bool isCircular = session_.scaffolds[scaffoldIndex].isCircular;
    return assignColor(static_cast<int>(scaffoldIndex), isCircular);
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
            break;
    }
}

void ContigView::resetAnimation() {
    // Clears all animation state values
    currentStep_ = 0;
    accumulator_ = 0.0f;
    playing_     = false;
    for (auto& ds : displayStates_) {
        ds.visible       = false;
        ds.fillFraction  = 0.0f;
        ds.basesAppended = 0;
    }
    if (session_.contigSteps.empty()) {
        for (auto& ds : displayStates_) {
            ds.visible      = true;
            ds.fillFraction = 1.0f;
        }
    }
}

void ContigView::seekToStep(size_t targetStep) {
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

// GENOME MAP TAB

void ContigView::renderGenomeMapTab(float availHeight) {
    if (session_.genomeSequence.empty()) {
        ImGui::TextDisabled("No genome sequence available.");
        ImGui::TextDisabled("Run the assembly pipeline with visualization enabled.");
        return;
    }

    const float totalWidth  = ImGui::GetContentRegionAvail().x;

    const float detailWidth = totalWidth * DETAIL_WIDTH_FRACTION;
    const float seqWidth    = totalWidth - detailWidth - ImGui::GetStyle().ItemSpacing.x;

    // Wrapping sequence view in a child window to enable scrolling
    ImGui::BeginChild("##genomemap",
                      ImVec2(seqWidth, availHeight),
                      false,
                      ImGuiWindowFlags_HorizontalScrollbar);

    ImDrawList* drawList    = ImGui::GetWindowDrawList();
    const std::string& genome = session_.genomeSequence;

    // Bases per row based on available width and cell size
    const int basesPerRow = std::max(1, static_cast<int>(
        (seqWidth - 16.0f) / CELL_WIDTH));
    const float rowHeight = CELL_HEIGHT + 2.0f;

    const size_t totalBases = genome.size();
    const size_t totalRows  = (totalBases + basesPerRow - 1) / basesPerRow;

    // Reserve vertical space for scrolling
    ImGui::Dummy(ImVec2(seqWidth - 16.0f,
                        static_cast<float>(totalRows) * rowHeight + 8.0f));

    ImVec2 origin = ImGui::GetItemRectMin();
    origin.x += 8.0f;
    origin.y += 4.0f;

    size_t segIdx = 0;

    for (size_t row = 0; row < totalRows; ++row) {
        size_t rowStart = row * static_cast<size_t>(basesPerRow);
        size_t rowEnd   = std::min(rowStart + static_cast<size_t>(basesPerRow), totalBases);

        for (size_t pos = rowStart; pos < rowEnd; ++pos) {

            // Advance to the segment that contains this position
            while (segIdx + 1 < genomeSegments_.size() &&
                   pos >= genomeSegments_[segIdx].startPos + genomeSegments_[segIdx].length)
                ++segIdx;

            const GenomeSegment& seg = genomeSegments_[segIdx];

            float cellX = origin.x + static_cast<float>(pos - rowStart) * CELL_WIDTH;
            float cellY = origin.y + static_cast<float>(row) * rowHeight;
            ImVec2 cellMin = { cellX,              cellY            };
            ImVec2 cellMax = { cellX + CELL_WIDTH, cellY + CELL_HEIGHT };

            // Choose color
            uint32_t color;
            if (seg.type == GenomeSegment::Type::Gap) {
                color = IM_COL32(80, 80, 80, 255); // dark gray for N gaps
            } else {
                color = scaffoldColor(seg.index);
                // Dim segments whose scaffold hasn't appeared in animation yet
                if (seg.index < session_.scaffolds.size()) {
                    const VisScaffold& vs = session_.scaffolds[seg.index];
                    bool anyVisible = false;
                    for (size_t ci : vs.contigIndices) {
                        if (ci < displayStates_.size() && displayStates_[ci].visible) {
                            anyVisible = true;
                            break;
                        }
                    }
                    if (!anyVisible)
                        color = IM_COL32(40, 40, 40, 255);
                }
            }

            drawList->AddRectFilled(cellMin, cellMax, color);

            // Grid lines between cells
            drawList->AddRect(cellMin, cellMax,
                              IM_COL32(20, 20, 20, 60), 0.0f, 0, 0.5f);

            // Selection highlight
            if (static_cast<int>(segIdx) == selectedSegment_) {
                drawList->AddRect(cellMin, cellMax,
                                  IM_COL32(255, 255, 255, 200), 0.0f, 0, 1.5f);
            }
        }
    }

    // Click detection
    for (size_t si = 0; si < genomeSegments_.size(); ++si) {
        const GenomeSegment& seg = genomeSegments_[si];

        size_t pos  = seg.startPos;
        size_t row  = pos / static_cast<size_t>(basesPerRow);
        size_t col  = pos % static_cast<size_t>(basesPerRow);

        float cellX = origin.x + static_cast<float>(col) * CELL_WIDTH;
        float cellY = origin.y + static_cast<float>(row) * rowHeight;

        // Width covers this segment's extent on its first row only
        size_t firstRowEnd = std::min(
            seg.startPos + seg.length,
            (row + 1) * static_cast<size_t>(basesPerRow));
        float buttonW = static_cast<float>(firstRowEnd - pos) * CELL_WIDTH;

        ImGui::SetCursorScreenPos({ cellX, cellY });
        ImGui::InvisibleButton(
            ("##seg" + std::to_string(si)).c_str(),
            ImVec2(buttonW, CELL_HEIGHT));

        if (ImGui::IsItemClicked()) {
            selectedSegment_ = static_cast<int>(si);
            selectedContig_  = -1; // clear animation tab selection
        }

        if (ImGui::IsItemHovered()) {
            ImGui::BeginTooltip();
            if (seg.type == GenomeSegment::Type::Gap) {
                ImGui::Text("Gap");
                ImGui::Text("Length:   %zu Ns", seg.length);
                ImGui::Text("Position: %zu", seg.startPos);
            } else {
                ImGui::Text("Scaffold %zu", seg.index);
                if (seg.index < session_.scaffolds.size())
                    ImGui::Text("Contigs:  %zu",
                        session_.scaffolds[seg.index].contigIndices.size());
                ImGui::Text("Length:   %zu bases", seg.length);
                ImGui::Text("Position: %zu", seg.startPos);
            }
            ImGui::EndTooltip();
        }
    }

    ImGui::EndChild();

    ImGui::SameLine();

    // Right detail panel
    ImGui::BeginChild("##mapdetail", ImVec2(detailWidth, availHeight), true);

    if (selectedSegment_ < 0 ||
        selectedSegment_ >= static_cast<int>(genomeSegments_.size())) {
        ImGui::TextDisabled("Click a segment\nto see details.");
        ImGui::EndChild();
        return;
    }

    const GenomeSegment& seg = genomeSegments_[selectedSegment_];

    if (seg.type == GenomeSegment::Type::Gap) {
        ImGui::TextColored(ImVec4(0.7f, 0.7f, 0.7f, 1.0f), "Gap Region");
        ImGui::Separator();
        ImGui::Text("Length:   %zu Ns", seg.length);
        ImGui::Text("Position: %zu", seg.startPos);
        ImGui::Spacing();
        ImGui::TextDisabled("Gap status: unknown");
        ImGui::TextDisabled("(gap filling not yet run)");

    } else {
        ImGui::TextColored(ImVec4(1.0f, 0.85f, 0.4f, 1.0f),
                           "Scaffold %zu", seg.index);
        ImGui::Separator();
        ImGui::Text("Length:   %zu bases", seg.length);
        ImGui::Text("Position: %zu", seg.startPos);

        if (seg.index < session_.scaffolds.size()) {
            const VisScaffold& vs = session_.scaffolds[seg.index];
            ImGui::Text("Circular: %s", vs.isCircular ? "yes" : "no");
            ImGui::Spacing();
            ImGui::TextDisabled("Contigs:");

            for (size_t i = 0; i < vs.contigIndices.size(); ++i) {
                size_t ci = vs.contigIndices[i];
                if (ci >= session_.contigs.size()) continue;
                const VisContig& c = session_.contigs[ci];

                ImGui::Spacing();
                ImGui::Text("  Contig %zu", ci);
                ImGui::Text("    Length: %zu bases", c.sequence.length());
                ImGui::Text("    Score:  %.4f", c.score);

                if (i < vs.gaps.size() && vs.gaps[i] == -1)
                    ImGui::TextColored(ImVec4(0.6f, 0.6f, 0.6f, 1.0f),
                                       "    [gap follows]");
            }
        }
    }

    ImGui::EndChild();
}

//SUB-RENDERERS

void ContigView::renderContigBar(ImDrawList* drawList, ImVec2 origin, size_t contigIdx) {
    const ContigDisplayState& ds = displayStates_[contigIdx];

    const float fullW   = barWidth(contigIdx);
    const float filledW = fullW * ds.fillFraction;

    // Background track always renders allowing unbuilt contigs to be seen
    drawList->AddRectFilled(
        origin,
        ImVec2(origin.x + fullW, origin.y + BAR_HEIGHT),
        IM_COL32(50, 50, 50, 180), 3.0f);

    if (!ds.visible) {
        ImGui::SetCursorScreenPos(origin);
        ImGui::InvisibleButton(
            ("##contig" + std::to_string(contigIdx)).c_str(),
            ImVec2(fullW, BAR_HEIGHT));
        return;
    }

    if (filledW > 0.0f) {
        drawList->AddRectFilled(
            origin,
            ImVec2(origin.x + filledW, origin.y + BAR_HEIGHT),
            ds.color, 3.0f);
    }

    if (static_cast<int>(contigIdx) == selectedContig_) {
        drawList->AddRect(
            ImVec2(origin.x - 1, origin.y - 1),
            ImVec2(origin.x + fullW + 1, origin.y + BAR_HEIGHT + 1),
            IM_COL32(255, 255, 255, 220), 3.0f, 0, 2.0f);
    }

    ImGui::SetCursorScreenPos(origin);
    ImGui::InvisibleButton(
        ("##contig" + std::to_string(contigIdx)).c_str(),
        ImVec2(fullW, BAR_HEIGHT));

    if (ImGui::IsItemClicked()) {
        selectedContig_  = static_cast<int>(contigIdx);
        selectedSegment_ = -1; // clear genome map selection
    }

    if (ImGui::IsItemHovered()) {
        ImGui::BeginTooltip();
        ImGui::Text("Contig %zu", contigIdx);
        ImGui::Text("Length: %zu bases", session_.contigs[contigIdx].sequence.length());
        ImGui::Text("Score:  %.3f", session_.contigs[contigIdx].score);
        ImGui::EndTooltip();
    }
}

void ContigView::renderBarPanel(float availWidth, float availHeight) {
    // Creating a scrollable child region for the contig bars
    ImGui::BeginChild("##barview",
                      ImVec2(availWidth, availHeight),
                      false,
                      ImGuiWindowFlags_HorizontalScrollbar);

    ImDrawList* drawList = ImGui::GetWindowDrawList();
    ImVec2 cursor        = ImGui::GetCursorScreenPos();
    float  x0            = cursor.x + 8.0f;

    auto advanceCursor = [&](float height) {
        ImGui::Dummy(ImVec2(BAR_MAX_WIDTH + 16.0f, height));
        cursor = ImGui::GetCursorScreenPos();
        cursor.x = x0;
    };

    // Rendering scaffolds
    for (const ScaffoldRow& row : scaffoldRows_) {
        const VisScaffold& vs = session_.scaffolds[row.scaffoldIndex];

        ImGui::SetCursorScreenPos(cursor);
        ImGui::TextColored(
            ImVec4(0.7f, 0.7f, 0.7f, 1.0f),
            "Scaffold %zu  (%zu contigs%s)",
            row.scaffoldIndex,
            row.sortedContigIndices.size(),
            row.isCircular ? ", circular" : "");
        advanceCursor(ImGui::GetTextLineHeight() + 4.0f);

        for (size_t ci : row.sortedContigIndices) {
            ImVec2 barOrigin = { x0, cursor.y };
            renderContigBar(drawList, barOrigin, ci);
            advanceCursor(BAR_HEIGHT + BAR_SPACING);

            for (size_t ei = 0; ei < vs.contigIndices.size(); ++ei) {
                if (vs.contigIndices[ei] == ci && ei + 1 < vs.gaps.size()) {
                    if (vs.gaps[ei] == -1) {
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
        advanceCursor(SCAFFOLD_GAP);
    }

    // Rendering unscaffolded contigs
    if (!unscaffoldedContigs_.empty()) {
        ImGui::SetCursorScreenPos(cursor);
        ImGui::TextColored(
            ImVec4(0.5f, 0.5f, 0.5f, 1.0f),
            "Unscaffolded  (%zu contigs)",
            unscaffoldedContigs_.size());
        advanceCursor(ImGui::GetTextLineHeight() + 4.0f);

        // Rendering each individual contig
        for (size_t ci : unscaffoldedContigs_) {
            ImVec2 barOrigin = { x0, cursor.y };
            renderContigBar(drawList, barOrigin, ci);
            advanceCursor(BAR_HEIGHT + BAR_SPACING);
        }
    }

    ImGui::EndChild();
}

void ContigView::renderDetailPanel(float panelWidth, float availHeight) {
    ImGui::BeginChild("##detail", ImVec2(panelWidth, availHeight), true);

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

    // Contig metadata
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

    // Sequence preview information
    ImGui::TextDisabled("Sequence preview:");

    constexpr size_t PREVIEW_LEN = 120;
    if (c.sequence.length() <= PREVIEW_LEN) {
        ImGui::TextWrapped("%s", c.sequence.c_str());
    } else {
        std::string preview = c.sequence.substr(0, PREVIEW_LEN) + "...";
        ImGui::TextWrapped("%s", preview.c_str());
    }

    ImGui::Spacing();

    // User selects copy sequence function
    if (ImGui::Button("Copy sequence"))
        ImGui::SetClipboardText(c.sequence.c_str());

    ImGui::EndChild();
}

void ContigView::renderTimeline() {
    ImGui::Separator();

    const size_t totalSteps = session_.contigSteps.size();

    if (ImGui::Button(playing_ ? "  Pause  " : "  Play   ")) {
        if (currentStep_ >= totalSteps)
            resetAnimation();
        playing_ = !playing_;
    }

    ImGui::SameLine();

    if (ImGui::Button("Restart"))
        resetAnimation();

    ImGui::SameLine();

    ImGui::Text("Step %zu / %zu", currentStep_, totalSteps);

    if (totalSteps == 0) {
        ImGui::SameLine();
        ImGui::TextDisabled("(no animation data)");
        return;
    }

    int sliderStep = static_cast<int>(currentStep_);
    int sliderMax  = static_cast<int>(totalSteps);

    ImGui::SetNextItemWidth(-180.0f);
    if (ImGui::SliderInt("##scrub", &sliderStep, 0, sliderMax)) {
        playing_ = false;
        seekToStep(static_cast<size_t>(sliderStep));
    }

    ImGui::SameLine();

    ImGui::TextDisabled("Speed:");
    for (int i = 0; i < NUM_PRESETS; ++i) {
        ImGui::SameLine();
        bool isSelected = (i == speedPresetIndex_);
        if (isSelected)
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.3f, 0.6f, 0.9f, 1.0f));
        if (ImGui::SmallButton(SPEED_LABELS[i]))
            speedPresetIndex_ = i;
        if (isSelected)
            ImGui::PopStyleColor();
    }
}

// RENDER - TOP LEVEL

void ContigView::render() {
    ImGuiIO& io = ImGui::GetIO();

    float menuBarHeight = ImGui::GetFrameHeight();
    ImGui::SetNextWindowPos(ImVec2(0, menuBarHeight));
    ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x,
                                    io.DisplaySize.y - menuBarHeight));

    ImGui::Begin("contig_view",
                 nullptr,
                 ImGuiWindowFlags_NoTitleBar  |
                 ImGuiWindowFlags_NoResize    |
                 ImGuiWindowFlags_NoMove      |
                 ImGuiWindowFlags_NoScrollbar |
                 ImGuiWindowFlags_NoBringToFrontOnFocus);

    // Reserve space at the bottom for the timeline
    const float timelineHeight = ImGui::GetTextLineHeightWithSpacing() * 3.5f
                               + ImGui::GetStyle().ItemSpacing.y * 2.0f
                               + 4.0f;     const float tabBarHeight   = ImGui::GetTextLineHeightWithSpacing() + 8.0f;

    const float availHeight    = ImGui::GetContentRegionAvail().y
                                 - timelineHeight - tabBarHeight;

    if (ImGui::BeginTabBar("##viewtabs")) {

        if (ImGui::BeginTabItem("Assembly Animation")) {
            const float totalWidth    = ImGui::GetContentRegionAvail().x;
            const float detailWidth   = totalWidth * DETAIL_WIDTH_FRACTION;
            const float barPanelWidth = totalWidth - detailWidth
                                        - ImGui::GetStyle().ItemSpacing.x;
            renderBarPanel(barPanelWidth, availHeight);
            ImGui::SameLine();
            renderDetailPanel(detailWidth, availHeight);
            ImGui::EndTabItem();
        }

        if (ImGui::BeginTabItem("Genome Map")) {
            renderGenomeMapTab(availHeight);
            ImGui::EndTabItem();
        }

        ImGui::EndTabBar();
    }

    renderTimeline();

    ImGui::End();
}