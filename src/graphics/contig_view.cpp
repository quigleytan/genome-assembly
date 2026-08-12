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

    const size_t scaffoldCount = session_.scaffolds.size();

    // gapFills holds exactly one entry per junction between consecutive
    // scaffolds (see VisSession::gapFills in vis_data.h), so segment
    // boundaries can be derived directly from scaffold/gap lengths instead
    // of scanning genomeSequence for 'N' runs. This also correctly surfaces
    // gaps that GapFiller resolved with real bridging sequence - those
    // contain no 'N' characters at all, so a content scan would have
    // silently merged them into the neighboring scaffold.
    if (scaffoldCount > 0 && session_.gapFills.size() == scaffoldCount - 1) {
        size_t pos = 0;
        for (size_t si = 0; si < scaffoldCount; ++si) {
            size_t scaffoldLen = 0;
            for (size_t ci : session_.scaffolds[si].contigIndices)
                if (ci < session_.contigs.size())
                    scaffoldLen += session_.contigs[ci].sequence.length();

            GenomeSegment seg;
            seg.type     = GenomeSegment::Type::Scaffold;
            seg.index    = si;
            seg.startPos = pos;
            seg.length   = scaffoldLen;
            genomeSegments_.push_back(seg);
            pos += scaffoldLen;

            if (si + 1 < scaffoldCount) {
                const VisGapFill& gf = session_.gapFills[si];

                GenomeSegment gapSeg;
                gapSeg.type         = GenomeSegment::Type::Gap;
                gapSeg.index        = si;
                gapSeg.startPos     = pos;
                gapSeg.length       = gf.filledLength;
                gapSeg.resolved     = gf.resolved;
                gapSeg.estimatedGap = gf.estimatedGap;
                genomeSegments_.push_back(gapSeg);
                pos += gf.filledLength;
            }
        }
        return;
    }

    // Fallback for sessions with no gap-fill data (e.g. a single scaffold,
    // or a mismatch): treat the whole genome as one Scaffold segment rather
    // than guessing junction positions.
    GenomeSegment seg;
    seg.type     = GenomeSegment::Type::Scaffold;
    seg.index    = 0;
    seg.startPos = 0;
    seg.length   = genome.length();
    genomeSegments_.push_back(seg);
}

size_t ContigView::findSegmentIndex(size_t pos) const {
    if (genomeSegments_.empty()) return 0;

    // genomeSegments_ is sorted by startPos and covers [0, genome length)
    // contiguously with no gaps between entries, so "last segment with
    // startPos <= pos" is always the (unique) segment containing pos.
    auto it = std::upper_bound(
        genomeSegments_.begin(), genomeSegments_.end(), pos,
        [](size_t value, const GenomeSegment& seg) { return value < seg.startPos; });

    if (it == genomeSegments_.begin()) return 0;
    --it;
    return static_cast<size_t>(std::distance(genomeSegments_.begin(), it));
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

void ContigView::addBases(size_t contigIndex, size_t n) {
    if (contigIndex >= displayStates_.size() || n == 0) return;

    ContigDisplayState& ds = displayStates_[contigIndex];
    ds.basesAppended += n;

    size_t totalBases = session_.contigs[contigIndex].sequence.length();
    if (totalBases > 0)
        ds.fillFraction = static_cast<float>(ds.basesAppended)
                        / static_cast<float>(totalBases);
}

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

        case TraversalStep::Type::BaseAppended:
            addBases(idx, step.count);
            break;

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
    currentStep_  = 0;
    accumulator_  = 0.0f;
    stepProgress_ = 0;
    playing_      = false;
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

    // Speed presets are bases/second. A BaseAppended step can now cover many
    // bases at once (see Recorder::baseAppended), so playback advances in
    // base-units rather than one discrete step per frame-tick: structural
    // events (start/finish) apply instantly and free, while a BaseAppended
    // step's count is drained incrementally via stepProgress_ so the bar
    // still visibly grows base-by-base instead of jumping in one step.
    const float basesPerSecond = SPEED_PRESETS[speedPresetIndex_];
    accumulator_ += deltaTime * basesPerSecond;

    while (accumulator_ > 0.0f && currentStep_ < session_.contigSteps.size()) {
        const TraversalStep& step = session_.contigSteps[currentStep_];

        if (step.type != TraversalStep::Type::BaseAppended) {
            applyStep(step);
            ++currentStep_;
            stepProgress_ = 0;
            continue;
        }

        const size_t remaining = step.count - stepProgress_;
        const size_t take = std::min(remaining, static_cast<size_t>(accumulator_));

        if (take == 0) break; // Not enough accumulated budget for one more base yet.

        addBases(step.contigIndex, take);
        stepProgress_ += take;
        accumulator_  -= static_cast<float>(take);

        if (stepProgress_ >= step.count) {
            stepProgress_ = 0;
            ++currentStep_;
        }
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

    // Only draw/hit-test rows currently scrolled into view. genomeSequence
    // can run to millions of bases, and this tab used to iterate every row
    // (draw pass) and every segment-row span (click-detection pass)
    // unconditionally every frame regardless of scroll position - fine for
    // a few thousand bases, unusable well before a few million.
    constexpr size_t ROW_MARGIN = 2; // extra rows drawn beyond the visible edge, to avoid pop-in
    const float scrollY = ImGui::GetScrollY();
    size_t firstVisibleRow = (scrollY > 0.0f) ? static_cast<size_t>(scrollY / rowHeight) : 0;
    firstVisibleRow = (firstVisibleRow > ROW_MARGIN) ? firstVisibleRow - ROW_MARGIN : 0;

    const size_t rowsInView   = static_cast<size_t>(availHeight / rowHeight) + 2 * ROW_MARGIN + 2;
    const size_t lastVisibleRow = (totalRows == 0) ? 0
        : std::min(totalRows - 1, firstVisibleRow + rowsInView);

    if (totalRows > 0) {
        size_t rowStartPos = firstVisibleRow * static_cast<size_t>(basesPerRow);
        size_t segIdx = findSegmentIndex(std::min(rowStartPos, totalBases - 1));

        for (size_t row = firstVisibleRow; row <= lastVisibleRow; ++row) {
            size_t rowStart = row * static_cast<size_t>(basesPerRow);
            size_t rowEnd   = std::min(rowStart + static_cast<size_t>(basesPerRow), totalBases);
            if (rowStart >= rowEnd) break;

            size_t pos = rowStart;
            while (pos < rowEnd) {

                // Advance to the segment that contains this position
                while (segIdx + 1 < genomeSegments_.size() &&
                       pos >= genomeSegments_[segIdx].startPos + genomeSegments_[segIdx].length)
                    ++segIdx;

                const GenomeSegment& seg = genomeSegments_[segIdx];
                size_t spanEnd = std::min(rowEnd, seg.startPos + seg.length);
                if (spanEnd <= pos) { ++pos; continue; } // defensive: never stall on a malformed segment

                // Choose color
                uint32_t color;
                if (seg.type == GenomeSegment::Type::Gap) {
                    // Amber = bridged via local reassembly, dark gray = N-padded (unresolved)
                    color = seg.resolved
                        ? IM_COL32(214, 168, 60, 255)
                        : IM_COL32(80, 80, 80, 255);
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

                bool isSelected = (static_cast<int>(segIdx) == selectedSegment_);

                // Draw each cell in this segment's span within the row
                for (size_t cellPos = pos; cellPos < spanEnd; ++cellPos) {
                    float cellX = origin.x + static_cast<float>(cellPos - rowStart) * CELL_WIDTH;
                    float cellY = origin.y + static_cast<float>(row) * rowHeight;
                    ImVec2 cellMin = { cellX,              cellY            };
                    ImVec2 cellMax = { cellX + CELL_WIDTH, cellY + CELL_HEIGHT };

                    drawList->AddRectFilled(cellMin, cellMax, color);
                    drawList->AddRect(cellMin, cellMax,
                                      IM_COL32(20, 20, 20, 60), 0.0f, 0, 0.5f);
                    if (isSelected)
                        drawList->AddRect(cellMin, cellMax,
                                          IM_COL32(255, 255, 255, 200), 0.0f, 0, 1.5f);
                }

                // One click/hover region for this segment's span within the row
                float spanX = origin.x + static_cast<float>(pos - rowStart) * CELL_WIDTH;
                float spanY = origin.y + static_cast<float>(row) * rowHeight;
                float spanW = static_cast<float>(spanEnd - pos) * CELL_WIDTH;

                ImGui::SetCursorScreenPos({ spanX, spanY });
                ImGui::InvisibleButton(
                    ("##seg" + std::to_string(segIdx) + "_r" + std::to_string(row)).c_str(),
                    ImVec2(spanW, CELL_HEIGHT));

                if (ImGui::IsItemClicked()) {
                    selectedSegment_ = static_cast<int>(segIdx);
                    selectedContig_  = -1;
                }

                if (ImGui::IsItemHovered()) {
                    ImGui::BeginTooltip();
                    if (seg.type == GenomeSegment::Type::Gap) {
                        ImGui::Text(seg.resolved ? "Gap (resolved)" : "Gap (unresolved)");
                        ImGui::Text("Length:    %zu bases", seg.length);
                        ImGui::Text("Estimated: %zu bases", seg.estimatedGap);
                        ImGui::Text("Position:  %zu", seg.startPos);
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

                pos = spanEnd;
            }
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
        ImGui::TextColored(
            seg.resolved ? ImVec4(0.90f, 0.72f, 0.30f, 1.0f) : ImVec4(0.7f, 0.7f, 0.7f, 1.0f),
            seg.resolved ? "Gap (Resolved)" : "Gap (Unresolved)");
        ImGui::Separator();
        ImGui::Text("Length:    %zu bases", seg.length);
        ImGui::Text("Estimated: %zu bases", seg.estimatedGap);
        ImGui::Text("Position:  %zu", seg.startPos);
        ImGui::Spacing();

        if (seg.resolved)
            ImGui::TextWrapped(
                "Bridged via bounded local reassembly over the de Bruijn "
                "graph - a real connecting path was found between the "
                "flanking contigs, so this region is not N-padded.");
        else
            ImGui::TextWrapped(
                "No connecting path was found within the search bound. "
                "Padded with N's, sized to the k-mer frequency drop "
                "estimate at the flanking contig boundaries.");

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

                if (i < vs.gaps.size() && vs.gaps[i] == -1) {
                    bool haveGapInfo = seg.index < session_.gapFills.size();
                    bool resolved    = haveGapInfo && session_.gapFills[seg.index].resolved;
                    ImGui::TextColored(
                        resolved ? ImVec4(0.90f, 0.72f, 0.30f, 1.0f) : ImVec4(0.6f, 0.6f, 0.6f, 1.0f),
                        "    [gap follows%s]",
                        haveGapInfo ? (resolved ? " - resolved" : " - unresolved") : "");
                }
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
        }

        // Trailing gap indicator: a scaffold's last contig (in walk order)
        // always carries an UNKNOWN_GAP marker (see ScaffoldEntry) - this is
        // the junction GapEstimator/GapFiller resolved between this scaffold
        // and whichever one follows it in the pseudo-genome.
        if (!vs.gaps.empty() && vs.gaps.back() == -1) {
            bool   haveGapInfo  = row.scaffoldIndex < session_.gapFills.size();
            bool   resolved     = haveGapInfo && session_.gapFills[row.scaffoldIndex].resolved;
            size_t estimatedGap = haveGapInfo ? session_.gapFills[row.scaffoldIndex].estimatedGap : 0;
            size_t filledLength = haveGapInfo ? session_.gapFills[row.scaffoldIndex].filledLength : 0;

            uint32_t gapColor = resolved
                ? IM_COL32(214, 168, 60, 220)
                : IM_COL32(100, 100, 100, 180);

            ImVec2 gapOrigin = { x0, cursor.y };
            drawList->AddRectFilled(
                gapOrigin,
                ImVec2(gapOrigin.x + GAP_BAR_WIDTH, gapOrigin.y + BAR_HEIGHT * 0.5f),
                gapColor);

            ImGui::SetCursorScreenPos(gapOrigin);
            ImGui::InvisibleButton(
                ("##gap" + std::to_string(row.scaffoldIndex)).c_str(),
                ImVec2(GAP_BAR_WIDTH, BAR_HEIGHT * 0.5f));

            if (ImGui::IsItemHovered()) {
                ImGui::BeginTooltip();
                if (haveGapInfo) {
                    ImGui::Text(resolved ? "Gap: resolved via local reassembly"
                                          : "Gap: unresolved (N-padded)");
                    ImGui::Text("Estimated: %zu bases", estimatedGap);
                    ImGui::Text("Filled:    %zu bases", filledLength);
                } else {
                    ImGui::Text("Gap: no resolution data available");
                }
                ImGui::EndTooltip();
            }

            advanceCursor(BAR_HEIGHT * 0.5f + BAR_SPACING);
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