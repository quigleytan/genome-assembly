/*
 * ContigView.h
 * Summary:
 * - Renders the contig/scaffold bar visualization and animation playback.
 * - Owned by VisualizerApp. Called every frame via update() and render().
 * - All rendering is done through ImGui's DrawList API — no direct OpenGL calls.
 * Important notes:
 * - Bars are sorted by contig length descending within each scaffold group.
 * - Animation interpolates bar width between ContigStarted and ContigFinished
 *   so playback is readable at any speed regardless of sequence length.
 * - selectedContig_ drives the detail panel — -1 means nothing selected.
 */

#ifndef CONTIG_VIEW_H
#define CONTIG_VIEW_H

#include <vector>
#include <string>
#include <cstdint>

#include "imgui.h"
#include "VisData.h"

class ContigView {

public:

    /**
     * @brief Constructs a ContigView for the given session.
     * Precomputes display state, sort order, colors, and bar widths.
     * @param session Fully loaded VisSession. Must outlive this ContigView.
     */
    explicit ContigView(const VisSession& session);

    /**
     * @brief Advances animation state by deltaTime seconds.
     * No-op if not playing or if the step list is exhausted.
     * @param deltaTime Seconds since last frame (from VisualizerApp's frame timer).
     */
    void update(float deltaTime);

    /**
     * @brief Renders the full ContigView panel using ImGui.
     * Must be called between ImGui::NewFrame() and ImGui::Render().
     */
    void render();

private:

    // DATA

    const VisSession& session_;

    // PER-CONTIG DISPLAY STATE

    struct ContigDisplayState {
        bool     visible       = false; // True after ContigStarted step
        float    fillFraction  = 0.0f;  // 0.0 = empty, 1.0 = fully built
        uint32_t color         = 0;     // IM_COL32 RGBA — set in constructor
        size_t   basesAppended = 0;     // Running count for interpolation
    };
    std::vector<ContigDisplayState> displayStates_;

    // One entry per scaffold — sorted contig index list
    struct ScaffoldRow {
        size_t              scaffoldIndex;
        std::vector<size_t> sortedContigIndices; // sorted longest-first
        bool                isCircular = false;
    };
    std::vector<ScaffoldRow> scaffoldRows_;

    // Unscaffolded contigs — sorted longest-first
    std::vector<size_t> unscaffoldedContigs_;

    // ANIMATION STATE

    size_t currentStep_  = 0;
    bool   playing_      = false;
    float  accumulator_  = 0.0f;  // Fractional steps carried between frames

    // Fixed speed presets — steps per second
    static constexpr float SPEED_PRESETS[]    = { 10.0f, 50.0f, 200.0f, 1000.0f };
    static constexpr const char* SPEED_LABELS[] = { "0.5x", "1x", "5x", "50x" };
    static constexpr int   NUM_PRESETS        = 4;
    int    speedPresetIndex_ = 1; // Default: 1x

    // INTERACTION STATE

    int selectedContig_ = -1; // Index into session_.contigs, -1 = none

    // LAYOUT CONSTANTS

    static constexpr float BAR_HEIGHT    = 22.0f;
    static constexpr float BAR_MAX_WIDTH = 580.0f;
    static constexpr float BAR_SPACING   = 5.0f;   // Vertical gap between bars
    static constexpr float SCAFFOLD_GAP  = 16.0f;  // Extra gap between scaffold groups
    static constexpr float GAP_BAR_WIDTH = 10.0f;  // Width of unknown-gap indicator
    static constexpr float DETAIL_WIDTH_FRACTION = 0.30f; // Detail panel takes 30% of window

    // Genome map tab state
    static constexpr float CELL_WIDTH     = 10.0f;  // px per base character
    static constexpr float CELL_HEIGHT    = 18.0f;  // px per row
    static constexpr int   BASES_PER_ROW  = 60;     // wraps at 60 chars like FASTA

    // Segment lookup - maps character position to scaffold or gap index
    struct GenomeSegment {
        enum class Type { Scaffold, Gap };
        Type   type;
        size_t index;        // scaffold index or gap index
        size_t startPos;     // position in genomeSequence
        size_t length;
    };
    std::vector<GenomeSegment> genomeSegments_; // built in buildDisplayData()
    int selectedSegment_ = -1;                  // -1 = none

    // Sub-renderer added in phase 2 for genome map tab
    void renderAnimationTab();
    void renderGenomeMapTab();
    void buildGenomeSegments();
    // PRECOMPUTED DATA

    size_t maxContigLength_ = 1; // Longest contig - denominator for bar widths

    // PRIVATE METHODS

    /**
     * @brief Assigns colors, builds scaffold rows, and sorts contig indices.
     * Called once from the constructor.
     */
    void buildDisplayData();

    /**
     * @brief Assigns an IM_COL32 color to a contig based on scaffold membership.
     * Uses golden-ratio hue spacing so adjacent scaffolds have distinct colors.
     * Circular contigs get a desaturated variant of the same hue.
     * Unscaffolded contigs are gray.
     */
    uint32_t assignColor(int scaffoldIndex, bool isCircular) const;

    uint32_t scaffoldColor(size_t scaffoldIndex) const;

    /**
     * @brief Converts HSV to an IM_COL32 RGBA value.
     * h in [0,1], s in [0,1], v in [0,1].
     */
    static uint32_t hsvToImCol32(float h, float s, float v, float a = 1.0f);

    /**
     * @brief Computes the pixel width for a contig bar.
     * Proportional to sequence length relative to maxContigLength_.
     */
    float barWidth(size_t contigIndex) const;

    /**
     * @brief Applies a single TraversalStep to the display state.
     * Called by update() for each step consumed this frame.
     */
    void applyStep(const TraversalStep& step);

    /**
     * @brief Resets all display state and replays from step 0.
     * Called when the user scrubs the timeline slider or presses restart.
     */
    void resetAnimation();

    /**
     * @brief Seeks animation to a specific step index.
     * Replays all steps from 0 to targetStep by calling applyStep.
     * Used when the user drags the timeline slider.
     */
    void seekToStep(size_t targetStep);

    // Sub-render methods for different UI sections

    /**
     * @brief Renders the scrollable bar view panel (left side).
     * @param availWidth Total width available for the bar view panel.
     * @param availHeight Total height available (excluding timeline bar).
     */
    void renderBarPanel(float availWidth, float availHeight);

    /**
     * @brief Renders one contig bar using ImGui DrawList.
     * @param drawList   ImGui draw list to write into.
     * @param origin     Top-left corner of the bar in screen space.
     * @param contigIdx  Index into session_.contigs.
     */
    void renderContigBar(ImDrawList* drawList, ImVec2 origin, size_t contigIdx);

    /**
     * @brief Renders the detail panel (right side).
     * Shows metadata for selectedContig_, or a placeholder if none selected.
     * @param panelWidth Width of the detail panel.
     * @param availHeight Height available.
     */
    void renderDetailPanel(float panelWidth, float availHeight);

    /**
     * @brief Renders the animation timeline bar at the bottom.
     * Contains play/pause, step counter, scrub slider, and speed presets.
     */
    void renderTimeline();
};

#endif // CONTIG_VIEW_H