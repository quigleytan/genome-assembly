/*
 * gap_filling.h
 * Created by Tanner Quigley on 3/27/26
 * Summary:
 * - Attempts to bridge an UNKNOWN_GAP between two scaffolded contigs by
 *   performing a bounded local search ("local reassembly") over the
 *   already-built DeBruijnGraph for a path from the upstream contig's end
 *   node to the downstream contig's start node.
 * - Falls back to an N-run sized by GapEstimator when no such path exists
 *   within the search bound - e.g. the region truly was not sequenced, or
 *   the resolution strategy pruned the connecting branch during scaffolding.
 * Important notes:
 * - "Local" refers to the search being depth-bounded (sized off of
 *   GapEstimator's estimate) rather than a full Eulerian re-walk - it reuses
 *   the graph already built from all reads instead of re-assembling from
 *   raw reads.
 * - A successful reassembly is not necessarily biologically unique; the
 *   shortest bridging path found by breadth-first search is used.
 */

#ifndef GAP_FILLING_H
#define GAP_FILLING_H

#include <cstddef>
#include <optional>
#include <string>

#include "assembler/core/kmer_types.h"
#include "assembler/core/de_bruijn_graph.h"
#include "assembler/construction/gap_estimation.h"

// SEARCH CONFIGURATION
//
// Kept as a free struct (not nested in GapFiller) for the same reason as
// GapEstimationConfig - see the comment in gap_estimation.h.

struct GapFillConfig {
    size_t searchPadding  = 20;   // Extra depth beyond the gap estimate, to absorb estimation error.
    size_t maxSearchDepth = 1000; // Hard cap on local reassembly BFS depth, regardless of estimate.
};

class GapFiller {

public:

    struct Result {
        std::string sequence;         // Bridging bases (excludes both flanking contigs' own
                                       // sequence); empty when resolved is false.
        bool        resolved = false; // True if a connecting path was found via local reassembly.
        size_t      gapEstimate = 0;  // GapEstimator's estimate - sizes the N-fallback when
                                       // resolved is false, and bounds the search depth.
    };

    /**
     * @brief Constructor for GapFiller.
     * @param graph     Fully constructed DeBruijnGraph (the same one ContigAssembler walked).
     * @param estimator GapEstimator sharing the same k-mer coverage data, used to size
     *                  both the BFS search bound and the N-run fallback.
     * @param config    Tunable search parameters.
     */
    GapFiller(const DeBruijnGraph& graph, const GapEstimator& estimator, GapFillConfig config = GapFillConfig());

    /**
     * @brief Attempts to resolve the gap between two contigs.
     *
     * First estimates the gap length via estimator_, then performs a
     * bounded BFS in graph_ from upstreamEndNode to downstreamStartNode,
     * capped at min(estimate + searchPadding, maxSearchDepth) steps. If a
     * path is found, it is returned as the bridging sequence with
     * resolved = true. Otherwise the returned Result carries only the
     * estimate, so the caller can pad with that many N's instead.
     *
     * @param upstreamSeq         Full sequence of the upstream contig (for estimation).
     * @param upstreamEndNode     Encoded k-1 mer where the upstream contig ends.
     * @param downstreamSeq       Full sequence of the downstream contig (for estimation).
     * @param downstreamStartNode Encoded k-1 mer where the downstream contig starts.
     * @return Result describing the outcome (see struct docs above).
     */
    [[nodiscard]] Result fillGap(const std::string& upstreamSeq, NodeId upstreamEndNode,
                                  const std::string& downstreamSeq, NodeId downstreamStartNode) const;

private:

    const DeBruijnGraph& graph_;
    const GapEstimator&  estimator_;
    GapFillConfig config_;

    /**
     * @brief Bounded breadth-first search for a path from `from` to `to` in graph_.
     *
     * Each graph edge represents appending one base (the last character of
     * the neighbor's decoded k-1 mer label). Every node is visited at most
     * once, so the search is O(local nodes/edges within maxDepth) - no
     * combinatorial blowup from revisiting nodes along cycles or repeats.
     *
     * @param from     Node to search from (upstream contig's end node).
     * @param to       Node to search for (downstream contig's start node).
     * @param maxDepth Maximum number of edges to traverse before giving up.
     * @return The bridging base sequence (excluding the seed k-1 mer) if
     *         `to` was reached within maxDepth steps; std::nullopt otherwise.
     */
    [[nodiscard]] std::optional<std::string> localReassembly(NodeId from, NodeId to, size_t maxDepth) const;
};

#endif //GAP_FILLING_H
