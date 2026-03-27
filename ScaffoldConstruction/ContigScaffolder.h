/*
 * ContigScaffolder.h
 * Created by Tanner Quigley
 * Summary:
 * - Orders contigs into scaffolds using shared boundary node relationships.
 * - Supports multiple branch resolution strategies for comparison.
 * Important notes:
 * - Operates on contigs produced by ContigTraversal.
 * - Scaffolds are ordered lists of contig indices with gap estimates.
 * - Ambiguous branch points are handled according to ResolutionStrategy.
 * - KmerFrequency and OverlapQuality scoring require an optional KmerTable.
 *   If none is provided, those metrics are skipped and weight falls to Length.
 */

#ifndef CONTIG_SCAFFOLDER_H
#define CONTIG_SCAFFOLDER_H

#include <vector>
#include <string>

#include "KmerTypes.h"
#include "ScaffoldConstruction/DeBruijnGraph.h"
#include "ScaffoldConstruction/ContigTraversal.h"
#include "DataProcessing/KmerTable.h"
#include "DataProcessing/OpenAddressingTable.h"

// SCORING WEIGHT CONFIGURATIONS

struct ScoringWeights {
    double lengthWeight         = 1.0;
    double kmerFrequencyWeight  = 0.0;
    double overlapQualityWeight = 0.0;

    // Convenience presets
    static ScoringWeights lengthOnly()    { return {1.0, 0.0, 0.0}; }
    static ScoringWeights frequencyOnly() { return {0.0, 1.0, 0.0}; }
    static ScoringWeights combined()      { return {0.4, 0.4, 0.2}; }
};

// RESOLUTION STRATEGY HANDLING

struct ResolutionStrategy {
    bool skipAmbiguous = true;
    ScoringWeights weights;

    // Convenience presets
    static ResolutionStrategy skip() {
        ResolutionStrategy s;
        s.skipAmbiguous = true;
        s.weights = ScoringWeights::lengthOnly();
        return s;
    }

    static ResolutionStrategy greedy() {
        ResolutionStrategy s;
        s.skipAmbiguous = false;
        s.weights = ScoringWeights::lengthOnly();
        return s;
    }

    static ResolutionStrategy scored() {
        ResolutionStrategy s;
        s.skipAmbiguous = false;
        s.weights = ScoringWeights::combined();
        return s;
    }
};

// SCAFFOLD ENTRY STRUCT

struct ScaffoldEntry {
    static constexpr int UNKNOWN_GAP = -1;
    static constexpr int DIRECT_OVERLAP = 0;

    static constexpr int NOT_APPLICABLE = -1;

    size_t contigIndex; // Index into the contigs vector
    int gapAfter;       // Bases between this contig and the next.
                        // UNKNOWN_GAP (-1) = unresolved, DIRECT_OVERLAP (0) = direct overlap
    int score = NOT_APPLICABLE;          // NOT_APPLICABLE (-1)
};

// SCAFFOLD STRUCT

struct Scaffold {
    std::vector<ScaffoldEntry> entries;
    bool isCircular = false;

    [[nodiscard]] size_t contigCount() const { return entries.size(); }
};

// SCAFFOLDER CLASS

class ContigScaffolder {

private:

    const std::vector<ContigTraversal::Contig>& contigs_;
    const DeBruijnGraph& graph_;
    const KmerTable* kmerTable_;
    size_t maxContigLength_ = 0;
    ResolutionStrategy strategy_;
    std::vector<Scaffold> scaffolds_;

    // Connection maps
    OpenAddressingTable<NodeId, std::vector<size_t>> endNodeMap_;
    OpenAddressingTable<NodeId, std::vector<size_t>> startNodeMap_;

    /**
     * @brief Populates endNodeMap_ and startNodeMap_ from the contig list.
     *
     * Iterates through all contigs and inserts each contig's index into
     * endNodeMap_ keyed by its endNode, and into startNodeMap_ keyed by
     * its startNode.
     */
    void buildConnectionMap();

    /**
     * @brief Computes a length-based score for a contig.
     *
     * Longer contigs typically have greater k-mer support and fewer unresolved
     * ambiguities, making them more likely to represent correct assemblies.
     * This function assigns higher scores to longer contigs by dividing the contig's length
     * by the max length of all contigs.
     *
     * @param contig Contig to be evaluated.
     * @return Net score based on contig length from [0, 1], to be weighted.
     */
    [[nodiscard]] double computeLengthScore(const ContigTraversal::Contig& contig) const;

    /**
     * @brief Computes a normalized frequency-based score for a contig.
     *
     * Evaluates the level of support for a contig by examining the frequencies (counts)
     * of its constituent k-mers in the k-mer table.
     * The average k-mer frequency is computed across the contig, then normalized
     * by an upper bound defined as (mean + 2 * standard deviation) of the observed
     * k-mer counts. This reduces the influence of extreme high-frequency k-mers
     * (e.g., repeats) while still favoring consistently well-supported contigs.
     *
     * @param contig Contig to be evaluated.
     * @return Normalized net frequency score from [0, 1], to be weighted.
     */
    [[nodiscard]] double computeFrequencyScore(const ContigTraversal::Contig& contig) const;

    /**
     * @brief Computes an overlap consistency score for a contig.
     *
     * Measures how well the suffix of the contig's sequence
     * aligns with the expected sequence encoded by its terminal graph node.
     * The expected suffix is decoded from the contig's end node and compared
     * base-by-base with the actual suffix of the contig sequence.
     *
     * @param contig Contig to be evaluated.
     * @return Net overlap consistency score from [0, 1], to be weighted.
     */
    [[nodiscard]] double computeOverlapScore(const ContigTraversal::Contig& contig) const;

    /**
     * @brief Scores a candidate contig for branch resolution.
     *
     * Computes a weighted score using strategy_.weights. Metrics requiring
     * kmerTable_ are skipped if kmerTable_ is null, and their weight is
     * redistributed to the Length metric.
     *
     * @param contigIndex Index of the candidate contig.
     * @return Weighted score — higher means more likely to be the correct next contig.
     */
    [[nodiscard]] double scoreContig(size_t contigIndex) const;

    /**
     * @brief Resolves the next contig index from a given boundary node.
     *
     * If strategy_.skipAmbiguous is true and more than one contig starts
     * at boundaryNode, returns npos. Otherwise scores all unvisited candidates
     * using scoreContig() and returns the highest scoring one.
     *
     * @param boundaryNode The endNode of the current contig.
     * @param visited      Tracks which contigs are already assigned.
     * @return Index into contigs_ of the next contig, or npos if unresolvable.
     */
    [[nodiscard]] size_t resolveNext(NodeId boundaryNode, const std::vector<bool>& visited) const;

    /**
     * @brief Walks the connection map from startIndex to build one scaffold.
     *
     * Follows endNode -> startNode links until a dead end, ambiguous branch
     * point (if skipAmbiguous), or already-visited contig is reached.
     * Each step appends a ScaffoldEntry with gapAfter = DIRECT_OVERLAP for
     * connected contigs and gapAfter = UNKNOWN_GAP at unresolved terminations.
     *
     * @param startIndex Index of the contig to begin the scaffold from.
     * @param visited    Tracks which contigs are already assigned.
     * @return Populated Scaffold.
     */
    [[nodiscard]] Scaffold walkScaffold(size_t startIndex,
                                        std::vector<bool>& visited) const;

    /**
     * @brief Returns true if a contig is a valid scaffold entry point.
     *
     * A contig qualifies if its startNode has no entries in endNodeMap_,
     * meaning no other contig ends where this one begins. These are the
     * natural starting points for linear scaffolds.
     *
     * @param contigIndex Index of the contig to check.
     * @return True if no other contig ends at this contig's startNode.
     */
    [[nodiscard]] bool isScaffoldStart(size_t contigIndex) const;

public:

    /**
     * @brief Constructor for ContigScaffolder.
     *
     * @param contigs    Contigs produced by ContigTraversal.
     * @param graph      DeBruijnGraph used to build the contigs.
     * @param strategy   Branch resolution strategy (default: ResolutionStrategy::skip()).
     * @param kmerTable  Optional KmerTable for frequency-based scoring. Pass nullptr to disable.
     */
    ContigScaffolder(const std::vector<ContigTraversal::Contig>& contigs,
                     const DeBruijnGraph& graph,
                     ResolutionStrategy strategy = ResolutionStrategy::skip(),
                     const KmerTable* kmerTable = nullptr);

    /**
     * @brief Builds all scaffolds from the contig connection map.
     *
     * Phase 1: Builds a connection map from contig boundary nodes.
     * Phase 2: Walks from all scaffold entry points, producing linear scaffolds.
     * Phase 3: Handles any remaining unvisited contigs as isolated or circular
     *          scaffolds with no external entry point.
     */
    void buildScaffolds();

    /**
     * @brief Returns the list of scaffolds.
     * @return Const reference to the scaffold vector.
     */
    [[nodiscard]] const std::vector<Scaffold>& getScaffolds() const;

    /**
     * @brief Reports scaffolding statistics to stdout.
     *
     * Strategy:           Resolution strategy used (skip / greedy / scored).
     * Total scaffolds:    Number of scaffolds produced.
     * Circular scaffolds: Number of scaffolds forming a closed loop.
     * Longest scaffold:   Length in bases of the longest scaffold.
     * Unresolved gaps:    Number of gaps with UNKNOWN_GAP size.
     * N50:                Scaffold-level N50.
     */
    void printStats() const;
};

#endif //CONTIG_SCAFFOLDER_H