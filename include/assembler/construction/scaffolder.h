/*
 * scaffolder.h
 * Created by Tanner Quigley 3/15/26
 * Summary:
 * - Orders contigs into scaffolds using shared boundary node relationships.
 * - Supports multiple branch resolution strategies for comparison.
 * Important notes:
 * - Operates on contigs produced by ContigAssembler.
 * - Scaffolds are ordered lists of contig indices with gap estimates.
 * - Ambiguous branch points are handled according to ResolutionStrategy.
 * - KmerFrequency and OverlapQuality scoring require an optional kmer_table.
 *   If none is provided, those metrics are skipped and weight falls to Length.
 */

#ifndef SCAFFOLDER_H
#define SCAFFOLDER_H

#include <vector>

#include "assembler/core/kmer_types.h"
#include "assembler/core/vis_data.h"
#include "assembler/core/de_bruijn_graph.h"
#include "assembler/core/kmer_table.h"
#include "assembler/core/open_addressing_table.h"
#include "assembler/construction/contig_assembler.h"

// SCORING WEIGHT CONFIGURATIONS

struct ScoringWeights {
    double lengthWeight         = 1.0;
    double kmerFrequencyWeight  = 0.0;
    double overlapQualityWeight = 0.0;

    static ScoringWeights lengthOnly()    { return {1.0, 0.0, 0.0}; }
    static ScoringWeights frequencyOnly() { return {0.0, 1.0, 0.0}; }
    static ScoringWeights combined()      { return {0.4, 0.4, 0.2}; }
};

// RESOLUTION STRATEGY

struct ResolutionStrategy {
    bool skipAmbiguous = true;
    ScoringWeights weights;

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

// SCAFFOLD ENTRY

struct ScaffoldEntry {
    static constexpr int UNKNOWN_GAP    = -1;
    static constexpr int DIRECT_OVERLAP =  0;

    size_t contigIndex;          // Index into the contigs vector
    int    gapAfter;             // UNKNOWN_GAP (-1) or DIRECT_OVERLAP (0)
    double score = 0.0;          // Score assigned by scoreContig() - 0.0 if unscored
};

// SCAFFOLD


struct Scaffold {
    std::vector<ScaffoldEntry> entries;
    bool isCircular = false;

    [[nodiscard]] size_t contigCount() const { return entries.size(); }
};

// SCAFFOLDER CLASS

class Scaffolder {

private:

    const std::vector<ContigAssembler::Contig>& contigs_;
    const DeBruijnGraph& graph_;
    const KmerTable*     kmerTable_;
    size_t               maxContigLength_ = 0;
    ResolutionStrategy   strategy_;
    std::vector<Scaffold> scaffolds_;

    // Connection maps: node → list of contig indices starting/ending there
    HashTable<NodeId, std::vector<size_t>> endNodeMap_;
    HashTable<NodeId, std::vector<size_t>> startNodeMap_;

    void buildConnectionMap();

    [[nodiscard]] double computeLengthScore(const ContigAssembler::Contig& contig) const;
    [[nodiscard]] double computeFrequencyScore(const ContigAssembler::Contig& contig) const;
    [[nodiscard]] double computeOverlapScore(const ContigAssembler::Contig& contig) const;

    [[nodiscard]] double scoreContig(size_t contigIndex) const;

    [[nodiscard]] size_t resolveNext(NodeId boundaryNode,
                                     const std::vector<bool>& visited) const;

    [[nodiscard]] Scaffold walkScaffold(size_t startIndex,
                                        std::vector<bool>& visited) const;

    [[nodiscard]] bool isScaffoldStart(size_t contigIndex) const;

public:

    /**
     * @brief Constructor for ContigScaffolder.
     *
     * @param contigs    Contigs produced by ContigAssembler.
     * @param graph      DeBruijnGraph used to build the contigs.
     * @param strategy   Branch resolution strategy (default: skip).
     * @param kmerTable  Optional KmerTable for frequency scoring. nullptr = disabled.
     */
    Scaffolder(const std::vector<ContigAssembler::Contig>& contigs,
                     const DeBruijnGraph& graph,
                     ResolutionStrategy strategy  = ResolutionStrategy::skip(),
                     const KmerTable*   kmerTable = nullptr);

    /**
     * @brief Builds all scaffolds from the contig connection map.
     *
     * Phase 1: walks from all linear entry points.
     * Phase 2: handles remaining unvisited contigs (isolated/circular).
     * Must be called before getScaffolds() or toVisSession().
     */
    void buildScaffolds();

    /**
     * @brief Returns the list of scaffolds.
     * @return Const reference to the scaffold vector.
     */
    [[nodiscard]] const std::vector<Scaffold>& getScaffolds() const;

    /**
     * @brief Converts scaffolding results into VisContig and VisScaffold
     *        structs and populates them into an existing VisSession.
     *
     * Called by the pipeline after buildScaffolds() when visualization
     * is requested. Populates session.contigs and session.scaffolds.
     *
     * @param session  VisSession to populate. Must already have session.k set.
     * @param nodeLen  k-1 - used to decode start/end node labels.
     */
    void toVisSession(VisSession& session, size_t nodeLen) const;

    /**
     * @brief Reports scaffolding statistics to for debugging.
     */
    void printStats() const;
};

#endif //SCAFFOLDER_H