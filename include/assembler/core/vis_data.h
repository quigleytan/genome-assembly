/*
 * vis_data.h
 * Summary:
 * - Shared data structures used by both the assembly pipeline (writing)
 *   and the visualizer (reading). No external dependencies beyond STL
 *   and kmer_types.h, so both sides can include this without issues.
 * - NodeId raw values are stored for identity comparisons, but all
 *   human-readable labels are pre-decoded strings, so the Visualizer does
 *   not need to include kmer_encoding.
 * Important notes:
 * - VisSession is the top-level container passed to data_exporter/data_loader.
 * - contigSteps drives the contig assembly animation in contig_view.
 */

#ifndef VIS_DATA_H
#define VIS_DATA_H

#include <vector>
#include <string>
#include <cstdint>

#include "assembler/core/kmer_types.h"

// ANIMATION STEPS

struct TraversalStep {
    enum class Type {
        // Contig animation events
        ContigStarted,    // A new walkContig() call has begun
        BaseAppended,     // A run of consecutive bases was added to the current contig
        ContigFinished,   // walkContig() returned: contig is complete

    };

    Type   type;
    NodeId from = 0;        // Source node (edge steps)
    NodeId to   = 0;        // Destination node (edge steps)

    // Contig step fields
    size_t contigIndex = 0; // Which contig this step belongs to
    char   base        = 0; // The base appended (BaseAppended steps only)
    size_t count       = 0; // Number of bases this run covers (BaseAppended only).
                            // Recorder coalesces consecutive same-contig
                            // base-append calls into one step with a running
                            // count instead of one step per base, so replaying
                            // a contig's growth is O(contigs) not O(bases).
    std::string sequence;   // Full sequence snapshot (ContigFinished only,
                            // avoids storing partial strings at every step)
};

// CONTIG SCAFFOLDING STRUCTS

struct VisContig {
    std::string sequence;           // Full assembled base sequence
    std::string startLabel;         // Decoded k-1 mer label for startNode
    std::string endLabel;           // Decoded k-1 mer label for endNode
    NodeId      startNode = 0;      // Raw encoded value (for identity comparisons)
    NodeId      endNode   = 0;
    bool        isCircular = false;
    double      score      = 0.0;   // Scaffolder score (0.0 if not scored)
    int         scaffoldIndex = -1; // Which scaffold this contig belongs to (-1 = unscaffolded)
};

// SCAFFOLD REPRESENTATION

struct VisScaffold {
    std::vector<size_t> contigIndices; // Indices into VisSession::contigs
    std::vector<int>    gaps;          // Gap after each contig (-1 = unknown, 0 = direct overlap)
    bool isCircular = false;
};

// GAP RESOLUTION
//
// One entry per junction between consecutive scaffolds in genomeSequence's
// assembly order, i.e. VisSession::gapFills.size() == scaffolds.size() - 1
// (0 if there are 0 or 1 scaffolds). gapFills[i] describes the junction
// between scaffolds[i] and scaffolds[i+1] - produced by GapEstimator/
// GapFiller, not by Scaffolder itself.

struct VisGapFill {
    bool   resolved     = false; // True if bridged via local reassembly over the graph.
    size_t estimatedGap = 0;     // GapEstimator's estimate, in bases.
    size_t filledLength = 0;     // Bases actually inserted into genomeSequence at this
                                  // junction - the reassembled bridge length when resolved,
                                  // or the N-run length (== estimatedGap) otherwise.
};

// CONTAINER

struct VisSession {
    // Assembly metadata
    size_t      k = 0;
    std::string sourceFile;     // Path to the FASTA/FASTQ file used
    std::string strategyName;   // "skip", "greedy", or "scored"
    std::string genomeSequence; // Reconstructed genome sequence

    // Contig and scaffold data
    std::vector<VisContig>   contigs;
    std::vector<VisScaffold> scaffolds;
    std::vector<VisGapFill>  gapFills;   // Inter-scaffold junctions, see VisGapFill above.

    // Animation step lists
    std::vector<TraversalStep> contigSteps;  // Phase 1 animation
};

#endif // VIS_DATA_H