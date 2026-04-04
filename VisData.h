/*
 * VisData.h
 * Summary:
 * - Shared data structures used by both the assembly pipeline (writing)
 *   and the Visualizer (reading). No external dependencies beyond STL
 *   and KmerTypes.h, so both sides can include this without issues.
 * - NodeId raw values are stored for identity comparisons, but all
 *   human-readable labels are pre-decoded strings so the Visualizer
 *   does not need to link KmerEncoding.
 * Important notes:
 * - VisSession is the top-level container passed to DataExporter/DataLoader.
 * - contigSteps drives the contig assembly animation in ContigView.
 * - eulerSteps is reserved for phase 2 (GraphView / Eulerian animation).
 */

#ifndef VIS_DATA_H
#define VIS_DATA_H

#include <vector>
#include <string>
#include <cstdint>
#include "KmerTypes.h"

// Animation steps

struct TraversalStep {
    enum class Type {
        // Contig animation events
        ContigStarted,    // A new walkContig() call has begun
        BaseAppended,     // One base was added to the current contig sequence
        ContigFinished,   // walkContig() returned; contig is complete

        // Eulerian animation events (reserved, unused in phase 1)
        EdgeConsumed,     // Hierholzer consumed one edge
        NodeCommitted,    // A node was committed to the Eulerian path
    };

    Type   type;
    NodeId from = 0;        // Source node (edge steps)
    NodeId to   = 0;        // Destination node (edge steps)

    // Contig step fields
    size_t contigIndex = 0; // Which contig this step belongs to
    char   base        = 0; // The base appended (BaseAppended steps only)
    std::string sequence;   // Full sequence snapshot (ContigFinished only,
                            // avoids storing partial strings at every step)
};

// CONTIG SCAFFOLDING STRUCTS

struct VisContig {
    std::string sequence;       // Full assembled base sequence
    std::string startLabel;     // Decoded k-1 mer label for startNode
    std::string endLabel;       // Decoded k-1 mer label for endNode
    NodeId      startNode = 0;  // Raw encoded value (for identity comparisons)
    NodeId      endNode   = 0;
    bool        isCircular = false;
    double      score      = 0.0;  // Scaffolder score (0.0 if not scored)
    int         scaffoldIndex = -1; // Which scaffold this contig belongs to (-1 = unscaffolded)
};

// Scaffold representation

struct VisScaffold {
    std::vector<size_t> contigIndices; // Indices into VisSession::contigs
    std::vector<int>    gaps;          // Gap after each contig (-1 = unknown, 0 = direct overlap)
    bool isCircular = false;
};

// EULERIAN TRAVERSAL STRUCTS

struct VisNode {
    NodeId      id = 0;
    std::string label;       // Decoded k-1 mer string
    float       x = 0.0f;
    float       y = 0.0f;
    size_t      inDegree  = 0;
    size_t      outDegree = 0;
};

struct VisEdge {
    NodeId from         = 0;
    NodeId to           = 0;
    size_t multiplicity = 1; // How many times this edge appears in the graph
};

// Container

struct VisSession {
    // Assembly metadata
    size_t      k = 0;
    std::string sourceFile;   // Path to the FASTA/FASTQ file used
    std::string strategyName; // "skip", "greedy", "scored", or "eulerian"
    std::string genomeSequence; // Reconstructed genome sequence

    // Contig and scaffold data (phase 1)
    std::vector<VisContig>   contigs;
    std::vector<VisScaffold> scaffolds;

    // Graph data (phase 2 — left empty in phase 1 runs)
    std::vector<VisNode> nodes;
    std::vector<VisEdge> edges;

    // Animation step lists
    std::vector<TraversalStep> contigSteps;  // Phase 1 animation
    std::vector<TraversalStep> eulerSteps;   // Phase 2 animation (reserved)
};

#endif // VIS_DATA_H