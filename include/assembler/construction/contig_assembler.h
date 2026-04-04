/*
* ContigTraversal.h
 * Created by Tanner Quigley on 3/18/2026
 * Summary:
 * - Traverses a DeBruijnGraph to extract contiguous DNA fragments (contigs).
 * - Contigs are produced in two phases: branch points and source nodes first,
 *   then isolated cycles unreachable from external entry points.
 * - All contigs include k-1 overlapping bases at boundaries for scaffolding.
 * Important notes:
 * - Operates on an already-constructed DeBruijnGraph passed by reference.
 * - Contig sequences are produced in overlap mode — boundary nodes contribute
 *   their full k-1 chars to adjacent contigs.
 * - recorder_ is optional (default nullptr). When null, all recording calls
 *   are no-ops — zero overhead for normal non-visualization pipeline runs.
 * TODO: Build in handling for starting traversal at "sinks".
 */

#ifndef CONTIG_TRAVERSAL_H
#define CONTIG_TRAVERSAL_H

#include <vector>
#include <string>

#include "assembler/core/recorder.h"
#include "assembler/core/de_bruijn_graph.h"
#include "assembler/core/open_addressing_table.h"

class ContigAssembler {

public:

    struct Contig {
        std::string sequence;    // Assembled base sequence including k-1 overlap at boundaries.
        NodeId startNode;        // Encoded k-1 mer where the walk began.
        NodeId endNode;          // Encoded k-1 mer where the walk terminated.
        bool isCircular = false; // True if the walk returned to startNode forming a closed loop.
    };

private:

    DeBruijnGraph& graph_;                                     // Graph to be traversed.
    HashTable<NodeId, std::vector<NodeId>> adjCopy_; // Table of all nodes and neighbor lists.
    std::vector<Contig> contigs_;                              // List of all contigs derived from graph_.
    Recorder* recorder_;                                       // Optional — null in normal pipeline runs.

    /**
     * @brief Initializes adjCopy_, table of ID's with neighbor lists as values.
     * Pre-allocates adjCopy_ size to avoid rehashing.
     */
    void initializeAdjacency();

    /**
     * @brief Boolean test for node ambiguity.
     * @param node Node to be tested.
     * @return True if the node has inDegree > 1 or outDegree > 1.
     */
    [[nodiscard]] bool isAmbiguous(NodeId node) const;

    /**
     * @brief Walks the graph from startNode and assembles a single contig.
     *
     * Emits the full k-1 char seed of startNode, then appends the last character
     * of each subsequent node until a unitig boundary or dead end is reached.
     * Detects circular contigs if the walk returns to startNode, trimming the
     * k-2 overlap introduced by the circular join.
     * If recorder_ is non-null, emits ContigStarted, BaseAppended, and
     * ContigFinished steps for the animation.
     *
     * @param startNode   The encoded k-1 mer to begin the walk from.
     * @param contigIndex Index this contig will occupy in contigs_ — passed
     *                    to the recorder so steps reference the right contig.
     * @return Contig struct containing the assembled sequence, start/end nodes,
     *         and circularity flag.
     */
    [[nodiscard]] Contig walkContig(NodeId startNode, size_t contigIndex);

    /**
     * @brief Processes unreached cycles from phase 1 of contig construction.
     *
     * Stage 2 of contig construction: occurs after the main traversal processes
     * all ambiguous and source nodes. Isolated cycles refer to any unambiguous
     * nodes with remaining edges that were inaccessible due to no external
     * entry/exit points.
     */
    void handleIsolatedCycles();

public:

    /**
     * @brief Constructor for ContigTraversal.
     *
     * @param g        Input graph to be traversed.
     * @param recorder Optional recorder for visualization. Pass nullptr (default)
     *                 for normal pipeline runs — all recording calls become no-ops.
     */
    explicit ContigAssembler(DeBruijnGraph& g, Recorder* recorder = nullptr);

    /**
     * @brief Constructs all contigs from the graph in two phases.
     *
     * Phase 1: walks from all branch points and source nodes, consuming their
     * outgoing edges and recording each resulting contig.
     * Phase 2: calls handleIsolatedCycles() to walk any remaining isolated cycles
     * unreachable from phase 1.
     */
    void computeContigs();

    /**
     * @brief Returns the list of contigs for graph_.
     * @return Const reference to the contig list.
     */
    [[nodiscard]] const std::vector<Contig>& getContigs() const;

    /**
     * @brief Reports contig statistics to stdout.
     *
     * Total contigs:    Number of contiguous fragments derived.
     * Circular contigs: Number of contigs that form a closed loop.
     * Total bases:      Sum of the length of each contig.
     * N50:              Length at which 50% of total assembled bases are
     *                   contained in contigs of that length or longer.
     */
    void printStats() const;
};

#endif