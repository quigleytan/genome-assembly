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
 * TODO: Build in handling for starting traversal at "sinks".
 */

#ifndef CONTIG_TRAVERSAL_H
#define CONTIG_TRAVERSAL_H

#include <vector>
#include <string>

#include "KmerTypes.h"
#include "Recorder.h"
#include "ScaffoldConstruction/DeBruijnGraph.h"
#include "DataProcessing/OpenAddressingTable.h"

class ContigTraversal {

public: // PUBLIC CLASS MEMBERS

    struct Contig {
        std::string sequence;    // Assembled base sequence including k-1 overlap at boundaries.
        NodeId startNode;        // Encoded k-1 mer where the walk began.
        NodeId endNode;          // Encoded k-1 mer where the walk terminated.
        bool isCircular = false; // True if the walk returned to startNode forming a closed loop.
    };

private:

    DeBruijnGraph& graph_;                                     // Graph to be traversed.
    OpenAddressingTable<NodeId, std::vector<NodeId>> adjCopy_; // Table of all nodes and neighbor lists.
    std::vector<Contig> contigs_;                              // List of all contigs derived from graph_.

    /**
     * @brief Initializes adjCopy_, table of ID's with neighbor lists as values.
     * Pre-allocates adjCopy_ size to avoid rehashing.
     */
    void initializeAdjacency();

    /**
     * @brief Boolean test for node ambiguity
     * @param node Node to be tested.
     * @return bool True if the node has inDegree > 1 or outDegree > 1, false otherwise.
     */
    [[nodiscard]] bool isAmbiguous(NodeId node) const;

    /**
     * @brief Walks the graph from startNode and assembles a single contig.
     *
     * Emits the full k-1 char seed of startNode, then appends the last character
     * of each subsequent node until a unitig boundary or dead end is reached.
     * Detects circular contigs if the walk returns to startNode, trimming the
     * k-2 overlap introduced by the circular join.
     *
     * @param startNode The encoded k-1 mer to begin the walk from.
     * @return Contig struct containing the assembled sequence, start/end nodes,
     *         and circularity flag.
     */
    [[nodiscard]] Contig walkContig(NodeId startNode);

    /**
     * @brief Processes unreached cycles from phase 1 of contig construction.
     *
     * Stage 2 of contig construction: occurs after the main traversal processes all ambiguous and source nodes.
     * Isolated cycles refer to any unambiguous with remaining edges that were inaccessible due to no external
     * entry/exit points. Completes contig walks for each node found.
     */
    void handleIsolatedCycles();

public:

    /**
     * @brief Constructor for Contig Traversal class.
     *
     * Initializes adjCopy_ to 2 * the number of nodes in the argument graph to avoid extra allocation later.
     *
     * @param g Input graph to be traversed.
     */
    explicit ContigTraversal(DeBruijnGraph& g, Recorder* recorder = nullptr);

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
     * @return List of contig structs.
     */
    [[nodiscard]] const std::vector<Contig>& getContigs() const;

    /**
     * @brief Reports the information about the contigs.
     *
     * Total contigs:    Number of contiguous fragments derived.
     * Circular contigs: Number of contigs that form a closed loop.
     * Total bases:      Sum of the length of each contig.
     * N50:              Measure of contiguity, the length at which 50% of total assembled bases are contained in
     *                   contigs of said length or longer.
     */
    void printStats() const;
};

#endif