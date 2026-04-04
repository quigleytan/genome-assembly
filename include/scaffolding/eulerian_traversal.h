/*
 * EulerianTraversal.h
 * Created by Tanner Quigley on 2/15/2026
 * Summary:
 * - Computes an Eulerian path or circuit through a DeBruijnGraph using
 *   Hierholzer's algorithm, visiting every edge exactly once.
 * - Reconstructs a DNA sequence string from the resulting node path.
 * Important notes:
 * - Requires a valid Eulerian graph — exactly one node with outDegree - inDegree == 1
 *   for a path, or all nodes balanced for a circuit.
 * - computePath() must be called before getPath() or reconstructGenome().
 */

#ifndef EULERIAN_PATH_H
#define EULERIAN_PATH_H

#include <vector>
#include <string>

#include "kmer_types.h"
#include "scaffolding/de_bruijn_graph.h"
#include "encoding/open_addressing_table.h"


class EulerianTraversal {

private:

    DeBruijnGraph& graph_;                                     // Graph to be traversed.
    HashTable<NodeId, std::vector<NodeId>> adjCopy_; // Table of all nodes and neighbor lists.
    std::vector<NodeId> path_;                                 // List of nodes visited.

    /**
     * @brief Initializes adjCopy_, table of ID's with neighbor lists as values.
     * Sorting the neighbor lists ensures a deterministic traversal order.
     */
    void initializeAdjacency();

    /**
     * @brief Identifies the correct start node for the Eulerian traversal.
     *
     * For a path, returns the unique node where outDegree - inDegree == 1.
     * For a circuit, all nodes are balanced so returns an arbitrary node.
     *
     * @return NodeId of the start node.
     * @throws std::runtime_error if the graph is not Eulerian.
     */
    [[nodiscard]] NodeId findStartNode() const;

public:

    /**
     * @brief Constructor for EulerianTraversal.
     * @param g The DeBruijnGraph to traverse.
     */
    EulerianTraversal(DeBruijnGraph& g);

    /**
     * @brief Computes an Eulerian path or circuit using Hierholzer's algorithm.
     *
     * Traverses the De Bruijn graph visiting every edge exactly once. Begins
     * from the start node determined by findStartNode() — a node with
     * outDegree - inDegree == 1 for a path, or any node for a circuit.
     *
     * At each step, if the current node has remaining edges, one edge is
     * consumed and the next node is pushed onto the stack. If no edges remain,
     * the current node is committed to the path and the algorithm backtracks.
     * The path is reversed at the end since nodes are committed in reverse
     * order during backtracking.
     *
     * Results are stored in path_ and can be retrieved via getPath() or
     * converted to a DNA string via reconstructGenome().
     *
     * @throws std::runtime_error if the graph does not contain an Eulerian
     *         path or circuit.
     */
    void computePath();

    /**
     * @brief Returns the path (reconstructed sequence).
     * @return Const reference to the vector of node IDs in traversal order.
     */
    [[nodiscard]] const std::vector<NodeId>& getPath() const;

    /**
     * @brief Reconstructs a DNA string from the computed Eulerian path.
     *
     * Decodes the sequence of k-1 mer node IDs stored in path_ back into a
     * contiguous DNA string. The first node contributes its full k-1 chars as
     * the seed, and each subsequent node contributes only its last character
     * to avoid double-counting the k-2 base overlap between adjacent nodes.
     *
     * For circuits, the final node is excluded from the loop since it is
     * identical to the first node and its character is already represented
     * in the seed.
     *
     * @param isCircular True if the path is a circuit (path.front == path.back),
     *                  false if it is a linear Eulerian path.
     * @return Reconstructed DNA sequence string.
     * @throws std::runtime_error if called before computePath().
     */
    [[nodiscard]] std::string reconstructGenome(bool isCircular) const;
};

#endif