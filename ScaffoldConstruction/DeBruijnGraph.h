/*
 * DeBruijnGraph.h
 * Created by Tanner Quigley on 2/15/2026
 * Summary:
 * - Representation of a De Bruijn graph to be used for Eulerian walks for genome assembly.
 * - Graph nodes are k-1 mers, transitions are represented as full kmers.
 * Important notes:
 * - IMPORTANT: All kmer and node lookups should already be encoded as an __uint128_t.
 * - Graph stored in a linear probing, open addressing hash table, with an encoded k-1mer (NodeId) as the key and
 *   the struct NodeData in value.
 */

#ifndef DE_BRUIJN_GRAPH_H
#define DE_BRUIJN_GRAPH_H

#include <cstdint>
#include <iostream>
#include <vector>
#include "KmerTypes.h"
#include "DataProcessing/OpenAddressingTable.h"

class DeBruijnGraph {

private:

    const size_t k_;
    const __uint128_t kMask_;

    size_t nodeCount_ = 0;
    size_t edgeCount_ = 0;

    // Note: outdegree is tracked by neighbors.size
    struct NodeData {
        std::vector<NodeId> neighbors_; // Outgoing edges only
        size_t inDegree_ = 0;             // Number of edges entering the node

        [[nodiscard]] const std::vector<NodeId>& getNeighbors() const { return neighbors_; }
        [[nodiscard]] size_t getInDegree() const { return inDegree_; }
        [[nodiscard]] size_t getOutDegree() const { return neighbors_.size(); }

        std::vector<NodeId>& getNeighbors() { return neighbors_; }

        void addNeighbor(NodeId neighbor) { neighbors_.push_back(neighbor); }
        void incrementInDegree() { ++inDegree_; }
    };

    OpenAddressingTable<NodeId, NodeData> table_;

    /**
     * @brief Extracts k-1 prefix and suffix to create nodes.
     *
     * "Chops" the kmer into its prefix and suffix, allowing for storage in the table.
     * Applies a bitmask to keep the size of nodes at k-1.
     *
     * @param kmer The kmer in which the prefix and suffix will be pulled from.
     * @return Returns a tuple containing k-1 prefix and suffix.
     */
    [[nodiscard]] std::pair<NodeId, NodeId> chop(NodeId kmer) const;

public:

    /**
     * @brief Constructor for De Bruijn graph class.
     *
     * Initializes a De Bruijn graph for k-1 mers using size k given upon construction.
     * Contains a hash table of nodes containing neighbor and out/indegree information.
     * All kmers to be added to the graph MUST BE ENCODED prior to insertion.
     *
     * @param k The k size for the graph: k must be > 1.
     * @param initialSize Default at 101, used to pre-size the table to reduce rehash collisions.
     */
    explicit DeBruijnGraph(size_t k, size_t initialSize = 101);

    /**
     * @brief Returns the size of K-mers in the graph.
     * @return k_ The k size for the graph.
     */
    [[nodiscard]] size_t getK() const;

    /**
     * @brief Returns the total node counts in the graph.
     * @return Number of nodes (k-1 mers).
     */
    [[nodiscard]] size_t getNodeCount() const;

    /**
     * @brief Returns the total edge counts in the graph.
     * @return Number of edges (kmer transitions).
     */
    [[nodiscard]] size_t getEdgeCount() const;

    /**
     * @brief Checks table to the desired node.
     * @param queriedNode Encoded node to be found in the table.
     * @return Returns a pointer to the node requested, nullptr if not found.
     */
    [[nodiscard]] const NodeData *findNode(NodeId queriedNode) const;

    /**
     * @brief Adds the requested kmer to the graph.
     *
     * Adds both the prefix and suffix of inputKmer to the table and adds the suffix to
     * the prefix's list of neighbors (individual adjacency list).
     *
     * @param inputKmer The kmer whose derived k-1mers will be added to the table.
     */
    void addKmer(NodeId inputKmer);

    /**
     * @brief Returns set of all k-1 nodes.
     *
     * Iterates through the graph and adds each node's key (encoded base as an uint_64) found
     * to a list of nodes.
     *
     * @return Vector of all nodes in the graph (size k-1).
     */
    [[nodiscard]] std::vector<NodeId> getAllNodes() const;

    void printGraph() const;
};

#endif //DE_BRUIJN_GRAPH_H