/*
 * de_bruijn_graph.h
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

#include "assembler/core/kmer_types.h"
#include "assembler/core/open_addressing_table.h"

class DeBruijnGraph {

public:

    /**
     * @brief A single outgoing edge, compacted by multiplicity.
     *
     * A de Bruijn graph node (a k-1 mer) has at most 4 distinct outgoing
     * edges - one per possible next base (A/C/G/T). Prior to this struct,
     * addKmer() pushed one raw vector entry per k-mer OCCURRENCE, so a
     * node touched by high-coverage reads could accumulate millions of
     * duplicate entries for the same handful of neighbors. Edge collapses
     * repeat occurrences of the same (from -> to) transition into a single
     * entry with a running count, bounding every node's neighbor list to
     * at most 4 entries regardless of coverage depth.
     *
     * getOutDegree()/getInDegree() still sum to the same multiplicity-aware
     * totals as before (needed for EulerianAssembler's balance check), and
     * traversal (ContigAssembler/EulerianAssembler) consumes one unit of
     * weight per step instead of popping a duplicate - so results are
     * unchanged, only the storage is smaller.
     */
    struct Edge {
        NodeId to;
        size_t weight = 0;
    };

private:

    const size_t k_;
    const __uint128_t kMask_;

    size_t nodeCount_ = 0;
    size_t edgeCount_ = 0;

    struct NodeData {
        std::vector<Edge> neighbors_;      // Outgoing edges only - at most 4 entries.
        size_t inDegree_ = 0;              // Total incoming edge occurrences (multiplicity-aware).
        size_t uniqueInDegree_ = 0;        // Count of distinct predecessor nodes (topology-aware).

        [[nodiscard]] const std::vector<Edge>& getNeighbors() const { return neighbors_; }

        // Multiplicity-aware degree: total edge OCCURRENCES (coverage-weighted).
        // Required by EulerianAssembler, whose balance check and per-occurrence
        // traversal are only correct when they account for read coverage depth.
        [[nodiscard]] size_t getInDegree() const { return inDegree_; }
        [[nodiscard]] size_t getOutDegree() const {
            size_t total = 0;
            for (const auto& e : neighbors_) total += e.weight;
            return total;
        }

        // Topology-aware degree: count of DISTINCT neighbor nodes, ignoring
        // how many times each edge was observed. This is what "is this node
        // a branch point" should mean for unitig extraction (ContigAssembler) -
        // a node covered 100,000x by a single non-branching path is still a
        // single non-branching path, not a branch.
        [[nodiscard]] size_t getUniqueOutDegree() const { return neighbors_.size(); }
        [[nodiscard]] size_t getUniqueInDegree() const { return uniqueInDegree_; }

        std::vector<Edge>& getNeighbors() { return neighbors_; }

        // Increments the existing edge's weight if this transition has been
        // seen before; otherwise adds a new (bounded to <=4) entry. Returns
        // true iff this was a newly created distinct edge, so the caller can
        // update the target's uniqueInDegree_ accordingly.
        bool addNeighbor(NodeId neighbor) {
            for (auto& e : neighbors_) {
                if (e.to == neighbor) { ++e.weight; return false; }
            }
            neighbors_.push_back({neighbor, 1});
            return true;
        }
        void incrementInDegree() { ++inDegree_; }
        void incrementUniqueInDegree() { ++uniqueInDegree_; }
    };

    HashTable<NodeId, NodeData> table_;

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