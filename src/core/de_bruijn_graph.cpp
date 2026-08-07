#include "assembler/core/de_bruijn_graph.h"
#include "assembler/construction/kmer_encoding.h"

// PRIVATE HELPER FUNCTION

std::pair<NodeId, NodeId> DeBruijnGraph::chop(NodeId kmer) const {
    __uint128_t prefix = kmer >> 2; // Naturally discards
    __uint128_t suffix = kmer & kMask_;

    return {prefix, suffix};
}

// PUBLIC

DeBruijnGraph::DeBruijnGraph(size_t k, size_t expectedNodes)
    : k_(KmerEncoding::validateK(k)),
      kMask_(KmerEncoding::bitmask(k_ - 1)),
      table_(expectedNodes * 2) {}

// GETTERS

size_t DeBruijnGraph::getK() const {
    return k_;
}

size_t DeBruijnGraph::getNodeCount() const {
    return nodeCount_;
}

size_t DeBruijnGraph::getEdgeCount() const {
    return edgeCount_;
}

// GRAPH FUNCTIONALITY

const DeBruijnGraph::NodeData *DeBruijnGraph::findNode(NodeId node) const{
    return table_.find(node);
}

void DeBruijnGraph::addKmer(NodeId inputKmer) {
    auto [prefix, suffix] = chop(inputKmer);

    auto prefixNew = table_.insert(prefix).second;
    auto suffixNew = table_.insert(suffix).second;

    if (prefixNew) nodeCount_++;
    if (suffixNew) nodeCount_++;

    NodeData* from = table_.find(prefix);
    NodeData* to   = table_.find(suffix);

    bool newEdge = from->addNeighbor(suffix);
    to->incrementInDegree();
    if (newEdge) to->incrementUniqueInDegree();
    edgeCount_++;
}

std::vector<NodeId> DeBruijnGraph::getAllNodes() const {
    // Initializes the return vector.
    std::vector<NodeId> nodes;
    nodes.reserve(nodeCount_);
    // Populates the vector with all the item keys.
    for (auto it = table_.begin(); it != table_.end(); ++it) {
        nodes.push_back((*it).key);
    }
    return nodes;
}

void DeBruijnGraph::printGraph() const {
    std::cout << "De Bruijn Graph\n";
    std::cout << "--------------------------------------\n";

    auto nodes = getAllNodes();

    for (NodeId node : nodes) {
        std::cout << KmerEncoding::decode(node, k_ - 1)
                  << " | in: " << findNode(node)->getInDegree()
                  << " | out: " << findNode(node)->getOutDegree()
                  << " | -> ";
        auto neighbors = findNode(node)->getNeighbors();
        for (const auto& edge : neighbors) {
            std::cout << KmerEncoding::decode(edge.to, k_ - 1)
                      << "(x" << edge.weight << ") ";
        }
        std::cout << "\n";
    }
    std::cout << "--------------------------------------\n";
}
