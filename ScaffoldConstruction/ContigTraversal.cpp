#include "ContigTraversal.h"

#include <algorithm>
#include <iostream>

#include "DataProcessing/KmerEncoding.h"

// PRIVATE HELPER FUNCTIONS

void ContigTraversal::initializeAdjacency() {
    adjCopy_ = OpenAddressingTable<NodeId, std::vector<NodeId>>(
        graph_.getNodeCount() * 2);

    // Pass 1: Insert all nodes
    for (NodeId node : graph_.getAllNodes()) {
        adjCopy_.insert(node);
    }

    // Pass 2: Assign and sort neighbor lists
    for (NodeId node : graph_.getAllNodes()) {
        const auto* data = graph_.findNode(node);
        auto* neighborRef = adjCopy_.find(node);
        *neighborRef = data->getNeighbors();
        std::sort(neighborRef->begin(), neighborRef->end());
    }
}

bool ContigTraversal::isAmbiguous(NodeId node) const {
    const auto* data = graph_.findNode(node);
    return data->getInDegree() > 1 || data->getOutDegree() > 1;
}

ContigTraversal::Contig ContigTraversal::walkContig(NodeId startNode) {
    Contig result;
    result.startNode = startNode;
    result.isCircular = false;

    const size_t nodeLen = graph_.getK() - 1;
    result.sequence = KmerEncoding::decode(startNode, nodeLen);

    NodeId current = startNode;

    while (true) {
        auto* neighbors = adjCopy_.find(current);
        if (!neighbors || neighbors->empty()) {
            result.endNode = current;
            break;
        }

        NodeId next = neighbors->back();
        neighbors->pop_back();
        result.sequence += KmerEncoding::decode(next, nodeLen).back();

        if (next == startNode) {
            result.endNode = startNode;
            result.isCircular = true;
            size_t overlap = nodeLen - 1; // k-2
            if (result.sequence.length() > overlap)
                result.sequence.resize(result.sequence.length() - overlap);
            break;
        }

        const auto* nextData = graph_.findNode(next);
        if (nextData->getInDegree() != 1 || nextData->getOutDegree() != 1) {
            result.endNode = next;
            break;
        }
        current = next;
    }
    return result;
}

void ContigTraversal::handleIsolatedCycles() {
    for (NodeId node : graph_.getAllNodes()) {
        // Checks for ambiguity - skips if ambiguous.
        if (isAmbiguous(node)) continue;

        auto* neighbors = adjCopy_.find(node); // Retrieves the neighbor list for this node.

        while (neighbors && !neighbors->empty()) {
            Contig contig = walkContig(node);

            if (!contig.sequence.empty())
                contigs_.push_back(std::move(contig));
        }
    }
}

// PUBLIC

ContigTraversal::ContigTraversal(DeBruijnGraph& g, Recorder* recorder)
    : graph_(g), adjCopy_(g.getNodeCount() * 2) {}

void ContigTraversal::computeContigs() {
    // Preparing data structures.
    contigs_.clear();
    initializeAdjacency();

    // Phase 1 - General contig walks.
    for (NodeId node : graph_.getAllNodes()) {
        const auto* data = graph_.findNode(node);

        bool isBranchPoint = data->getInDegree() > 1 || data->getOutDegree() > 1;
        bool isSource      = data->getInDegree() == 0 && data->getOutDegree() >= 1;

        if (!isBranchPoint && !isSource) continue; // Skips nodes that are not unitig boundaries.

        auto* neighbors = adjCopy_.find(node);

        while (neighbors && !neighbors->empty()) {
            Contig contig = walkContig(node);

            if (!contig.sequence.empty()) {
                contigs_.push_back(std::move(contig));
            }
        }
    }

    // Phase 2 - Closed loop handling.
    handleIsolatedCycles();
}

const std::vector<ContigTraversal::Contig>& ContigTraversal::getContigs() const {
    return contigs_;
}

void ContigTraversal::printStats() const {
    if (contigs_.empty()) {
        std::cout << "No contigs found\n";
        return;
    }

    // Initializing reporting variables.
    size_t totalLength = 0;
    size_t circularCount = 0;
    std::vector<size_t> lengths;
    lengths.reserve(contigs_.size());

    // Length calculation.
    for (const auto& contig : contigs_) {
        size_t len = contig.sequence.length();
        lengths.push_back(len);
        totalLength += len;
        if (contig.isCircular) ++circularCount;
    }

    // Sorting lengths vector for N50 calculation.
    std::sort(lengths.rbegin(), lengths.rend());

    size_t half = (totalLength + 1) / 2;
    size_t accumulated = 0;
    size_t n50 = 0;
    for (size_t len : lengths) {
        accumulated += len;
        if (accumulated >= half) { n50 = len; break; }
    }

    std::cout << "--------------------------------------\n";
    std::cout << "Total contigs:    " << contigs_.size() << "\n";
    std::cout << "Circular contigs: " << circularCount << "\n";
    std::cout << "Total bases:      " << totalLength << "\n";
    std::cout << "Shortest contig:  " << lengths.back() << " bases\n";
    std::cout << "Longest contig:   " << lengths.front() << " bases\n";
    std::cout << "N50:              " << n50 << " bases\n";
    std::cout << "--------------------------------------\n";
}
