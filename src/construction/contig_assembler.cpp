#include "assembler/construction/contig_assembler.h"

#include <algorithm>
#include <iostream>

#include "assembler/construction/kmer_encoding.h"

// PRIVATE HELPER FUNCTIONS

void ContigAssembler::initializeAdjacency() {
    adjCopy_ = HashTable<NodeId, std::vector<DeBruijnGraph::Edge>>(
        graph_.getNodeCount() * 2);

    // Pass 1: Insert all nodes
    for (NodeId node : graph_.getAllNodes()) {
        adjCopy_.insert(node);
    }

    // Pass 2: Assign and sort neighbor lists (at most 4 entries per node,
    // so this copy/sort is now O(nodeCount) instead of O(total occurrences)).
    for (NodeId node : graph_.getAllNodes()) {
        const auto* data = graph_.findNode(node);
        auto* neighborRef = adjCopy_.find(node);
        *neighborRef = data->getNeighbors();
        std::sort(neighborRef->begin(), neighborRef->end(),
                  [](const DeBruijnGraph::Edge& a, const DeBruijnGraph::Edge& b) {
                      return a.to < b.to;
                  });
    }
}

bool ContigAssembler::isAmbiguous(NodeId node) const {
    const auto* data = graph_.findNode(node);
    // Topology-aware: a node covered many times by one non-branching path
    // is still non-branching. Multiplicity is handled separately via each
    // edge's weight during traversal (see walkContig).
    return data->getUniqueInDegree() > 1 || data->getUniqueOutDegree() > 1;
}

ContigAssembler::Contig ContigAssembler::walkContig(NodeId startNode, size_t contigIndex) {
    Contig result;
    result.startNode  = startNode;
    result.isCircular = false;

    const size_t nodeLen = graph_.getK() - 1;
    result.sequence = KmerEncoding::decode(startNode, nodeLen);

    // Record the walk started at startNode with the full k-1 char seed.
    if (recorder_)
        recorder_->contigStarted(contigIndex, startNode,
                                 KmerEncoding::decode(startNode, nodeLen));

    NodeId current = startNode;

    while (true) {
        auto* neighbors = adjCopy_.find(current);
        if (!neighbors || neighbors->empty()) {
            result.endNode = current;

            // Record the walk ended at a dead end.
            if (recorder_)
                recorder_->contigFinished(contigIndex, current,
                                          KmerEncoding::decode(current, nodeLen),
                                          result.sequence, false);
            break;
        }

        // Consume one unit of weight from the last edge; only remove the
        // entry once its full multiplicity has been walked.
        DeBruijnGraph::Edge& edge = neighbors->back();
        NodeId next = edge.to;
        if (--edge.weight == 0) neighbors->pop_back();

        char appendedBase = KmerEncoding::decode(next, nodeLen).back();
        result.sequence += appendedBase;

        // Record one base was appended.
        if (recorder_)
            recorder_->baseAppended(contigIndex, appendedBase, next);

        if (next == startNode) {
            result.endNode    = startNode;
            result.isCircular = true;
            size_t overlap = nodeLen; // k-1 chars to trim
            if (result.sequence.length() > overlap)
                result.sequence.resize(result.sequence.length() - overlap);

            // Record the walk returned to the start.
            if (recorder_)
                recorder_->contigFinished(contigIndex, startNode,
                                          KmerEncoding::decode(startNode, nodeLen),
                                          result.sequence, true);
            break;
        }

        const auto* nextData = graph_.findNode(next);
        if (nextData->getUniqueInDegree() != 1 || nextData->getUniqueOutDegree() != 1) {
            result.endNode = next;

            // Record the walk ended at a branch point.
            if (recorder_)
                recorder_->contigFinished(contigIndex, next,
                                          KmerEncoding::decode(next, nodeLen),
                                          result.sequence, false);
            break;
        }

        current = next;
    }

    return result;
}

void ContigAssembler::handleIsolatedCycles() {
    for (NodeId node : graph_.getAllNodes()) {
        if (isAmbiguous(node)) continue;

        auto* neighbors = adjCopy_.find(node);

        while (neighbors && !neighbors->empty()) {
            size_t contigIndex = contigs_.size(); // index this contig will occupy
            Contig contig = walkContig(node, contigIndex);

            if (!contig.sequence.empty())
                contigs_.push_back(std::move(contig));
        }
    }
}

// PUBLIC

ContigAssembler::ContigAssembler(DeBruijnGraph& g, Recorder* recorder)
    : graph_(g),
      adjCopy_(g.getNodeCount() * 2),
      recorder_(recorder) {}

void ContigAssembler::computeContigs() {
    contigs_.clear();
    initializeAdjacency();

    // Phase 1 - branch points and source nodes
    for (NodeId node : graph_.getAllNodes()) {
        const auto* data = graph_.findNode(node);

        bool isBranchPoint = data->getUniqueInDegree() > 1 || data->getUniqueOutDegree() > 1;
        bool isSource      = data->getUniqueInDegree() == 0 && data->getUniqueOutDegree() >= 1;

        if (!isBranchPoint && !isSource) continue;

        auto* neighbors = adjCopy_.find(node);

        while (neighbors && !neighbors->empty()) {
            size_t contigIndex = contigs_.size(); // index this contig will occupy
            Contig contig = walkContig(node, contigIndex);

            if (!contig.sequence.empty())
                contigs_.push_back(std::move(contig));
        }
    }

    // Stage 2 - isolated cycles
    handleIsolatedCycles();
}

const std::vector<ContigAssembler::Contig>& ContigAssembler::getContigs() const {
    return contigs_;
}

void ContigAssembler::printStats() const {
    if (contigs_.empty()) {
        std::cout << "No contigs found\n";
        return;
    }

    size_t totalLength   = 0;
    size_t circularCount = 0;
    std::vector<size_t> lengths;
    lengths.reserve(contigs_.size());

    for (const auto& contig : contigs_) {
        size_t len = contig.sequence.length();
        lengths.push_back(len);
        totalLength += len;
        if (contig.isCircular) ++circularCount;
    }

    std::sort(lengths.rbegin(), lengths.rend());

    size_t half        = (totalLength + 1) / 2;
    size_t accumulated = 0;
    size_t n50         = 0;
    for (size_t len : lengths) {
        accumulated += len;
        if (accumulated >= half) { n50 = len; break; }
    }

    std::cout << "--------------------------------------\n";
    std::cout << "Total contigs:    " << contigs_.size()  << "\n";
    std::cout << "Circular contigs: " << circularCount    << "\n";
    std::cout << "Total bases:      " << totalLength      << "\n";
    std::cout << "Shortest contig:  " << lengths.back()   << " bases\n";
    std::cout << "Longest contig:   " << lengths.front()  << " bases\n";
    std::cout << "N50:              " << n50              << " bases\n";
    std::cout << "--------------------------------------\n";
}