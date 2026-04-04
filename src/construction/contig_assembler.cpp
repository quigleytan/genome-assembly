#include "assembler/construction/contig_assembler.h"

#include <algorithm>
#include <iostream>

#include "assembler/construction/kmer_encoding.h"

// PRIVATE HELPER FUNCTIONS

void ContigAssembler::initializeAdjacency() {
    adjCopy_ = HashTable<NodeId, std::vector<NodeId>>(
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

bool ContigAssembler::isAmbiguous(NodeId node) const {
    const auto* data = graph_.findNode(node);
    return data->getInDegree() > 1 || data->getOutDegree() > 1;
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

        NodeId next = neighbors->back();
        neighbors->pop_back();

        char appendedBase = KmerEncoding::decode(next, nodeLen).back();
        result.sequence += appendedBase;

        // Record one base was appended.
        if (recorder_)
            recorder_->baseAppended(contigIndex, appendedBase, next);

        if (next == startNode) {
            result.endNode    = startNode;
            result.isCircular = true;
            size_t overlap = nodeLen - 1; // k-2 chars to trim
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
        if (nextData->getInDegree() != 1 || nextData->getOutDegree() != 1) {
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

    // Phase 1 — branch points and source nodes
    for (NodeId node : graph_.getAllNodes()) {
        const auto* data = graph_.findNode(node);

        bool isBranchPoint = data->getInDegree() > 1 || data->getOutDegree() > 1;
        bool isSource      = data->getInDegree() == 0 && data->getOutDegree() >= 1;

        if (!isBranchPoint && !isSource) continue;

        auto* neighbors = adjCopy_.find(node);

        while (neighbors && !neighbors->empty()) {
            size_t contigIndex = contigs_.size(); // index this contig will occupy
            Contig contig = walkContig(node, contigIndex);

            if (!contig.sequence.empty())
                contigs_.push_back(std::move(contig));
        }
    }

    // Stage 2 — isolated cycles
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