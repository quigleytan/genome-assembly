//
// Created by quigl on 3/27/2026.
//

#include "assembler/construction/gap_filling.h"

#include <algorithm>
#include <queue>

#include "assembler/core/open_addressing_table.h"
#include "assembler/construction/kmer_encoding.h"

// PRIVATE

std::optional<std::string> GapFiller::localReassembly(NodeId from, NodeId to, size_t maxDepth) const {
    if (from == to) return std::string();

    struct Frontier {
        NodeId      node;
        std::string path;
        size_t      depth;
    };

    std::queue<Frontier>   frontier;
    HashTable<NodeId, bool> visited(101);

    frontier.push({from, "", 0});
    visited.insert(from);

    const size_t nodeLen = graph_.getK() - 1;

    while (!frontier.empty()) {
        Frontier current = std::move(frontier.front());
        frontier.pop();

        if (current.depth >= maxDepth) continue;

        const auto* nodeData = graph_.findNode(current.node);
        if (!nodeData) continue;

        for (const auto& edge : nodeData->getNeighbors()) {
            auto [_, inserted] = visited.insert(edge.to);
            if (!inserted) continue; // Already reached via an equal-or-shorter path.

            char appendedBase = KmerEncoding::decode(edge.to, nodeLen).back();
            std::string nextPath = current.path + appendedBase;

            if (edge.to == to)
                return nextPath;

            frontier.push({edge.to, std::move(nextPath), current.depth + 1});
        }
    }

    return std::nullopt;
}

// PUBLIC

GapFiller::GapFiller(const DeBruijnGraph& graph, const GapEstimator& estimator, GapFillConfig config)
    : graph_(graph), estimator_(estimator), config_(config) {}

GapFiller::Result GapFiller::fillGap(const std::string& upstreamSeq, NodeId upstreamEndNode,
                                     const std::string& downstreamSeq, NodeId downstreamStartNode) const {
    Result result;
    result.gapEstimate = estimator_.estimateGap(upstreamSeq, downstreamSeq);

    const size_t searchDepth = std::min(result.gapEstimate + config_.searchPadding, config_.maxSearchDepth);

    std::optional<std::string> path = localReassembly(upstreamEndNode, downstreamStartNode, searchDepth);
    if (path.has_value()) {
        result.sequence = std::move(*path);
        result.resolved = true;
    }

    return result;
}
