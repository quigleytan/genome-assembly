#include "EulerianTraversal.h"

#include <stack>
#include <algorithm>

#include "DataProcessing/KmerEncoding.h"

// PRIVATE HELPER FUNCTIONS

void EulerianTraversal::initializeAdjacency() {

    for (NodeId node : graph_.getAllNodes()) {
        adjCopy_.insert(node);
    }

    for (NodeId node : graph_.getAllNodes()) {
        const auto* data = graph_.findNode(node);
        auto* neighborRef = adjCopy_.find(node);  // find() returns pointer, stable after inserts done
        *neighborRef = data->getNeighbors();
        std::sort(neighborRef->begin(), neighborRef->end());
    }
}

NodeId EulerianTraversal::findStartNode() const {

    auto nodes = graph_.getAllNodes();
    NodeId startNode = nodes.front(); // Protects in case of cycle.

    int startCount = 0;
    int endCount = 0;

    for (NodeId node : nodes) {

        const auto* data = graph_.findNode(node);
        int diff = static_cast<int>(data->getOutDegree()) - static_cast<int>(data->getInDegree());

        if (diff == 1) {
            startNode = node;
            startCount++;
        } else if (diff == -1) {
            endCount++;
        } else if (diff != 0) {
            throw std::runtime_error("Graph is not Eulerian");
        }
    }
    // Checking Eulerian conditions
    if (!((startCount == 1 && endCount == 1) || (startCount == 0 && endCount == 0)))
        throw std::runtime_error("Graph does not contain an Eulerian path");

    return startNode; // Returns a valid node in both cycle and path.
}

// PUBLIC

EulerianTraversal::EulerianTraversal(DeBruijnGraph& g) : graph_(g), adjCopy_(g.getNodeCount() * 2) {}

void EulerianTraversal::computePath() {

    path_.clear();
    initializeAdjacency();

    std::stack<NodeId> stack;

    NodeId start = findStartNode();
    stack.push(start);

    // Main traversal loop.
    while (!stack.empty()) {

        // Accesses the current node and retrieve its remaining unused edges.
        NodeId currentID = stack.top();
        auto* neighbors = adjCopy_.find(currentID);

        // Checks if the current node has neighbors, then takes an edge and moves to the next node.
        if (!neighbors->empty()) {
            NodeId next = neighbors->back();
            neighbors->pop_back();   // Deletes the edge immediately.
            stack.push(next);
        } else { // If no edges remain, adds the node to the path and backtrack.
            path_.push_back(currentID);
            stack.pop();
        }
    }

    std::reverse(path_.begin(), path_.end());
}

const std::vector<NodeId>& EulerianTraversal::getPath() const {
    return path_;
}

std::string EulerianTraversal::reconstructGenome(bool isCircuit) const {
    if (path_.empty())
        throw std::runtime_error("Path is empty — call computePath() first");

    const size_t nodeLen = graph_.getK() - 1;
    std::string genome = KmerEncoding::decode(path_.front(), nodeLen);

    // Circuit case - exclude the last node (identical to the first).
    size_t end = isCircuit ? path_.size() - 1 : path_.size();

    for (size_t i = 1; i < end; ++i) {
        genome += KmerEncoding::decode(path_[i], nodeLen).back();
    }
    return genome;
}
