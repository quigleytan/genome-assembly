/*
 * Recorder.h
 * Summary:
 * - Lightweight header-only interface for recording traversal animation steps
 *   into a VisSession during a ContigTraversal or EulerianTraversal run.
 * - ContigTraversal holds an optional Recorder* (default nullptr).
 *   When null, every method is a no-op — zero overhead for normal pipeline runs.
 * - The Recorder does not own the VisSession. The caller is responsible for
 *   keeping the VisSession alive for the duration of the traversal.
 * Usage:
 *   VisSession session;
 *   Recorder recorder(&session);
 *   ContigTraversal ct(graph, &recorder);
 *   ct.computeContigs();
 *   // session.contigSteps now contains the full animation sequence
 */

#ifndef RECORDER_H
#define RECORDER_H

#include <string>
#include "VisData.h"

class Recorder {
private:

    VisSession* session_;

public:

    /**
     * @brief Constructs a Recorder writing into the provided VisSession.
     * @param session Non-owning pointer to the VisSession to record into.
     *                Must remain valid for the lifetime of this Recorder.
     */
    explicit Recorder(VisSession* session) : session_(session) {}

    /**
     * @brief Returns true if this recorder is active (has a valid session).
     * Callers can guard optional recording blocks with this check.
     */
    [[nodiscard]] bool isActive() const { return session_ != nullptr; }

    // ─────────────────────────────────────────────
    // CONTIG ANIMATION EVENTS
    // Called by ContigTraversal during walkContig()
    // ─────────────────────────────────────────────

    /**
     * @brief Records that a new contig walk has begun from startNode.
     *
     * @param contigIndex Index this contig will occupy in VisSession::contigs
     *                    (caller assigns the index before recording begins).
     * @param startNode   Encoded k-1 mer where the walk starts.
     * @param startLabel  Decoded string label for startNode.
     */
    void contigStarted(size_t contigIndex, NodeId startNode, const std::string& startLabel) {
        if (!session_) return;
        TraversalStep step;
        step.type         = TraversalStep::Type::ContigStarted;
        step.contigIndex  = contigIndex;
        step.from         = startNode;
        step.sequence     = startLabel; // Reuse sequence field for label display
        session_->contigSteps.push_back(std::move(step));
    }

    /**
     * @brief Records that one base was appended to the current contig.
     *
     * Called once per node consumed inside walkContig()'s main loop.
     *
     * @param contigIndex Index of the contig being built.
     * @param base        The character appended (last char of the next node's decoded label).
     * @param toNode      The node just consumed.
     */
    void baseAppended(size_t contigIndex, char base, NodeId toNode) {
        if (!session_) return;
        TraversalStep step;
        step.type        = TraversalStep::Type::BaseAppended;
        step.contigIndex = contigIndex;
        step.base        = base;
        step.to          = toNode;
        session_->contigSteps.push_back(std::move(step));
    }

    /**
     * @brief Records that a contig walk has completed.
     *
     * Stores the full finished sequence as a snapshot so the visualizer can
     * display it immediately without replaying every BaseAppended step.
     *
     * @param contigIndex Index of the completed contig.
     * @param endNode     Encoded k-1 mer where the walk terminated.
     * @param endLabel    Decoded string label for endNode.
     * @param sequence    The complete assembled sequence string.
     * @param isCircular  Whether the walk returned to startNode.
     */
    void contigFinished(size_t contigIndex, NodeId endNode,
                        const std::string& endLabel,
                        const std::string& sequence,
                        bool isCircular) {
        if (!session_) return;
        TraversalStep step;
        step.type        = TraversalStep::Type::ContigFinished;
        step.contigIndex = contigIndex;
        step.to          = endNode;
        step.sequence    = sequence;
        step.base        = isCircular ? 'C' : 'L'; // C = circular, L = linear
                                                    // Reuses base field to avoid
                                                    // adding a bool to TraversalStep
        session_->contigSteps.push_back(std::move(step));
    }

// EULERIAN TRAVERSAL EVENTS

    /**
     * @brief Records that Hierholzer's algorithm consumed one directed edge.
     * @param from Source node of the consumed edge.
     * @param to   Destination node of the consumed edge.
     */
    void edgeConsumed(NodeId from, NodeId to) {
        if (!session_) return;
        TraversalStep step;
        step.type = TraversalStep::Type::EdgeConsumed;
        step.from = from;
        step.to   = to;
        session_->eulerSteps.push_back(std::move(step));
    }

    /**
     * @brief Records that a node was committed to the Eulerian path.
     * @param node The node being committed during backtracking.
     */
    void nodeCommitted(NodeId node) {
        if (!session_) return;
        TraversalStep step;
        step.type = TraversalStep::Type::NodeCommitted;
        step.from = node;
        session_->eulerSteps.push_back(std::move(step));
    }

};

#endif // RECORDER_H