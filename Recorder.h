/*
 * Recorder.h
 * Summary:
 * - Lightweight header-only interface for recording traversal animation steps
 *   into a VisSession during a ContigTraversal.
 * - ContigTraversal holds an optional Recorder* (default nullptr).
 * - The Recorder does not own the VisSession. The caller pipeline is responsible
 *   for keeping the VisSession alive for the duration of the traversal.
 */

#ifndef RECORDER_H
#define RECORDER_H

#include <string>

#include "VisData.h"

class Recorder {

public:

    /**
     * @brief Constructs a Recorder writing into the provided VisSession.
     * @param session Non-owning pointer to the VisSession to record into.
     */
    explicit Recorder(VisSession* session) : session_(session) {}

    /**
     * @brief Returns true if this recorder is active (has a valid session).
     */
    [[nodiscard]] bool isActive() const { return session_ != nullptr; }

    // CONTIG ANIMATION EVENTS

    /**
     * @brief Records that a new contig walk has begun from startNode.
     *
     * @param contigIndex Caller-assigned index to be occupied in VisSession::contigs.
     * @param startNode   Encoded k-1 mer where the walk starts.
     * @param startLabel  Decoded string label for startNode.
     */
    void contigStarted(size_t contigIndex, NodeId startNode, const std::string& startLabel) {
        if (!session_) return;
        TraversalStep step;
        step.type         = TraversalStep::Type::ContigStarted;
        step.contigIndex  = contigIndex;
        step.from         = startNode;
        step.sequence     = startLabel;
        session_->contigSteps.push_back(std::move(step));
    }

    /**
     * @brief Records that a single base was appended to the current contig.
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
        session_->contigSteps.push_back(std::move(step));
    }

private:

    VisSession* session_;

};
#endif // RECORDER_H