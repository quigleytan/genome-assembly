/*
 * gap_estimation.h
 * Created by Tanner Quigley on 3/27/26
 * Summary:
 * - Estimates the number of missing bases between two contigs that the
 *   Scaffolder joined with an UNKNOWN_GAP (no shared boundary node found
 *   in the de Bruijn graph).
 * - Estimation is based on k-mer frequency drop-off near each contig's
 *   terminal boundary: read coverage thins out approaching a true,
 *   unsequenced gap, so a lower boundary-to-genome-average frequency
 *   ratio implies a larger unresolved region.
 * Important notes:
 * - Requires the KmerTable used to build the DeBruijnGraph (raw k-mer
 *   counts), not the graph itself - the estimate is coverage based, not
 *   topology based.
 * - This is a heuristic, not a physical measurement: without paired-end/
 *   mate-pair insert size data, the true gap length is unknowable from
 *   single-end reads alone. Treat estimates as a sizing hint for
 *   GapFiller's search bound and N-padding, not a precise distance.
 */

#ifndef GAP_ESTIMATION_H
#define GAP_ESTIMATION_H

#include <cstddef>
#include <string>

#include "assembler/core/kmer_table.h"

// ESTIMATION CONFIGURATION
//
// Kept as a free struct (not nested in GapEstimator) because a nested class's
// default member initializers cannot be used by a default constructor-argument
// declared earlier in the same enclosing class - see ScoringWeights/
// ResolutionStrategy in scaffolder.h for the same pattern.

struct GapEstimationConfig {
    size_t windowSize  = 5;     // Number of boundary k-mers sampled at each contig end.
    size_t minGap      = 1;     // Floor for a detected (non-adjacent) gap, in bases.
    size_t maxGap      = 500;   // Ceiling clamp for pathological low-coverage estimates.
    double scaleFactor = 100.0; // Bases attributed to a full (ratio 0) frequency drop.
};

class GapEstimator {

public:

    /**
     * @brief Constructor for GapEstimator.
     *
     * Scans kmerTable once to compute the genome-wide mean k-mer frequency,
     * which serves as the "expected coverage" baseline that boundary
     * frequencies are compared against.
     *
     * @param kmerTable Raw k-mer counts from the reads used to build the graph.
     * @param k         K-mer size used to populate kmerTable.
     * @param config    Tunable estimation parameters.
     */
    GapEstimator(const KmerTable& kmerTable, size_t k, GapEstimationConfig config = GapEstimationConfig());

    /**
     * @brief Estimates the number of unresolved bases between two contigs.
     *
     * Samples the last windowSize k-mers of upstreamSeq and the first
     * windowSize k-mers of downstreamSeq, averages their frequencies from
     * kmerTable_, and compares that boundary frequency against the
     * genome-wide mean. A larger drop (lower boundary/mean ratio) maps to
     * a larger estimated gap via scaleFactor, then clamps to
     * [config_.minGap, config_.maxGap].
     *
     * @param upstreamSeq   Full sequence of the contig ending the gap.
     * @param downstreamSeq Full sequence of the contig starting after the gap.
     * @return Estimated gap length in bases, in [config_.minGap, config_.maxGap].
     */
    [[nodiscard]] size_t estimateGap(const std::string& upstreamSeq,
                                      const std::string& downstreamSeq) const;

    /**
     * @brief Returns the genome-wide mean k-mer frequency computed at construction.
     */
    [[nodiscard]] double meanFrequency() const { return meanFrequency_; }

private:

    const KmerTable& kmerTable_;
    size_t k_;
    GapEstimationConfig config_;
    double meanFrequency_;

    [[nodiscard]] double computeMeanFrequency() const;

    /**
     * @brief Averages k-mer frequency over a boundary window of a sequence.
     *
     * @param seq       Sequence to sample from.
     * @param fromStart True to sample the first windowSize k-mers, false to
     *                  sample the last windowSize k-mers.
     * @return Average frequency of the sampled window (0.0 if seq is shorter
     *         than k_, since no complete k-mer can be formed).
     */
    [[nodiscard]] double computeBoundaryFrequency(const std::string& seq, bool fromStart) const;
};

#endif //GAP_ESTIMATION_H
