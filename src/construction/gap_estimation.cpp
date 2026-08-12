//
// Created by quigl on 3/27/2026.
//

#include "assembler/construction/gap_estimation.h"

#include <algorithm>
#include <cmath>

#include "assembler/construction/kmer_encoding.h"

// PRIVATE

double GapEstimator::computeMeanFrequency() const {
    if (kmerTable_.getNumItems() == 0) return 0.0;

    size_t total = 0;
    for (const auto& entry : kmerTable_)
        total += entry.value;

    return static_cast<double>(total) / static_cast<double>(kmerTable_.getNumItems());
}

double GapEstimator::computeBoundaryFrequency(const std::string& seq, bool fromStart) const {
    if (seq.length() < k_) return 0.0;

    const size_t maxKmers    = seq.length() - k_ + 1;
    const size_t sampleCount = std::min(config_.windowSize, maxKmers);

    // sampleCount overlapping k-mers span exactly (sampleCount - 1) + k_ bases.
    const size_t windowLen = sampleCount + k_ - 1;
    const std::string window = fromStart
        ? seq.substr(0, windowLen)
        : seq.substr(seq.length() - windowLen);

    NodeId kmer = KmerEncoding::encode(window.substr(0, k_));
    double total = 0.0;
    const size_t* count = kmerTable_.find(kmer);
    total += count ? static_cast<double>(*count) : 0.0;

    for (size_t i = k_; i < window.length(); ++i) {
        kmer  = KmerEncoding::roll(kmer, window[i], k_);
        count = kmerTable_.find(kmer);
        total += count ? static_cast<double>(*count) : 0.0;
    }

    return total / static_cast<double>(sampleCount);
}

// PUBLIC

GapEstimator::GapEstimator(const KmerTable& kmerTable, size_t k, GapEstimationConfig config)
    : kmerTable_(kmerTable),
      k_(k),
      config_(config),
      meanFrequency_(computeMeanFrequency()) {}

size_t GapEstimator::estimateGap(const std::string& upstreamSeq, const std::string& downstreamSeq) const {
    const double upstreamFreq   = computeBoundaryFrequency(upstreamSeq, false);
    const double downstreamFreq = computeBoundaryFrequency(downstreamSeq, true);
    const double edgeFreq       = (upstreamFreq + downstreamFreq) / 2.0;

    // Ratio of 1.0 = boundary coverage matches (or exceeds) the genome-wide
    // mean - no drop detected, so the contigs are likely nearly adjacent.
    // Ratio near 0.0 = boundary coverage has thinned out sharply, implying
    // a larger unresolved region between the two contigs.
    const double ratio = (meanFrequency_ > 0.0)
        ? std::clamp(edgeFreq / meanFrequency_, 0.0, 1.0)
        : 1.0;

    const double rawGap = (1.0 - ratio) * config_.scaleFactor;
    const size_t gap    = static_cast<size_t>(std::llround(rawGap));

    return std::clamp(gap, config_.minGap, config_.maxGap);
}
