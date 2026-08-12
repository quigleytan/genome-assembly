#include <iostream>
#include <string>

#include "assembler/core/kmer_table.h"
#include "assembler/core/de_bruijn_graph.h"
#include "assembler/construction/kmer_encoding.h"
#include "assembler/construction/gap_estimation.h"
#include "assembler/construction/gap_filling.h"

bool GapEstimatorNoDropTests();
bool GapEstimatorDropDetectionTests();
bool GapFillerLocalReassemblySuccessTests();
bool GapFillerSearchBoundTests();
bool GapFillerUnreachableTargetTests();

int main() {

    bool passed = true;
    bool result;

    result = GapEstimatorNoDropTests();
    passed &= result;
    if (result) std::cout << "Gap Estimator No Drop Tests Passed" << std::endl;

    result = GapEstimatorDropDetectionTests();
    passed &= result;
    if (result) std::cout << "Gap Estimator Drop Detection Tests Passed" << std::endl;

    result = GapFillerLocalReassemblySuccessTests();
    passed &= result;
    if (result) std::cout << "Gap Filler Local Reassembly Success Tests Passed" << std::endl;

    result = GapFillerSearchBoundTests();
    passed &= result;
    if (result) std::cout << "Gap Filler Search Bound Tests Passed" << std::endl;

    result = GapFillerUnreachableTargetTests();
    passed &= result;
    if (result) std::cout << "Gap Filler Unreachable Target Tests Passed" << std::endl;

    return passed ? 0 : 1;
}

// -----------------------------------------------------------------------
// Test: When both contig boundaries sit on high-frequency k-mers (at or
// above the genome-wide mean), no coverage drop is detected, so the
// estimate should collapse to the configured floor (minGap).
// -----------------------------------------------------------------------
bool GapEstimatorNoDropTests() {
    bool passed = true;
    const size_t k = 3;

    KmerTable table(100, k);
    // "AAA" and "AAC" are high-coverage (10x); "ACG"/"CGT" are low-coverage (1x).
    for (int i = 0; i < 10; ++i) table.insert(KmerEncoding::encode("AAA"));
    for (int i = 0; i < 10; ++i) table.insert(KmerEncoding::encode("AAC"));
    table.insert(KmerEncoding::encode("ACG"));
    table.insert(KmerEncoding::encode("CGT"));

    GapEstimationConfig config;
    config.windowSize = 2;
    config.minGap     = 1;
    config.maxGap     = 50;
    config.scaleFactor = 50.0;

    GapEstimator estimator(table, k, config);

    // Both contigs' boundary windows land entirely on the high-frequency pair.
    std::string upstream   = "AAAC";
    std::string downstream = "AAAC";

    size_t gap = estimator.estimateGap(upstream, downstream);
    if (gap != config.minGap) {
        passed = false;
        std::cout << "[NoDrop] Expected gap to collapse to minGap (" << config.minGap
                  << "), got " << gap << "\n";
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: When both contig boundaries sit on low-frequency k-mers (well below
// the genome-wide mean), a sizable gap should be estimated, and it should be
// strictly larger than the no-drop baseline.
// -----------------------------------------------------------------------
bool GapEstimatorDropDetectionTests() {
    bool passed = true;
    const size_t k = 3;

    KmerTable table(100, k);
    for (int i = 0; i < 10; ++i) table.insert(KmerEncoding::encode("AAA"));
    for (int i = 0; i < 10; ++i) table.insert(KmerEncoding::encode("AAC"));
    table.insert(KmerEncoding::encode("ACG"));
    table.insert(KmerEncoding::encode("CGT"));
    // mean frequency = (10 + 10 + 1 + 1) / 4 = 5.5

    GapEstimationConfig config;
    config.windowSize  = 2;
    config.minGap      = 1;
    config.maxGap       = 50;
    config.scaleFactor = 50.0;

    GapEstimator estimator(table, k, config);

    // Upstream's trailing 2 kmers = ACG, CGT (both freq 1); downstream's
    // leading 2 kmers = ACG, CGT as well.
    std::string upstream   = "AAACGT"; // kmers: AAA, AAC, ACG, CGT
    std::string downstream = "ACGTTT"; // kmers: ACG, CGT, GTT, TTT

    size_t dropGap = estimator.estimateGap(upstream, downstream);

    std::string baselineUpstream   = "AAAC";
    std::string baselineDownstream = "AAAC";
    size_t baselineGap = estimator.estimateGap(baselineUpstream, baselineDownstream);

    if (dropGap <= baselineGap) {
        passed = false;
        std::cout << "[Drop] Expected drop gap (" << dropGap
                  << ") to exceed no-drop baseline (" << baselineGap << ")\n";
    }

    if (dropGap < config.minGap || dropGap > config.maxGap) {
        passed = false;
        std::cout << "[Drop] Gap estimate " << dropGap << " outside configured bounds\n";
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: Builds a small, repeat-free de Bruijn graph from a single sequence
// and verifies GapFiller finds the exact bridging path (via local BFS
// reassembly) between two non-adjacent nodes on that path.
// -----------------------------------------------------------------------
bool GapFillerLocalReassemblySuccessTests() {
    bool passed = true;

    const std::string seq = "AACCGGTTA"; // all 3-mers distinct - no branching
    const size_t k = 4;

    KmerTable kTable(seq.length(), k);
    KmerEncoding::encodeSequence(seq, k, kTable);

    DeBruijnGraph graph(k);
    for (const auto& entry : kTable)
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);

    GapEstimator estimator(kTable, k);
    GapFiller    filler(graph, estimator);

    NodeId upstreamEnd     = KmerEncoding::encode(seq.substr(0, 3)); // "AAC"
    NodeId downstreamStart = KmerEncoding::encode(seq.substr(4, 3)); // "GGT"

    GapFiller::Result result = filler.fillGap(seq, upstreamEnd, seq, downstreamStart);

    if (!result.resolved) {
        passed = false;
        std::cout << "[Reassembly] Expected gap to resolve via local reassembly\n";
        return passed;
    }

    std::string expected = seq.substr(3, 4); // "CGGT"
    if (result.sequence != expected) {
        passed = false;
        std::cout << "[Reassembly] Expected bridging sequence \"" << expected
                  << "\", got \"" << result.sequence << "\"\n";
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: The same connected path as above, but with maxSearchDepth capped
// below the number of steps actually required. GapFiller must respect the
// bound and report unresolved rather than searching past it.
// -----------------------------------------------------------------------
bool GapFillerSearchBoundTests() {
    bool passed = true;

    const std::string seq = "AACCGGTTA";
    const size_t k = 4;

    KmerTable kTable(seq.length(), k);
    KmerEncoding::encodeSequence(seq, k, kTable);

    DeBruijnGraph graph(k);
    for (const auto& entry : kTable)
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);

    GapEstimator estimator(kTable, k);

    GapFillConfig config;
    config.searchPadding  = 0;
    config.maxSearchDepth = 2; // The true path needs 4 steps.
    GapFiller filler(graph, estimator, config);

    NodeId upstreamEnd     = KmerEncoding::encode(seq.substr(0, 3)); // "AAC"
    NodeId downstreamStart = KmerEncoding::encode(seq.substr(4, 3)); // "GGT"

    GapFiller::Result result = filler.fillGap(seq, upstreamEnd, seq, downstreamStart);

    if (result.resolved) {
        passed = false;
        std::cout << "[SearchBound] Expected unresolved result when path exceeds maxSearchDepth\n";
    }

    if (!result.sequence.empty()) {
        passed = false;
        std::cout << "[SearchBound] Expected empty sequence for unresolved gap\n";
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: Searching for a node that never appears in the graph must exhaust
// the search and report unresolved, still carrying a valid gap estimate.
// -----------------------------------------------------------------------
bool GapFillerUnreachableTargetTests() {
    bool passed = true;

    const std::string seq = "AACCGGTTA";
    const size_t k = 4;

    KmerTable kTable(seq.length(), k);
    KmerEncoding::encodeSequence(seq, k, kTable);

    DeBruijnGraph graph(k);
    for (const auto& entry : kTable)
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);

    GapEstimator estimator(kTable, k);
    GapFiller    filler(graph, estimator);

    NodeId upstreamEnd = KmerEncoding::encode(seq.substr(0, 3)); // "AAC"
    NodeId unreachable  = KmerEncoding::encode("TTT");            // never appears in seq

    GapFiller::Result result = filler.fillGap(seq, upstreamEnd, seq, unreachable);

    if (result.resolved) {
        passed = false;
        std::cout << "[Unreachable] Expected no path to an absent target node\n";
    }

    if (result.gapEstimate == 0) {
        passed = false;
        std::cout << "[Unreachable] Expected a nonzero gap estimate fallback\n";
    }

    return passed;
}
