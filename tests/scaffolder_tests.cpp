#include <iostream>
#include <vector>

#include "assembler/core/de_bruijn_graph.h"
#include "assembler/core/kmer_table.h"
#include "assembler/core/vis_data.h"
#include "assembler/construction/kmer_encoding.h"
#include "assembler/construction/contig_assembler.h"
#include "assembler/construction/scaffolder.h"

bool LinearChainScaffoldTests();
bool AmbiguousBranchSkipTests();
bool AmbiguousBranchGreedyTests();
bool CircularScaffoldTests();
bool IsolatedContigScaffoldTests();
bool VisSessionExportTests();
bool EndToEndPipelineScaffoldTests();

// Builds a Contig directly from k-1 mer boundary strings, avoiding the need
// to round-trip through a DeBruijnGraph/ContigAssembler for scenarios where
// only Scaffolder's own logic (connection map, resolution strategy, circular
// detection) is under test.
static ContigAssembler::Contig makeContig(const std::string& sequence,
                                          const std::string& startNodeSeq,
                                          const std::string& endNodeSeq) {
    ContigAssembler::Contig contig;
    contig.sequence   = sequence;
    contig.startNode  = KmerEncoding::encode(startNodeSeq);
    contig.endNode    = KmerEncoding::encode(endNodeSeq);
    contig.isCircular = false;
    return contig;
}

int main() {

    bool passed = true;
    bool result;

    result = LinearChainScaffoldTests();
    passed &= result;
    if (result) std::cout << "Linear Chain Scaffold Tests Passed" << std::endl;

    result = AmbiguousBranchSkipTests();
    passed &= result;
    if (result) std::cout << "Ambiguous Branch Skip Tests Passed" << std::endl;

    result = AmbiguousBranchGreedyTests();
    passed &= result;
    if (result) std::cout << "Ambiguous Branch Greedy Tests Passed" << std::endl;

    result = CircularScaffoldTests();
    passed &= result;
    if (result) std::cout << "Circular Scaffold Tests Passed" << std::endl;

    result = IsolatedContigScaffoldTests();
    passed &= result;
    if (result) std::cout << "Isolated Contig Scaffold Tests Passed" << std::endl;

    result = VisSessionExportTests();
    passed &= result;
    if (result) std::cout << "Vis Session Export Tests Passed" << std::endl;

    result = EndToEndPipelineScaffoldTests();
    passed &= result;
    if (result) std::cout << "End-to-End Pipeline Scaffold Tests Passed" << std::endl;

    return passed ? 0 : 1;
}

// -----------------------------------------------------------------------
// Test: Two unambiguous contigs sharing a single boundary node (A's end ==
// B's start, and no other contig touches that node) should be walked into
// a single two-entry scaffold, terminating with an unknown gap after B.
// -----------------------------------------------------------------------
bool LinearChainScaffoldTests() {
    bool passed = true;

    ContigAssembler::Contig a = makeContig("AAACG", "AAA", "ACG");
    ContigAssembler::Contig b = makeContig("ACGTT", "ACG", "GTT");
    std::vector<ContigAssembler::Contig> contigs = {a, b};

    DeBruijnGraph graph(4);
    Scaffolder scaffolder(contigs, graph); // default strategy: skip()
    scaffolder.buildScaffolds();

    const auto& scaffolds = scaffolder.getScaffolds();

    if (scaffolds.size() != 1) {
        passed = false;
        std::cout << "[Linear] Expected 1 scaffold, got " << scaffolds.size() << "\n";
        return passed;
    }

    const Scaffold& scaffold = scaffolds[0];

    if (scaffold.entries.size() != 2) {
        passed = false;
        std::cout << "[Linear] Expected 2 entries, got " << scaffold.entries.size() << "\n";
    }

    if (scaffold.isCircular) {
        passed = false;
        std::cout << "[Linear] Scaffold should not be circular\n";
    }

    if (scaffold.entries.size() == 2) {
        if (scaffold.entries[0].contigIndex != 0 || scaffold.entries[1].contigIndex != 1) {
            passed = false;
            std::cout << "[Linear] Entries out of order\n";
        }
        if (scaffold.entries[0].gapAfter != ScaffoldEntry::DIRECT_OVERLAP) {
            passed = false;
            std::cout << "[Linear] Expected direct overlap between A and B\n";
        }
        if (scaffold.entries[1].gapAfter != ScaffoldEntry::UNKNOWN_GAP) {
            passed = false;
            std::cout << "[Linear] Expected unknown gap after the final contig\n";
        }
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: A ends at a node where two contigs (B, shorter; C, longer) both
// start. With ResolutionStrategy::skip(), any ambiguous boundary halts the
// walk rather than guessing - so A, B, and C should each end up as their
// own single-entry scaffold instead of being chained together.
// -----------------------------------------------------------------------
bool AmbiguousBranchSkipTests() {
    bool passed = true;

    ContigAssembler::Contig a = makeContig("AAACG",     "AAA", "ACG");
    ContigAssembler::Contig b = makeContig("ACGTT",     "ACG", "GTT"); // short arm
    ContigAssembler::Contig c = makeContig("ACGTTTGG",  "ACG", "TGG"); // long arm
    std::vector<ContigAssembler::Contig> contigs = {a, b, c};

    DeBruijnGraph graph(4);
    Scaffolder scaffolder(contigs, graph, ResolutionStrategy::skip());
    scaffolder.buildScaffolds();

    const auto& scaffolds = scaffolder.getScaffolds();

    if (scaffolds.size() != 3) {
        passed = false;
        std::cout << "[Skip] Expected 3 scaffolds, got " << scaffolds.size() << "\n";
    }

    for (const auto& scaffold : scaffolds) {
        if (scaffold.entries.size() != 1) {
            passed = false;
            std::cout << "[Skip] Expected every scaffold to be a single contig\n";
        }
    }

    // Find A's scaffold specifically and confirm it recorded an unresolved gap.
    bool foundAScaffold = false;
    for (const auto& scaffold : scaffolds) {
        if (scaffold.entries.size() == 1 && scaffold.entries[0].contigIndex == 0) {
            foundAScaffold = true;
            if (scaffold.entries[0].gapAfter != ScaffoldEntry::UNKNOWN_GAP) {
                passed = false;
                std::cout << "[Skip] Expected A's ambiguous branch to record an unknown gap\n";
            }
        }
    }
    if (!foundAScaffold) {
        passed = false;
        std::cout << "[Skip] Did not find contig A in any scaffold\n";
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: Same ambiguous branch as above, but with ResolutionStrategy::greedy()
// (skipAmbiguous = false, length-only scoring). The longer contig (C) should
// win the branch and be chained onto A; the shorter contig (B) is left as
// its own isolated scaffold since nothing chose to walk into it.
// -----------------------------------------------------------------------
bool AmbiguousBranchGreedyTests() {
    bool passed = true;

    ContigAssembler::Contig a = makeContig("AAACG",     "AAA", "ACG");
    ContigAssembler::Contig b = makeContig("ACGTT",     "ACG", "GTT"); // short arm (len 5)
    ContigAssembler::Contig c = makeContig("ACGTTTGG",  "ACG", "TGG"); // long arm (len 8)
    std::vector<ContigAssembler::Contig> contigs = {a, b, c};

    DeBruijnGraph graph(4);
    Scaffolder scaffolder(contigs, graph, ResolutionStrategy::greedy());
    scaffolder.buildScaffolds();

    const auto& scaffolds = scaffolder.getScaffolds();

    if (scaffolds.size() != 2) {
        passed = false;
        std::cout << "[Greedy] Expected 2 scaffolds, got " << scaffolds.size() << "\n";
        return passed;
    }

    bool foundChain = false;
    bool foundIsolatedB = false;

    for (const auto& scaffold : scaffolds) {
        if (scaffold.entries.size() == 2 &&
            scaffold.entries[0].contigIndex == 0 &&
            scaffold.entries[1].contigIndex == 2) {
            foundChain = true;
            if (scaffold.entries[0].gapAfter != ScaffoldEntry::DIRECT_OVERLAP) {
                passed = false;
                std::cout << "[Greedy] Expected A->C to be a direct overlap\n";
            }
        }
        if (scaffold.entries.size() == 1 && scaffold.entries[0].contigIndex == 1) {
            foundIsolatedB = true;
        }
    }

    if (!foundChain) {
        passed = false;
        std::cout << "[Greedy] Expected A to chain onto the longer contig C\n";
    }
    if (!foundIsolatedB) {
        passed = false;
        std::cout << "[Greedy] Expected the shorter contig B to be left isolated\n";
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: Three contigs forming a pure cycle (A -> B -> C -> A) with no
// external entry point. None qualify as a linear "scaffold start" in phase 1
// (each has a predecessor), so all three are picked up by phase 2 as a
// single circular scaffold.
// -----------------------------------------------------------------------
bool CircularScaffoldTests() {
    bool passed = true;

    ContigAssembler::Contig a = makeContig("AAACG",  "AAA", "ACG");
    ContigAssembler::Contig b = makeContig("ACGTT",  "ACG", "GTT");
    ContigAssembler::Contig c = makeContig("GTTAAA", "GTT", "AAA");
    std::vector<ContigAssembler::Contig> contigs = {a, b, c};

    DeBruijnGraph graph(4);
    Scaffolder scaffolder(contigs, graph); // skip() - no ambiguity present here
    scaffolder.buildScaffolds();

    const auto& scaffolds = scaffolder.getScaffolds();

    if (scaffolds.size() != 1) {
        passed = false;
        std::cout << "[Circular] Expected 1 scaffold, got " << scaffolds.size() << "\n";
        return passed;
    }

    const Scaffold& scaffold = scaffolds[0];

    if (scaffold.entries.size() != 3) {
        passed = false;
        std::cout << "[Circular] Expected 3 entries, got " << scaffold.entries.size() << "\n";
    }

    if (!scaffold.isCircular) {
        passed = false;
        std::cout << "[Circular] Scaffold should be circular\n";
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: A single contig with no shared boundary nodes forms its own
// one-entry, non-circular scaffold.
// -----------------------------------------------------------------------
bool IsolatedContigScaffoldTests() {
    bool passed = true;

    ContigAssembler::Contig a = makeContig("TTTTGGG", "TTT", "GGG");
    std::vector<ContigAssembler::Contig> contigs = {a};

    DeBruijnGraph graph(4);
    Scaffolder scaffolder(contigs, graph);
    scaffolder.buildScaffolds();

    const auto& scaffolds = scaffolder.getScaffolds();

    if (scaffolds.size() != 1) {
        passed = false;
        std::cout << "[Isolated] Expected 1 scaffold, got " << scaffolds.size() << "\n";
        return passed;
    }

    if (scaffolds[0].entries.size() != 1) {
        passed = false;
        std::cout << "[Isolated] Expected a single-entry scaffold\n";
    }

    if (scaffolds[0].isCircular) {
        passed = false;
        std::cout << "[Isolated] Scaffold should not be circular\n";
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: toVisSession() should populate VisContig/VisScaffold entries that
// mirror the underlying scaffold structure - correct contig count, correct
// scaffold membership indices, decoded node labels, and gap values.
// -----------------------------------------------------------------------
bool VisSessionExportTests() {
    bool passed = true;

    ContigAssembler::Contig a = makeContig("AAACG", "AAA", "ACG");
    ContigAssembler::Contig b = makeContig("ACGTT", "ACG", "GTT");
    std::vector<ContigAssembler::Contig> contigs = {a, b};

    DeBruijnGraph graph(4);
    Scaffolder scaffolder(contigs, graph);
    scaffolder.buildScaffolds();

    VisSession session;
    session.k = 4;
    scaffolder.toVisSession(session, /*nodeLen=*/3);

    if (session.contigs.size() != 2) {
        passed = false;
        std::cout << "[VisSession] Expected 2 contigs, got " << session.contigs.size() << "\n";
        return passed;
    }

    if (session.contigs[0].startLabel != "AAA" || session.contigs[0].endLabel != "ACG") {
        passed = false;
        std::cout << "[VisSession] Contig A labels decoded incorrectly\n";
    }

    if (session.contigs[0].scaffoldIndex != 0 || session.contigs[1].scaffoldIndex != 0) {
        passed = false;
        std::cout << "[VisSession] Both contigs should belong to scaffold 0\n";
    }

    if (session.scaffolds.size() != 1) {
        passed = false;
        std::cout << "[VisSession] Expected 1 scaffold, got " << session.scaffolds.size() << "\n";
        return passed;
    }

    const VisScaffold& vs = session.scaffolds[0];

    if (vs.contigIndices != std::vector<size_t>{0, 1}) {
        passed = false;
        std::cout << "[VisSession] Scaffold contig indices incorrect\n";
    }

    if (vs.gaps != std::vector<int>{ScaffoldEntry::DIRECT_OVERLAP, ScaffoldEntry::UNKNOWN_GAP}) {
        passed = false;
        std::cout << "[VisSession] Scaffold gap values incorrect\n";
    }

    if (vs.isCircular) {
        passed = false;
        std::cout << "[VisSession] Scaffold should not be circular\n";
    }

    return passed;
}

// -----------------------------------------------------------------------
// Test: Full pipeline sanity check - KmerTable -> DeBruijnGraph ->
// ContigAssembler -> Scaffolder on a repeat-free sequence. ContigAssembler
// should produce a single contig spanning the whole sequence, and Scaffolder
// should wrap it in a single non-circular scaffold.
// -----------------------------------------------------------------------
bool EndToEndPipelineScaffoldTests() {
    bool passed = true;

    std::string sequence = "ACGTTTGA";
    int k = 4;

    KmerTable kTable(sequence.length(), k);
    KmerEncoding::encodeSequence(sequence, k, kTable);

    DeBruijnGraph graph(k);
    for (const auto& entry : kTable)
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);

    ContigAssembler contigAssembler(graph);
    contigAssembler.computeContigs();

    Scaffolder scaffolder(contigAssembler.getContigs(), graph);
    scaffolder.buildScaffolds();

    const auto& scaffolds = scaffolder.getScaffolds();

    if (scaffolds.size() != 1) {
        passed = false;
        std::cout << "[EndToEnd] Expected 1 scaffold, got " << scaffolds.size() << "\n";
        return passed;
    }

    if (scaffolds[0].entries.size() != 1) {
        passed = false;
        std::cout << "[EndToEnd] Expected a single-contig scaffold for a repeat-free sequence\n";
    }

    if (scaffolds[0].isCircular) {
        passed = false;
        std::cout << "[EndToEnd] Scaffold should not be circular\n";
    }

    return passed;
}
