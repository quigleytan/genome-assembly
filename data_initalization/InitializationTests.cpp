#include "dna_sequence.h"
#include "sequence_reader.h"
#include "encoding/kmer_encoding.h"
#include "encoding/kmer_table.h"
#include <fstream>
#include <sstream>
#include <iostream>

bool testDnaSequence();
bool testSequenceReader();
bool testMultiReadFastq();
bool testMalformedFastq();

int main() {

    if (testDnaSequence()) {
        std::cout << "DNA Sequence Tests Passed" << std::endl;
    }

    if (testSequenceReader()) {
        std::cout << "Sequence Reader Tests Passed" << std::endl;
    }

    if (testMultiReadFastq()) {
        std::cout << "Multi-Read FASTQ Tests Passed" << std::endl;
    }

    if (testMalformedFastq()) {
        std::cout << "Malformed FASTQ Tests Passed" << std::endl;
    }

    return 0;
}

bool testDnaSequence() {

    bool passed = true;

    DNASequence sequenceOne("Sequence 1", "ACGTACGT");

    if (sequenceOne.getSequence() != "ACGTACGT") {
        passed = false;
        std::cout << "Sequence initialization failed" << std::endl;
    }

    if (sequenceOne.getComplement() != "TGCATGCA") {
        passed = false;
        std::cout << "Reverse complement creation failed" << std::endl;
    }

    if (sequenceOne.getName() != "Sequence 1") {
        passed = false;
        std::cout << "Sequence name initialization failed" << std::endl;
    }

    if (sequenceOne.getLength() != 8) {
        passed = false;
        std::cout << "Sequence length calculation failed" << std::endl;
    }

    if (sequenceOne.getGCCount() != 4) {
        passed = false;
        std::cout << "GC count calculation failed" << std::endl;
    }

    if (sequenceOne.getGCPercent() != 0.5) {
        passed = false;
        std::cout << "GC percentage calculation failed" << std::endl;
    }

    return passed;
};

bool testSequenceReader() {
    bool passed = true;

    std::ifstream testFile("../Data/small_test.fna");

    auto genomeOpt = SequenceReader::readFasta(testFile);
    DNASequence genome = *genomeOpt;

    testFile.close();

    if (genome.getSequence() != "ACGTATTGAC") {
        passed = false;
    }

    if (genome.getName() != "TestSequence") {}

    return passed;
};

bool testMultiReadFastq() {
    bool passed = true;

    // Total bases across 3 reads of length 8 = 24
    // Using k=3: each read yields 6 k-mers, 18 total insertions
    const size_t k = 3;
    const size_t totalBases = 24;

    std::ifstream testFile("../Data/test_reads.fastq");

    if (!testFile.is_open()) {
        std::cout << "Could not open test_reads.fastq" << std::endl;
        return false;
    }

    KmerTable kTable(totalBases, k);
    SequenceReader::encodeAllReads(testFile, k, kTable);
    testFile.close();

    // ACG appears in Read1 (ACGTACGT, twice) and Read3 (ACGTTTGA, once) = 3 times
    const size_t* acgCount = kTable.find(KmerEncoding::encode("ACG"));
    if (!acgCount || *acgCount != 3) {
        passed = false;
        std::cout << "ACG count incorrect: expected 3, got "
                  << (acgCount ? *acgCount : 0) << std::endl;
    }

    // TGC appears in Read2 (TGCATGCA, twice) = 2 times
    const size_t* tgcCount = kTable.find(KmerEncoding::encode("TGC"));
    if (!tgcCount || *tgcCount != 2) {
        passed = false;
        std::cout << "TGC count incorrect: expected 2, got "
                  << (tgcCount ? *tgcCount : 0) << std::endl;
    }

    // AAA should not appear in any read
    const size_t* aaaCount = kTable.find(KmerEncoding::encode("AAA"));
    if (aaaCount != nullptr) {
        passed = false;
        std::cout << "AAA should not be present in table" << std::endl;
    }

    // Table should have more than 0 items
    if (kTable.getNumItems() == 0) {
        passed = false;
        std::cout << "kmer_table is empty after encoding reads" << std::endl;
    }

    return passed;
}

bool testMalformedFastq() {
    bool passed = true;

    const size_t k = 3;

    // Missing '+' line — should throw
    std::istringstream missingPlus(
        "@Read1\n"
        "ACGTACGT\n"
        "IIIIIIII\n"   // quality where '+' should be
    );
    try {
        KmerTable kTable(8, k);
        SequenceReader::encodeAllReads(missingPlus, k, kTable);
        passed = false;
        std::cout << "Missing '+' line: expected exception, none thrown" << std::endl;
    } catch (const std::runtime_error&) {
        // Expected
    }

    // Mismatched quality length — should throw
    std::istringstream badQuality(
        "@Read1\n"
        "ACGTACGT\n"
        "+\n"
        "III\n"   // Too short
    );
    try {
        KmerTable kTable(8, k);
        SequenceReader::encodeAllReads(badQuality, k, kTable);
        passed = false;
        std::cout << "Bad quality length: expected exception, none thrown" << std::endl;
    } catch (const std::runtime_error&) {
        // Expected
    }

    // Empty file — should throw from encodeAllReads
    std::istringstream emptyFile("");
    try {
        KmerTable kTable(1, k);
        SequenceReader::encodeAllReads(emptyFile, k, kTable);
        passed = false;
        std::cout << "Empty file: expected exception, none thrown" << std::endl;
    } catch (const std::runtime_error&) {
        // Expected
    }

    return passed;
}