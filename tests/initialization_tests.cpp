#include <iostream>
#include <sstream>

#include "assembler/core/dna_sequence.h"
#include "assembler/io/sequence_reader.h"
#include "assembler/exceptions/DNASequenceException.h"
#include "assembler/construction/kmer_encoding.h"
#include "assembler/core/kmer_table.h"

bool DnaSequenceTests();
bool DnaSequenceExceptionTests();
bool FastaReaderTests();
bool FastqReaderTests();
bool MultiReadFastqTests();
bool MalformedFastqTests();

int main() {

    bool passed = true;
    bool result;

    result = DnaSequenceTests();
    passed &= result;
    if (result) std::cout << "DNA Sequence Tests Passed" << std::endl;

    result = DnaSequenceExceptionTests();
    passed &= result;
    if (result) std::cout << "DNA Sequence Exception Tests Passed" << std::endl;

    result = FastaReaderTests();
    passed &= result;
    if (result) std::cout << "FASTA Reader Tests Passed" << std::endl;

    result = FastqReaderTests();
    passed &= result;
    if (result) std::cout << "FASTQ Reader Tests Passed" << std::endl;

    result = MultiReadFastqTests();
    passed &= result;
    if (result) std::cout << "Multi-Read FASTQ Tests Passed" << std::endl;

    result = MalformedFastqTests();
    passed &= result;
    if (result) std::cout << "Malformed FASTQ Tests Passed" << std::endl;

    return passed ? 0 : 1;
}

bool DnaSequenceTests() {
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

    // All-AT sequence: GC count/percent should be zero, not merely falsy.
    DNASequence atOnly("AT Only", "ATATAT");
    if (atOnly.getGCCount() != 0 || atOnly.getGCPercent() != 0.0) {
        passed = false;
        std::cout << "GC calculation incorrect for an all-AT sequence" << std::endl;
    }

    // All-GC sequence: complement should fully swap and GC percent should be 1.0.
    DNASequence gcOnly("GC Only", "GCGC");
    if (gcOnly.getComplement() != "CGCG" || gcOnly.getGCPercent() != 1.0) {
        passed = false;
        std::cout << "Complement/GC calculation incorrect for an all-GC sequence" << std::endl;
    }

    // Equality/inequality operators compare on sequence content only.
    DNASequence sequenceTwo("Sequence 2", "ACGTACGT");
    DNASequence sequenceThree("Sequence 3", "TTTTTTTT");

    if (!(sequenceOne == sequenceTwo)) {
        passed = false;
        std::cout << "Equality operator failed for identical sequences" << std::endl;
    }

    if (!(sequenceOne != sequenceThree)) {
        passed = false;
        std::cout << "Inequality operator failed for differing sequences" << std::endl;
    }

    return passed;
}

bool DnaSequenceExceptionTests() {
    bool passed = true;

    // Empty sequence should throw.
    try {
        DNASequence empty("Empty", "");
        passed = false;
        std::cout << "Empty sequence: expected exception, none thrown" << std::endl;
    } catch (const DNASequenceException&) {
        // Expected
    }

    // Invalid base should throw.
    try {
        DNASequence invalid("Invalid", "ACGTX");
        passed = false;
        std::cout << "Invalid base: expected exception, none thrown" << std::endl;
    } catch (const DNASequenceException&) {
        // Expected
    }

    // Lowercase bases are not accepted - only uppercase A/C/G/T are valid.
    try {
        DNASequence lowercase("Lowercase", "acgt");
        passed = false;
        std::cout << "Lowercase bases: expected exception, none thrown" << std::endl;
    } catch (const DNASequenceException&) {
        // Expected
    }

    return passed;
}

bool FastaReaderTests() {
    bool passed = true;

    // Standard single-record FASTA, sequence split across multiple lines.
    {
        std::istringstream in(
            ">TestSequence\n"
            "ACGTA\n"
            "TTGAC\n"
        );

        auto genomeOpt = SequenceReader::readFasta(in);
        if (!genomeOpt) {
            passed = false;
            std::cout << "[FASTA] Expected a sequence, got nullopt" << std::endl;
        } else {
            if (genomeOpt->getSequence() != "ACGTATTGAC") {
                passed = false;
                std::cout << "[FASTA] Sequence lines were not concatenated correctly" << std::endl;
            }
            if (genomeOpt->getName() != "TestSequence") {
                passed = false;
                std::cout << "[FASTA] Header name parsed incorrectly" << std::endl;
            }
        }
    }

    // Empty stream should return std::nullopt, not throw.
    {
        std::istringstream in("");
        auto genomeOpt = SequenceReader::readFasta(in);
        if (genomeOpt.has_value()) {
            passed = false;
            std::cout << "[FASTA] Expected nullopt for an empty stream" << std::endl;
        }
    }

    // Missing '>' on the header line should throw.
    {
        std::istringstream in("TestSequence\nACGT\n");
        try {
            SequenceReader::readFasta(in);
            passed = false;
            std::cout << "[FASTA] Missing '>' header: expected exception, none thrown" << std::endl;
        } catch (const std::runtime_error&) {
            // Expected
        }
    }

    // A second '>' record is not supported and should throw.
    {
        std::istringstream in(">First\nACGT\n>Second\nTTTT\n");
        try {
            SequenceReader::readFasta(in);
            passed = false;
            std::cout << "[FASTA] Multiple records: expected exception, none thrown" << std::endl;
        } catch (const std::runtime_error&) {
            // Expected
        }
    }

    return passed;
}

bool FastqReaderTests() {
    bool passed = true;

    // Standard single-record FASTQ.
    {
        std::istringstream in(
            "@Read1\n"
            "ACGTACGT\n"
            "+\n"
            "IIIIIIII\n"
        );

        auto readOpt = SequenceReader::readFastq(in);
        if (!readOpt) {
            passed = false;
            std::cout << "[FASTQ] Expected a sequence, got nullopt" << std::endl;
        } else {
            if (readOpt->getSequence() != "ACGTACGT") {
                passed = false;
                std::cout << "[FASTQ] Sequence parsed incorrectly" << std::endl;
            }
            if (readOpt->getName() != "Read1") {
                passed = false;
                std::cout << "[FASTQ] Header name parsed incorrectly" << std::endl;
            }
        }
    }

    // Empty stream should return std::nullopt, not throw.
    {
        std::istringstream in("");
        auto readOpt = SequenceReader::readFastq(in);
        if (readOpt.has_value()) {
            passed = false;
            std::cout << "[FASTQ] Expected nullopt for an empty stream" << std::endl;
        }
    }

    // Missing '@' on the header line should throw.
    {
        std::istringstream in("Read1\nACGT\n+\nIIII\n");
        try {
            SequenceReader::readFastq(in);
            passed = false;
            std::cout << "[FASTQ] Missing '@' header: expected exception, none thrown" << std::endl;
        } catch (const std::runtime_error&) {
            // Expected
        }
    }

    return passed;
}

bool MultiReadFastqTests() {
    bool passed = true;

    // Total bases across 3 reads of length 8 = 24
    // Using k=3: each read yields 6 k-mers, 18 total insertions
    const size_t k = 3;
    const size_t totalBases = 24;

    std::istringstream in(
        "@Read1\n"
        "ACGTACGT\n"
        "+\n"
        "IIIIIIII\n"
        "@Read2\n"
        "TGCATGCA\n"
        "+\n"
        "IIIIIIII\n"
        "@Read3\n"
        "ACGTTTGA\n"
        "+\n"
        "IIIIIIII\n"
    );

    KmerTable kTable(totalBases, k);
    SequenceReader::encodeAllReads(in, k, kTable);

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

bool MalformedFastqTests() {
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

    // Truncated record (sequence line missing entirely) — should throw
    std::istringstream truncated("@Read1\n");
    try {
        KmerTable kTable(8, k);
        SequenceReader::encodeAllReads(truncated, k, kTable);
        passed = false;
        std::cout << "Truncated record: expected exception, none thrown" << std::endl;
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
