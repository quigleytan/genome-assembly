#include <iostream>

#include "assembler/construction/kmer_encoding.h"
#include "assembler/core/kmer_table.h"
#include "assembler/core/dna_sequence.h"

bool KmerEncodingBasicTests();
bool KmerEncodingRollTests();
bool KmerEncodingValidateKTests();
bool KmerTableTests();

int main() {

    bool passed = true;
    bool result;

    result = KmerEncodingBasicTests();
    passed &= result;
    if (result) std::cout << "Kmer Encoding Basic Tests Passed" << std::endl;

    result = KmerEncodingRollTests();
    passed &= result;
    if (result) std::cout << "Kmer Encoding Roll Tests Passed" << std::endl;

    result = KmerEncodingValidateKTests();
    passed &= result;
    if (result) std::cout << "Kmer Encoding ValidateK Tests Passed" << std::endl;

    result = KmerTableTests();
    passed &= result;
    if (result) std::cout << "Kmer Table Tests Passed" << std::endl;

    return passed ? 0 : 1;
}

bool KmerEncodingBasicTests() {
    bool passed = true;

    // Should encode as A = 0b00, C = 0b01, G = 0b10, T = 0b11
    if (KmerEncoding::encode("A") != 0b00) {
        passed = false;
        std::cout << "Improper encoding of A" << std::endl;
    }

    if (KmerEncoding::encode("C") != 0b01) {
        passed = false;
        std::cout << "Improper encoding of C" << std::endl;
    }

    if (KmerEncoding::encode("G") != 0b10) {
        passed = false;
        std::cout << "Improper encoding of G" << std::endl;
    }

    if (KmerEncoding::encode("T") != 0b11) {
        passed = false;
        std::cout << "Improper encoding of T" << std::endl;
    }

    if (KmerEncoding::encode("ACGGTGT") != 1723) {
        passed = false;
        std::cout << "Improper encoding of ACGGTGT" << std::endl;
    }

    // decode() should exactly invert encode() for a single base and for a longer k-mer
    if (KmerEncoding::decode(0b00, 1) != "A") {
        passed = false;
        std::cout << "Improper decoding of A" << std::endl;
    }

    if (KmerEncoding::decode(KmerEncoding::encode("ACGGTGT"), 7) != "ACGGTGT") {
        passed = false;
        std::cout << "Encode/decode roundtrip failed for ACGGTGT" << std::endl;
    }

    // A run of leading A's (encoded as 0 bits) must not be dropped on decode -
    // the k-mer length pads the output rather than trimming leading zeros.
    if (KmerEncoding::decode(KmerEncoding::encode("AAC"), 3) != "AAC") {
        passed = false;
        std::cout << "Leading 'A' bases dropped on decode" << std::endl;
    }

    return passed;
}

bool KmerEncodingRollTests() {
    bool passed = true;

    // Rolling one base at a time should match directly encoding the resulting window.
    // Sequence: ACGTAC, k=3 -> windows: ACG, CGT, GTA, TAC
    std::string sequence = "ACGTAC";
    size_t k = 3;

    NodeId kmer = KmerEncoding::encode(sequence.substr(0, k));
    if (kmer != KmerEncoding::encode("ACG")) {
        passed = false;
        std::cout << "Initial encode incorrect for ACG" << std::endl;
    }

    kmer = KmerEncoding::roll(kmer, sequence[3], k); // adds 'T' -> CGT
    if (kmer != KmerEncoding::encode("CGT")) {
        passed = false;
        std::cout << "Roll did not produce CGT" << std::endl;
    }

    kmer = KmerEncoding::roll(kmer, sequence[4], k); // adds 'A' -> GTA
    if (kmer != KmerEncoding::encode("GTA")) {
        passed = false;
        std::cout << "Roll did not produce GTA" << std::endl;
    }

    kmer = KmerEncoding::roll(kmer, sequence[5], k); // adds 'C' -> TAC
    if (kmer != KmerEncoding::encode("TAC")) {
        passed = false;
        std::cout << "Roll did not produce TAC" << std::endl;
    }

    // bitmask(k) should retain exactly the low 2*k bits.
    if (KmerEncoding::bitmask(3) != 0b111111) {
        passed = false;
        std::cout << "bitmask(3) incorrect" << std::endl;
    }

    // Rolling past k characters must not leak bits from beyond the window -
    // the result should be indistinguishable from directly encoding the new window.
    NodeId rolledOnceMore = KmerEncoding::roll(kmer, 'G', k); // TAC -> ACG
    if (rolledOnceMore != KmerEncoding::encode("ACG")) {
        passed = false;
        std::cout << "Roll did not properly mask overflow bits" << std::endl;
    }

    return passed;
}

bool KmerEncodingValidateKTests() {
    bool passed = true;

    // Valid k should be returned unchanged.
    if (KmerEncoding::validateK(3) != 3) {
        passed = false;
        std::cout << "validateK rejected a valid k value" << std::endl;
    }

    if (KmerEncoding::validateK(KmerEncoding::MAX_K_128) != KmerEncoding::MAX_K_128) {
        passed = false;
        std::cout << "validateK rejected the maximum valid k value" << std::endl;
    }

    // k < 2 should throw.
    try {
        KmerEncoding::validateK(1);
        passed = false;
        std::cout << "validateK(1): expected exception, none thrown" << std::endl;
    } catch (const std::invalid_argument&) {
        // Expected
    }

    // k > MAX_K_128 should throw.
    try {
        KmerEncoding::validateK(KmerEncoding::MAX_K_128 + 1);
        passed = false;
        std::cout << "validateK(MAX_K_128 + 1): expected exception, none thrown" << std::endl;
    } catch (const std::invalid_argument&) {
        // Expected
    }

    return passed;
}

bool KmerTableTests() {
    bool passed = true;

    DNASequence genome("Sequence 1", "ACGTACGT");

    KmerTable kTable(genome.getLength(), 3);
    KmerEncoding::encodeSequence(genome.getSequence(), 3, kTable);

    if (kTable.getK() != 3) {
        passed = false;
        std::cout << "Kmer table initialized with incorrect k" << std::endl;
    }

    // ACGTACGT, k=3 -> ACG, CGT, GTA, TAC, ACG, CGT: ACG and CGT each appear twice.
    const size_t* acgCount = kTable.find(KmerEncoding::encode("ACG"));
    if (!acgCount || *acgCount != 2) {
        passed = false;
        std::cout << "ACG count incorrect: expected 2, got "
                  << (acgCount ? *acgCount : 0) << std::endl;
    }

    const size_t* gtaCount = kTable.find(KmerEncoding::encode("GTA"));
    if (!gtaCount || *gtaCount != 1) {
        passed = false;
        std::cout << "GTA count incorrect: expected 1, got "
                  << (gtaCount ? *gtaCount : 0) << std::endl;
    }

    // A k-mer that never occurs should be reported as absent, not as a zero-count entry.
    const size_t* aaaCount = kTable.find(KmerEncoding::encode("AAA"));
    if (aaaCount != nullptr) {
        passed = false;
        std::cout << "Invalid kmer found in table" << std::endl;
    }

    // 6 total k-mer insertions, 4 distinct keys (ACG, CGT, GTA, TAC).
    if (kTable.getNumItems() != 4) {
        passed = false;
        std::cout << "Kmer table item count incorrect: expected 4, got "
                  << kTable.getNumItems() << std::endl;
    }

    return passed;
}
