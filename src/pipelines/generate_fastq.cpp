/*
 * generate_fastq.cpp
 * Summary:
 * - Generates a synthetic FASTQ file for benchmarking the genome assembler.
 * - Produces a random reference sequence of a given length, then samples
 *   overlapping 100bp reads at the specified coverage depth.
 * - Output is valid FASTQ with constant Phred-40 quality scores ('I').
 * - No external dependencies — pure C++17 stdlib.
 *
 * Usage:
 *   ./GenerateFastq                              # defaults: 3.5Mb, 100bp, 30x, seed=42
 *   ./GenerateFastq <size> <read_len> <cov> <seed>
 *
 * Example:
 *   ./GenerateFastq 3500000 100 30 42
 *   → writes data/results/synthetic_3mb.fastq
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <filesystem>
#include <stdexcept>

// ─────────────────────────────────────────────
// CONSTANTS
// ─────────────────────────────────────────────

static constexpr size_t DEFAULT_GENOME_SIZE = 500000 ;
static constexpr size_t DEFAULT_READ_LEN    = 100;
static constexpr size_t DEFAULT_COVERAGE    = 10;        // 10x
static constexpr size_t DEFAULT_SEED        = 42;
static constexpr char   OUTPUT_PATH[]       = "../data/results/synthetic_3mb.fastq";

static const char BASES[] = { 'A', 'C', 'G', 'T' };

// ─────────────────────────────────────────────
// REFERENCE GENERATOR
// ─────────────────────────────────────────────

static std::string generateReference(size_t length, std::mt19937& rng) {
    std::uniform_int_distribution<int> baseDist(0, 3);
    std::string ref;
    ref.reserve(length);
    for (size_t i = 0; i < length; ++i)
        ref += BASES[baseDist(rng)];
    return ref;
}

// ─────────────────────────────────────────────
// READ SAMPLER
// Samples reads uniformly from both strands of
// the reference. Forward strand reads are taken
// directly; reverse complement reads simulate
// paired-end sequencing without needing a mate.
// ─────────────────────────────────────────────

static char complement(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        default:  return 'N';
    }
}

static std::string reverseComplement(const std::string& s) {
    std::string rc(s.rbegin(), s.rend());
    for (char& c : rc) c = complement(c);
    return rc;
}

static void writeReads(std::ostream& out,
                       const std::string& reference,
                       size_t readLen,
                       size_t coverage,
                       std::mt19937& rng) {

    if (reference.length() < readLen)
        throw std::runtime_error("Reference shorter than read length");

    // Number of reads needed to achieve target coverage
    size_t numReads = (reference.length() * coverage) / readLen;

    std::uniform_int_distribution<size_t> posDist(0, reference.length() - readLen);
    std::uniform_int_distribution<int>    strandDist(0, 1);

    // Constant quality string — Phred 40 ('I') for all bases
    std::string quality(readLen, 'I');

    for (size_t i = 0; i < numReads; ++i) {
        size_t pos = posDist(rng);
        std::string seq = reference.substr(pos, readLen);

        // Half the reads come from the reverse complement strand
        if (strandDist(rng) == 1)
            seq = reverseComplement(seq);

        out << "@read_" << i << "\n"
            << seq      << "\n"
            << "+"      << "\n"
            << quality  << "\n";
    }
}

// ─────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────

int main(int argc, char** argv) {

    size_t genomeSize = DEFAULT_GENOME_SIZE;
    size_t readLen    = DEFAULT_READ_LEN;
    size_t coverage   = DEFAULT_COVERAGE;
    size_t seed       = DEFAULT_SEED;

    if (argc >= 2) genomeSize = std::stoull(argv[1]);
    if (argc >= 3) readLen    = std::stoull(argv[2]);
    if (argc >= 4) coverage   = std::stoull(argv[3]);
    if (argc >= 5) seed       = std::stoull(argv[4]);

    std::cout << "Synthetic FASTQ Generator\n";
    std::cout << "--------------------------------------\n";
    std::cout << "Genome size: " << genomeSize  << " bp\n";
    std::cout << "Read length: " << readLen     << " bp\n";
    std::cout << "Coverage:    " << coverage    << "x\n";
    std::cout << "Seed:        " << seed        << "\n";
    std::cout << "Output:      " << OUTPUT_PATH << "\n";
    std::cout << "--------------------------------------\n";

    // Ensure output directory exists
    std::filesystem::create_directories(
        std::filesystem::path(OUTPUT_PATH).parent_path());

    std::mt19937 rng(static_cast<uint32_t>(seed));

    // Generate reference
    std::cout << "Generating reference sequence... ";
    std::cout.flush();
    auto t0 = std::chrono::high_resolution_clock::now();
    std::string reference = generateReference(genomeSize, rng);
    auto t1 = std::chrono::high_resolution_clock::now();
    double refMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "done (" << refMs << " ms)\n";

    // Write reads
    size_t numReads = (genomeSize * coverage) / readLen;
    std::cout << "Writing " << numReads << " reads at " << coverage
              << "x coverage... ";
    std::cout.flush();

    std::ofstream out(OUTPUT_PATH);
    if (!out.is_open())
        throw std::runtime_error("Could not open output file: "
                                 + std::string(OUTPUT_PATH));

    auto t2 = std::chrono::high_resolution_clock::now();
    writeReads(out, reference, readLen, coverage, rng);
    out.close();
    auto t3 = std::chrono::high_resolution_clock::now();
    double writeMs = std::chrono::duration<double, std::milli>(t3 - t2).count();
    std::cout << "done (" << writeMs << " ms)\n";

    // Report file size
    auto fileSize = std::filesystem::file_size(OUTPUT_PATH);
    std::cout << "File size:   "
              << static_cast<double>(fileSize) / (1024.0 * 1024.0)
              << " MB\n";
    std::cout << "Total time:  " << (refMs + writeMs) << " ms\n";
    std::cout << "Done. Run BenchmarkAssembly to profile the assembler.\n";

    return 0;
}