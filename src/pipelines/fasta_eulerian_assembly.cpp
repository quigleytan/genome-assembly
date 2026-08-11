#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "assembler/core/kmer_table.h"
#include "assembler/core/de_bruijn_graph.h"
#include "assembler/core/dna_sequence.h"

#include "assembler/construction/kmer_encoding.h"
#include "assembler/construction/eulerian_assembler.h"

#include "assembler/io/sequence_reader.h"
#include "assembler/io/console_input.h"


// HELPER FUNCTION

/**
 * @brief Finds the rotation offset needed to align assembled to the original.
 *
 * Searches for original as a substring of assembled+assembled using KMP.
 * The match position gives the number of bases assembled must be rotated
 * left to align with the original start position.
 *
 * @param original  The original DNA sequence string.
 * @param assembled The assembled string to align.
 * @return Rotation offset in bases, 0 if already aligned, npos if not found.
 */
static size_t findRotationOffset(const std::string& original, const std::string& assembled) {
    if (original.length() != assembled.length()) return std::string::npos;
    if (original == assembled) return 0;

    size_t n = original.length();
    std::string doubled = assembled + assembled;

    std::vector<int> lps(n, 0);
    int len = 0, i = 1;
    while (i < static_cast<int>(n)) {
        if (original[i] == original[len]) { lps[i++] = ++len; }
        else if (len != 0)                { len = lps[len - 1]; }
        else                              { lps[i++] = 0; }
    }

    size_t j = 0;
    for (size_t k = 0; k < doubled.size(); ++k) {
        while (j > 0 && doubled[k] != original[j]) j = lps[j - 1];
        if (doubled[k] == original[j]) ++j;
        if (j == n) return k - n + 1;
    }
    return std::string::npos;
}

// PIPELINE FUNCTIONS

/**
 * @brief Stage 1 - Load a FASTA file into a DNASequence.
 * @param path Path to the FASTA file.
 * @return DNASequence parsed from the file.
 * @throws std::runtime_error if the file cannot be opened or is empty.
 */
static DNASequence loadGenome(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + path);

    auto genomeOpt = SequenceReader::readFasta(file);
    file.close();

    if (!genomeOpt)
        throw std::runtime_error("File was empty: " + path);

    return *genomeOpt;
}

/**
 * @brief Stage 2 - Build a De Bruijn graph from a circular DNA sequence.
 *
 * Circularizes the sequence by appending the first k-1 bases, encodes all
 * k-mers into a kmer_table, then inserts them into a DeBruijnGraph pre-sized
 * to avoid rehashing.
 *
 * @param sequence The raw (non-circularized) DNA sequence string.
 * @param k        K-mer size to use.
 * @return Populated DeBruijnGraph ready for Eulerian traversal.
 */
static DeBruijnGraph buildGraph(const std::string& sequence, int k) {
    std::string circularSequence = sequence + sequence.substr(0, k - 1);

    KmerTable kTable(circularSequence.length(), k);
    KmerEncoding::encodeSequence(circularSequence, k, kTable);

    DeBruijnGraph graph(k, kTable.getNumItems() * 2);
    for (const auto& entry : kTable)
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);

    std::cout << "Unique k-mers:   " << kTable.getNumItems() << "\n";
    std::cout << "Graph built:     " << graph.getNodeCount() << " nodes, "
              << graph.getEdgeCount() << " edges\n";

    return graph;
}

/**
 * @brief Stage 3 - Traverse the graph and reconstruct the genome string.
 *
 * Runs Hierholzer's algorithm, reconstructs the assembled string, then
 * normalizes it back to the original start position using a rotation
 * search via KMP.
 *
 * @param graph            Populated DeBruijnGraph to traverse.
 * @param originalSequence Original (non-circularized) sequence string.
 * @param k                K-mer size used to build the graph.
 * @return Assembled genome string normalized to the original start position.
 */
static std::string assembleGenome(DeBruijnGraph& graph,
                                  const std::string& originalSequence,
                                  int k) {
    EulerianAssembler eulerianPath(graph);
    eulerianPath.computePath();

    std::cout << "Path length:     " << eulerianPath.getPath().size() << " nodes\n";

    std::string assembled = eulerianPath.reconstructGenome(true);
    if (assembled.length() > originalSequence.length())
        assembled = assembled.substr(0, originalSequence.length());

    size_t offset = findRotationOffset(originalSequence, assembled);
    if (offset != std::string::npos && offset != 0)
        assembled = assembled.substr(offset) + assembled.substr(0, offset);

    // Only print offset if rotation was found
    if (offset != std::string::npos)
        std::cout << "Rotation offset: " << offset << " bases\n";
    return assembled;
}

/**
 * @brief Stage 4 - Report assembly results to stdout.
 * @param original  Original DNASequence loaded from file.
 * @param assembled Assembled string produced by assembleGenome().
 */
static void reportResults(const DNASequence& original, const std::string& assembled) {
    bool match = assembled == original.getSequence();
    std::cout << "Assembled length: " << assembled.length() << " bases\n";
    std::cout << "Reconstruction:   " << (match ? "SUCCESSFUL" : "FAILED") << "\n";
    std::cout << "--------------------------------------\n";
}

// MAIN FUNCTION

int main() {
    try {
        const std::string path = ConsoleInput::promptFilePath(
            "FASTA file path: ");

        DNASequence genome = loadGenome(path);
        std::cout << genome.getName() << "\n";
        std::cout << "Sequence length: " << genome.getLength() << " bases\n";
        std::cout << "--------------------------------------\n";

        const size_t maxK = std::min<size_t>(KmerEncoding::MAX_K_128, genome.getLength());

        do {
            size_t k = ConsoleInput::promptSizeT("K-mer size", 2, maxK);

            std::cout << "--------------------------------------\n";
            std::cout << "Kmer size: " << k << "\n";
            DeBruijnGraph graph = buildGraph(genome.getSequence(), static_cast<int>(k));
            std::string assembled = assembleGenome(graph, genome.getSequence(), static_cast<int>(k));
            reportResults(genome, assembled);
        } while (ConsoleInput::promptYesNo("Assemble again with a different k?"));

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}