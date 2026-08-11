#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "assembler/core/kmer_table.h"
#include "assembler/core/de_bruijn_graph.h"

#include "assembler/construction/contig_assembler.h"

#include "assembler/io/sequence_reader.h"
#include "assembler/io/console_input.h"

// PIPELINE STAGES

/**
 * @brief Pre-scans a FASTQ file to sum every read's sequence length, so the
 *        KmerTable can be sized correctly without asking the user to guess
 *        an estimate (the table does not support rehashing mid-run).
 *
 * @param path Path to the FASTQ file.
 * @return Total base count across all reads.
 */
static size_t countFastqBases(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + path);

    size_t totalBases = 0;
    std::string line;
    size_t lineIndex = 0;
    while (std::getline(file, line)) {
        if (lineIndex % 4 == 1) totalBases += line.length();
        ++lineIndex;
    }
    return totalBases;
}

/**
 * @brief Stage 1 - Load FASTQ reads into a KmerTable.
 *
 * The returned kmer_table must stay instantiated through buildGraph(),
 * since the file stream is fully consumed here and cannot be re-read.
 *
 * @param path       Path to the FASTQ file.
 * @param k          K-mer size for encoding.
 * @param totalBases Total base count across all reads, from countFastqBases().
 * @return Populated KmerTable.
 */
static KmerTable loadReads(const std::string& path, size_t k, size_t totalBases) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + path);

    KmerTable kTable(totalBases, k);
    SequenceReader::encodeAllReads(file, k, kTable);
    file.close();

    std::cout << "Unique k-mers:   " << kTable.getNumItems() << "\n";
    return kTable;
}

/**
 * @brief Stage 2 - Build a DeBruijnGraph from an encoded KmerTable.
 *
 * Pre-sized to 2x unique k-mer count to reduce rehashing.
 *
 * @param kTable Populated KmerTable from loadReads().
 * @param k      K-mer size - must match the k used in loadReads().
 * @return Populated DeBruijnGraph.
 */
static DeBruijnGraph buildGraph(const KmerTable& kTable, size_t k) {
    DeBruijnGraph graph(k, kTable.getNumItems() * 2);

    for (const auto& entry : kTable)
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);

    std::cout << "Graph built:     " << graph.getNodeCount() << " nodes, "
              << graph.getEdgeCount() << " edges\n";
    return graph;
}

/**
 * @brief Stage 3 - Run contig traversal and print stats.
 *
 * @param graph Populated DeBruijnGraph.
 */
static void assembleContigs(DeBruijnGraph& graph) {
    ContigAssembler ct(graph);
    ct.computeContigs();
    ct.printStats();
}

// MAIN

int main() {
    try {
        const std::string path = ConsoleInput::promptFilePath(
            "FASTQ file path: ");

        const size_t totalBases = countFastqBases(path);
        std::cout << "Total bases:     " << totalBases << "\n";
        std::cout << "--------------------------------------\n";

        const size_t maxK = std::min<size_t>(63, totalBases);

        do {
            size_t k = ConsoleInput::promptSizeT("K-mer size", 2, maxK);

            std::cout << "--------------------------------------\n";
            std::cout << "Kmer size: " << k << "\n";

            KmerTable kTable  = loadReads(path, k, totalBases);
            DeBruijnGraph graph = buildGraph(kTable, k);
            assembleContigs(graph);
        } while (ConsoleInput::promptYesNo("Assemble again with a different k?"));

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}