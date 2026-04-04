#include <fstream>
#include <iostream>
#include <vector>

#include "assembler/core/data_list.h"
#include "assembler/core/kmer_table.h"
#include "assembler/core/de_bruijn_graph.h"

#include "assembler/construction/contig_assembler.h"

#include "assembler/io/sequence_reader.h"

// CONFIGURATION

static constexpr size_t ESTIMATED_TOTAL_BASES = 100000;

// PIPELINE STAGES

/**
 * @brief Stage 1 — Load FASTQ reads into a KmerTable.
 *
 * The returned kmer_table must stay alive through buildGraph(),
 * since the file stream is fully consumed here and cannot be re-read.
 *
 * @param path Path to the FASTQ file.
 * @param k    K-mer size for encoding.
 * @return Populated KmerTable.
 */
static KmerTable loadReads(const std::string& path, size_t k) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + path);

    KmerTable kTable(ESTIMATED_TOTAL_BASES, k);
    SequenceReader::encodeAllReads(file, k, kTable);
    file.close();

    std::cout << "Unique k-mers:   " << kTable.getNumItems() << "\n";
    return kTable;
}

/**
 * @brief Stage 2 — Build a DeBruijnGraph from an encoded KmerTable.
 *
 * Pre-sized to 2x unique k-mer count to reduce rehashing.
 *
 * @param kTable Populated KmerTable from loadReads().
 * @param k      K-mer size — must match the k used in loadReads().
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
 * @brief Stage 3 — Run contig traversal and print stats.
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
        const std::string path = "../Data/" + sequence[7];
        std::cout << "File: " << sequence[7] << "\n";

        const std::vector<size_t> kValues = { 9, 11, 13, 15 };

        for (size_t i = 0; i < kValues.size(); ++i) {
            size_t k = kValues[i];
            std::cout << "======================================\n";
            std::cout << "TEST CASE [" << i + 1 << "]: k = " << k << "\n";
            std::cout << "======================================\n";

            KmerTable kTable  = loadReads(path, k);
            DeBruijnGraph graph = buildGraph(kTable, k);
            assembleContigs(graph);
        }

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}