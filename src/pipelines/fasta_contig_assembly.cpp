#include <algorithm>
#include <fstream>
#include <iostream>

#include "assembler/core/dna_sequence.h"
#include "assembler/core/kmer_table.h"
#include "assembler/core/de_bruijn_graph.h"

#include "assembler/construction/kmer_encoding.h"
#include "assembler/construction/contig_assembler.h"

#include "assembler/io/sequence_reader.h"
#include "assembler/io/console_input.h"


// PIPELINE STAGES

/**
 * @brief Loads a FASTA file into a DNASequence.
 * @param path Path to the FASTA file.
 * @return Parsed DNASequence.
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
 * @brief Encodes a DNA sequence into a DeBruijn graph.
 * @param seq Input DNA sequence string.
 * @param k K-mer size to use for encoding.
 * @param circular If true, circularizes the sequence before encoding.
 * @return Populated DeBruijnGraph ready for traversal.
 */
static DeBruijnGraph buildGraph(const std::string& seq, int k, bool circular = false) {
    std::string inputSequence = circular
        ? seq + seq.substr(0, k - 1)
        : seq;

    KmerTable kTable(inputSequence.length(), k);
    KmerEncoding::encodeSequence(inputSequence, k, kTable);

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
 * @brief Runs contig traversal on the graph and reports results.
 * @param graph Populated DeBruijnGraph to traverse.
 * @param sequence Original DNASequence for reference.
 */
static void assembleContigs(DeBruijnGraph& graph) {
    ContigAssembler ct(graph);
    ct.computeContigs();
    ct.printStats();
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
            DeBruijnGraph graph = buildGraph(genome.getSequence(), static_cast<int>(k), false);
            assembleContigs(graph);
        } while (ConsoleInput::promptYesNo("Assemble again with a different k?"));

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}