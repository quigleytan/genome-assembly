#include <fstream>
#include <iostream>

#include "../data_list.h"
#include "../data_initalization/dna_sequence.h"
#include "../data_initalization/sequence_reader.h"

#include "../encoding/kmer_encoding.h"
#include "../encoding/kmer_table.h"

#include "../scaffolding/de_bruijn_graph.h"
#include "../scaffolding/contig_traversal.h"

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
    ContigTraversal ct(graph);
    ct.computeContigs();
    ct.printStats();
}

// MAIN FUNCTION

int main() {
    try {
        const std::string path = "../Data/" + sequence[0];

        DNASequence genome = loadGenome(path);
        std::cout << genome.getName() << "\n";
        std::cout << "Sequence length: " << genome.getLength() << " bases\n";
        std::cout << "--------------------------------------\n";

        std::vector testCases = {5, 20, 40, 60};

        for (int i = 0 ; i < testCases.size(); i++) {
            int k = testCases[i];
            std::cout << "TEST CASE [" << i + 1 << "]: k = " << k << std::endl;
            DeBruijnGraph graph = buildGraph(genome.getSequence(), k, false);
            assembleContigs(graph);
        }

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}