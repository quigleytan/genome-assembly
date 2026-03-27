#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "../DataList.h"
#include "../DataInitialization/DNASequence.h"
#include "../DataInitialization/SequenceReader.h"

#include "../DataProcessing/KmerEncoding.h"
#include "../DataProcessing/KmerTable.h"

#include "../ScaffoldConstruction/DeBruijnGraph.h"
#include "../ScaffoldConstruction/ContigScaffolder.h"
#include "../ScaffoldConstruction/ContigTraversal.h"

// CONSTANTS

static constexpr size_t UNKNOWN_GAP_NS   = 10;
static constexpr size_t FASTA_LINE_WIDTH = 60;

static constexpr size_t INTER_SCAFFOLD_NS = 10; // Number of Ns between scaffolds

/**
 * @brief Writes a full pseudo-genome by concatenating scaffolds and
 *        inserting N's between them for unknown gaps.
 *
 * @param scaffolds  Vector of Scaffold objects from ContigScaffolder
 * @param contigs    Original vector of contigs
 * @param genomeName Name for the output genome (FASTA header)
 * @param outputPath Path to write the FASTA file
 */
void writeFullGenomeFasta(
    const std::vector<Scaffold>& scaffolds,
    const std::vector<ContigTraversal::Contig>& contigs,
    const std::string& genomeName,
    const std::string& outputPath)
{
    std::ofstream out(outputPath);
    if (!out.is_open()) {
        throw std::runtime_error("Could not write output file: " + outputPath);
    }

    // Build full genome sequence
    std::string fullGenome;
    for (size_t i = 0; i < scaffolds.size(); ++i) {
        const Scaffold& scaffold = scaffolds[i];

        // Append contigs in scaffold
        for (size_t j = 0; j < scaffold.entries.size(); ++j) {
            const ScaffoldEntry& entry = scaffold.entries[j];
            fullGenome += contigs[entry.contigIndex].sequence;
        }

        // Add N's between scaffolds (but not after the last scaffold)
        if (i + 1 < scaffolds.size()) {
            fullGenome += std::string(INTER_SCAFFOLD_NS, 'N');
        }
    }

    // Write FASTA header
    out << ">" << genomeName
        << " scaffolds=" << scaffolds.size()
        << " gapNs=" << INTER_SCAFFOLD_NS
        << "\n";

    // Wrap sequence at FASTA_LINE_WIDTH
    for (size_t pos = 0; pos < fullGenome.size(); pos += FASTA_LINE_WIDTH) {
        out << fullGenome.substr(pos, FASTA_LINE_WIDTH) << "\n";
    }

    out.close();
    std::cout << "Full pseudo-genome written to: " << outputPath << "\n";
}

// STAGE 1 — Load reads from FASTQ into KmerTable

/**
 * @brief Opens a FASTQ file and encodes all reads into a KmerTable.
 *
 * The KmerTable is sized using totalBases (sum of all read lengths) and k.
 * It must be kept alive for the duration of the pipeline run, as ContigScaffolder
 * holds a const pointer to it for scored strategy frequency calculations.
 *
 * @param path       Path to the FASTQ file.
 * @param k          K-mer size for encoding.
 * @param totalBases Estimated total bases across all reads — used to size the table.
 * @return Populated KmerTable.
 * @throws std::runtime_error if the file cannot be opened or contains no reads.
 */
static KmerTable loadReads(const std::string& path, size_t k, size_t totalBases) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + path);

    KmerTable kTable(totalBases, k);
    SequenceReader::encodeAllReads(file, k, kTable);
    file.close();

    std::cout << "Reads encoded:   " << kTable.getNumItems() << " unique k-mers\n";
    return kTable;
}

// STAGE 2 — Build De Bruijn graph from KmerTable

/**
 * @brief Populates a DeBruijnGraph from an already-encoded KmerTable.
 *
 * Pre-sizes the graph to 2x the number of unique k-mers to reduce rehashing.
 * Each k-mer is inserted once per occurrence count stored in the table.
 *
 * @param kTable Populated KmerTable from loadReads().
 * @param k      K-mer size — must match the k used to build kTable.
 * @return Populated DeBruijnGraph ready for contig traversal.
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

// STAGE 3 — Contig traversal

/**
 * @brief Runs ContigTraversal on the graph and returns the contig list.
 *
 * Prints contig stats to stdout as a summary before returning.
 *
 * @param graph Populated DeBruijnGraph.
 * @return Vector of Contig structs produced by the traversal.
 */
static std::vector<ContigTraversal::Contig> buildContigs(DeBruijnGraph& graph) {
    ContigTraversal ct(graph);
    ct.computeContigs();
    ct.printStats();
    return ct.getContigs();
}

// STAGE 4 — Scaffold and report

/**
 * @brief Runs ContigScaffolder with the given strategy and prints stats.
 *
 * A fresh DeBruijnGraph and ContigTraversal are produced per call since
 * each strategy run is fully independent (rebuild everything per strategy).
 * The KmerTable pointer is passed through for scored strategy frequency scoring.
 *
 * @param contigs    Contigs from Stage 3.
 * @param graph      DeBruijnGraph from Stage 2.
 * @param strategy   Resolution strategy to use.
 * @param kTable     KmerTable for frequency-based scoring (may be nullptr for non-scored).
 * @return Populated ContigScaffolder after buildScaffolds() is called.
 */
static ContigScaffolder buildScaffolds(
    const std::vector<ContigTraversal::Contig>& contigs,
    const DeBruijnGraph& graph,
    ResolutionStrategy strategy,
    const KmerTable* kTable)
{
    ContigScaffolder scaffolder(contigs, graph, strategy, kTable);
    scaffolder.buildScaffolds();
    scaffolder.printStats();
    return scaffolder;
}

// STAGE 5 — Write gap-aware FASTA output

/**
 * @brief Writes scaffolds to a gap-aware FASTA file.
 *
 * Each scaffold is written as a single FASTA entry. Contigs within a scaffold
 * are joined directly (DIRECT_OVERLAP) or separated by 100 N's (UNKNOWN_GAP).
 * Lines are wrapped at FASTA_LINE_WIDTH (60) characters.
 *
 * Output filename format: scaffolds_k{k}_{strategyName}.fna
 *
 * @param scaffolds    Scaffolds produced by ContigScaffolder.
 * @param contigs      Contig list required to recover sequences.
 * @param k            K value used for this run — used in the filename.
 * @param strategyName Human-readable strategy label ("skip", "greedy", "scored").
 */
static void writeScaffoldFasta(
    const std::vector<Scaffold>& scaffolds,
    const std::vector<ContigTraversal::Contig>& contigs,
    size_t k,
    const std::string& strategyName)
{
    const std::string filename =
        "../Data/Results/scaffolds_k" + std::to_string(k) + "_" + strategyName + ".fna";

    std::ofstream out(filename);
    if (!out.is_open())
        throw std::runtime_error("Could not write output file: " + filename);

    for (size_t i = 0; i < scaffolds.size(); ++i) {
        const Scaffold& scaffold = scaffolds[i];

        // Build the full scaffold sequence
        std::string scaffoldSeq;
        for (size_t j = 0; j < scaffold.entries.size(); ++j) {
            const ScaffoldEntry& entry = scaffold.entries[j];
            scaffoldSeq += contigs[entry.contigIndex].sequence;

            // Insert N's for unresolved gaps — but not after the final contig
            if (j + 1 < scaffold.entries.size() &&
                entry.gapAfter == ScaffoldEntry::UNKNOWN_GAP) {
                scaffoldSeq += std::string(UNKNOWN_GAP_NS, 'N');
            }
        }

        // FASTA header
        out << ">Scaffold_" << (i + 1)
            << " contigs=" << scaffold.entries.size()
            << (scaffold.isCircular ? " circular" : "")
            << " strategy=" << strategyName
            << " k=" << k
            << "\n";

        // Sequence wrapped at 60 chars per line
        for (size_t pos = 0; pos < scaffoldSeq.size(); pos += FASTA_LINE_WIDTH)
            out << scaffoldSeq.substr(pos, FASTA_LINE_WIDTH) << "\n";
    }

    out.close();
    std::cout << "Output written:  " << filename << "\n";
}

// MAIN


int main() {
    try {
        // FASTQ file path — update sequence index as needed
        const std::string path = "../Data/" + sequence[7];

        // Estimated total bases — adjust to match your reads file
        // A rough upper bound is fine; KmerTable sizes conservatively
        const size_t estimatedTotalBases = 100000;

        const std::vector<size_t> kValues = {4};

        const std::vector<std::pair<ResolutionStrategy, std::string>> strategies = {
            { ResolutionStrategy::skip(),   "skip"   },
            { ResolutionStrategy::greedy(), "greedy" },
            { ResolutionStrategy::scored(), "scored" },
        };

        for (size_t k : kValues) {
            std::cout << "======================================\n";
            std::cout << "K = " << k << "\n";
            std::cout << "======================================\n";

            for (const auto& [strategy, strategyName] : strategies) {
                std::cout << "--------------------------------------\n";
                std::cout << "Strategy: " << strategyName << "\n";
                std::cout << "--------------------------------------\n";

                // Each run is fully independent — rebuild everything
                KmerTable kTable  = loadReads(path, k, estimatedTotalBases);
                DeBruijnGraph graph = buildGraph(kTable, k);
                auto contigs        = buildContigs(graph);

                // Pass kTable pointer only for scored — skip/greedy don't need it
                const KmerTable* kTablePtr =
                    (strategyName == "scored") ? &kTable : nullptr;

                ContigScaffolder scaffolder =
                    buildScaffolds(contigs, graph, strategy, kTablePtr);

                writeScaffoldFasta(
                    scaffolder.getScaffolds(), contigs, k, strategyName);
            }
        }

        size_t k = 9;
        KmerTable kTable  = loadReads(path, k, estimatedTotalBases);
        DeBruijnGraph graph = buildGraph(kTable, k);
        auto contigs        = buildContigs(graph);

        // Pass kTable pointer only for scored — skip/greedy don't need it
        const KmerTable* kTablePtr = &kTable;

        ContigScaffolder scaffolder =
            buildScaffolds(contigs, graph, ResolutionStrategy::scored(), kTablePtr);
        writeFullGenomeFasta(scaffolder.getScaffolds(), contigs,
    "Escherichia_pseudo_genome",
    "../Data/Results/full_reconstructed_genome.fna"
);

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}