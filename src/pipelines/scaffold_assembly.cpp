#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "assembler/core/kmer_table.h"
#include "assembler/core/de_bruijn_graph.h"
#include "assembler/core/vis_data.h"
#include "assembler/core/recorder.h"
#include "assembler/core/data_list.h"

#include "assembler/construction/scaffolder.h"
#include "assembler/construction/contig_assembler.h"

#include "assembler/graphics/vis_exporter.h"

#include "assembler/io/sequence_reader.h"

// CONSTANTS


static constexpr size_t UNKNOWN_GAP_NS    = 10;
static constexpr size_t INTER_SCAFFOLD_NS = 10;
static constexpr size_t FASTA_LINE_WIDTH  = 60;


// STAGE 1 — Load reads

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

// STAGE 2 — Build graph

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

// Normal pipeline run - no recording
static std::vector<ContigAssembler::Contig> buildContigs(DeBruijnGraph& graph) {
    ContigAssembler ct(graph);
    ct.computeContigs();
    ct.printStats();
    return ct.getContigs();
}

// Visualization run - records animation steps into session
static std::vector<ContigAssembler::Contig> buildContigs(DeBruijnGraph& graph,
                                                          Recorder& recorder) {
    ContigAssembler ct(graph, &recorder);
    ct.computeContigs();
    ct.printStats();
    return ct.getContigs();
}


// STAGE 4 — Scaffold


static Scaffolder buildScaffolds(
    const std::vector<ContigAssembler::Contig>& contigs,
    const DeBruijnGraph& graph,
    ResolutionStrategy strategy,
    const KmerTable* kTable)
{
    Scaffolder scaffolder(contigs, graph, strategy, kTable);
    scaffolder.buildScaffolds();
    scaffolder.printStats();
    return scaffolder;
}

// STAGE 5a — Write per-scaffold FASTA

static void writeScaffoldFasta(
    const std::vector<Scaffold>& scaffolds,
    const std::vector<ContigAssembler::Contig>& contigs,
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

        std::string scaffoldSeq;
        for (size_t j = 0; j < scaffold.entries.size(); ++j) {
            const ScaffoldEntry& entry = scaffold.entries[j];
            scaffoldSeq += contigs[entry.contigIndex].sequence;

            if (j + 1 < scaffold.entries.size() &&
                entry.gapAfter == ScaffoldEntry::UNKNOWN_GAP)
                scaffoldSeq += std::string(UNKNOWN_GAP_NS, 'N');
        }

        out << ">Scaffold_" << (i + 1)
            << " contigs=" << scaffold.entries.size()
            << (scaffold.isCircular ? " circular" : "")
            << " strategy=" << strategyName
            << " k=" << k << "\n";

        for (size_t pos = 0; pos < scaffoldSeq.size(); pos += FASTA_LINE_WIDTH)
            out << scaffoldSeq.substr(pos, FASTA_LINE_WIDTH) << "\n";
    }

    out.close();
    std::cout << "Output written:  " << filename << "\n";
}


// STAGE 5b — Write full pseudo-genome FASTA

static void writeFullGenomeFasta(
    const std::vector<Scaffold>& scaffolds,
    const std::vector<ContigAssembler::Contig>& contigs,
    const std::string& genomeName,
    const std::string& outputPath)
{
    std::ofstream out(outputPath);
    if (!out.is_open())
        throw std::runtime_error("Could not write output file: " + outputPath);

    std::string fullGenome;
    for (size_t i = 0; i < scaffolds.size(); ++i) {
        for (const auto& entry : scaffolds[i].entries)
            fullGenome += contigs[entry.contigIndex].sequence;

        if (i + 1 < scaffolds.size())
            fullGenome += std::string(INTER_SCAFFOLD_NS, 'N');
    }

    out << ">" << genomeName
        << " scaffolds=" << scaffolds.size()
        << " gapNs=" << INTER_SCAFFOLD_NS << "\n";

    for (size_t pos = 0; pos < fullGenome.size(); pos += FASTA_LINE_WIDTH)
        out << fullGenome.substr(pos, FASTA_LINE_WIDTH) << "\n";

    out.close();
    std::cout << "Full pseudo-genome written to: " << outputPath << "\n";
}


// STAGE 6 — Write .visdata for visualizer


static void writeVisData(
    const std::vector<ContigAssembler::Contig>& contigs,
    const Scaffolder& scaffolder,
    VisSession& session,
    const std::string& sourcePath,
    const std::string& strategyName,
    size_t k)
{

    session.k            = k;
    session.sourceFile   = sourcePath;
    session.strategyName = strategyName;

    std::string genome;
    for (size_t i = 0; i < scaffolder.getScaffolds().size(); ++i) {
        for (const auto& entry : scaffolder.getScaffolds()[i].entries)
            genome += contigs[entry.contigIndex].sequence;
        if (i + 1 < scaffolder.getScaffolds().size())
            genome += std::string(INTER_SCAFFOLD_NS, 'N');
    }
    session.genomeSequence = std::move(genome);

    scaffolder.toVisSession(session, k - 1);

    const std::string visPath =
        "../Data/Results/assembly_k" + std::to_string(k) + "_" + strategyName + ".visdata";



    DataExporter::write(session, visPath);

    std::cout << "\nVisualization data written to: " << visPath << "\n";
    std::cout << "Launch with: Visualizer.exe " << visPath << "\n";
}

// MAIN

int main() {
    try {
        const std::string path           = "../Data/" + sequence[7];
        const size_t estimatedTotalBases = 100000;

        // Phase 1

        const std::vector<size_t> kValues = { 4 };

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

                KmerTable     kTable = loadReads(path, k, estimatedTotalBases);
                DeBruijnGraph graph  = buildGraph(kTable, k);
                auto          contigs = buildContigs(graph); // no recorder

                const KmerTable* kTablePtr =
                    (strategyName == "scored") ? &kTable : nullptr;

                Scaffolder scaffolder =
                    buildScaffolds(contigs, graph, strategy, kTablePtr);

                writeScaffoldFasta(scaffolder.getScaffolds(), contigs, k, strategyName);
            }
        }

        // Phase 2

        const size_t primaryK = 6;
        const std::string primaryStrategy = "scored";

        std::cout << "======================================\n";
        std::cout << "Primary run: k=" << primaryK
                  << ", strategy=" << primaryStrategy << "\n";
        std::cout << "======================================\n";

        // Set up session and recorder BEFORE traversal so steps are captured
        VisSession session;
        Recorder   recorder(&session);

        KmerTable     kTable  = loadReads(path, primaryK, estimatedTotalBases);
        DeBruijnGraph graph   = buildGraph(kTable, primaryK);

        // Use recorder overload so animation steps are captured
        auto contigs = buildContigs(graph, recorder);

        Scaffolder scaffolder = buildScaffolds(
            contigs, graph, ResolutionStrategy::scored(), &kTable);

        // Write full genome FASTA
        writeFullGenomeFasta(
            scaffolder.getScaffolds(), contigs,
            "Escherichia_pseudo_genome",
            "../Data/Results/full_reconstructed_genome.fna");

        // Write .visdata — populates session metadata and contig/scaffold
        // structs, then serializes the whole session to disk
        writeVisData(contigs, scaffolder, session, path, primaryStrategy, primaryK);

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}