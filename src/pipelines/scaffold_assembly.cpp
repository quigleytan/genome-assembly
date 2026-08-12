#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "assembler/core/kmer_table.h"
#include "assembler/core/de_bruijn_graph.h"
#include "assembler/core/vis_data.h"
#include "assembler/core/recorder.h"

#include "assembler/construction/kmer_encoding.h"
#include "assembler/construction/scaffolder.h"
#include "assembler/construction/contig_assembler.h"
#include "assembler/construction/gap_estimation.h"
#include "assembler/construction/gap_filling.h"
#include "assembler/io/sequence_reader.h"
#include "assembler/io/console_input.h"
#include "assembler/graphics/vis_exporter.h"

// CONSTANTS

static constexpr size_t UNKNOWN_GAP_NS    = 10;
static constexpr size_t FASTA_LINE_WIDTH  = 60;

// STAGE 0 - Determine table sizing

/**
 * @brief Pre-scans a FASTQ file to sum every read's sequence length, so the
 *        KmerTable can be sized correctly without asking the user to guess
 *        an estimate (the table does not support rehashing mid-run).
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
 * @brief Derives a display name for the assembled genome from the input
 *        file's base name (strips directory and extension).
 */
static std::string deriveGenomeName(const std::string& path) {
    size_t slash = path.find_last_of("/\\");
    std::string filename = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = filename.find_last_of('.');
    return (dot == std::string::npos) ? filename : filename.substr(0, dot);
}

// STAGE 1 - Load reads

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

// STAGE 2 - Build graph

static DeBruijnGraph buildGraph(const KmerTable& kTable, size_t k) {
    DeBruijnGraph graph(k, kTable.getNumItems() * 2);

    for (const auto& entry : kTable)
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);

    std::cout << "Graph built:     " << graph.getNodeCount() << " nodes, "
              << graph.getEdgeCount() << " edges\n";
    return graph;
}

// STAGE 3 - Contig traversal

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


// STAGE 4 - Scaffold


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

// STAGE 5a - Write per-scaffold FASTA

static void writeScaffoldFasta(
    const std::vector<Scaffold>& scaffolds,
    const std::vector<ContigAssembler::Contig>& contigs,
    size_t k,
    const std::string& strategyName)
{
    const std::string filename =
        "../data/results/scaffolds_k" + std::to_string(k) + "_" + strategyName + ".fna";

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


// STAGE 5b - Write full pseudo-genome FASTA

struct GapFillStats {
    size_t gapsResolved     = 0; // Bridged via local reassembly over the de Bruijn graph.
    size_t gapsPadded       = 0; // No connecting path found within the search bound - N-padded.
    size_t totalPaddedBases = 0;
};

/**
 * @brief Concatenates every scaffold's contigs into one pseudo-genome, resolving
 *        the junction between consecutive scaffolds via GapFiller.
 *
 * For each scaffold boundary, GapFiller first estimates the gap length from the
 * k-mer frequency drop at the flanking contig ends, then attempts a bounded
 * local reassembly over the graph. A resolved junction contributes its real
 * bridging sequence; an unresolved one contributes an N-run sized to the estimate.
 */
static std::string buildGenomeSequence(
    const std::vector<Scaffold>& scaffolds,
    const std::vector<ContigAssembler::Contig>& contigs,
    const GapFiller& filler,
    GapFillStats& stats)
{
    std::string genome;
    for (size_t i = 0; i < scaffolds.size(); ++i) {
        for (const auto& entry : scaffolds[i].entries)
            genome += contigs[entry.contigIndex].sequence;

        if (i + 1 < scaffolds.size()) {
            const ContigAssembler::Contig& upstream   = contigs[scaffolds[i].entries.back().contigIndex];
            const ContigAssembler::Contig& downstream = contigs[scaffolds[i + 1].entries.front().contigIndex];

            GapFiller::Result result = filler.fillGap(
                upstream.sequence, upstream.endNode,
                downstream.sequence, downstream.startNode);

            if (result.resolved) {
                genome += result.sequence;
                ++stats.gapsResolved;
            } else {
                genome += std::string(result.gapEstimate, 'N');
                ++stats.gapsPadded;
                stats.totalPaddedBases += result.gapEstimate;
            }
        }
    }
    return genome;
}

static void writeFullGenomeFasta(
    const std::vector<Scaffold>& scaffolds,
    const std::vector<ContigAssembler::Contig>& contigs,
    const std::string& genomeName,
    size_t k,
    const std::string& strategyName,
    const GapFiller& filler)
{
    const std::string outputPath =
        "../data/results/genome_k" + std::to_string(k) + "_" + strategyName + ".fna";

    std::ofstream out(outputPath);
    if (!out.is_open())
        throw std::runtime_error("Could not write output file: " + outputPath);

    GapFillStats stats;
    std::string fullGenome = buildGenomeSequence(scaffolds, contigs, filler, stats);

    out << ">" << genomeName
        << " scaffolds=" << scaffolds.size()
        << " gapsResolved=" << stats.gapsResolved
        << " gapsPadded=" << stats.gapsPadded << "\n";

    for (size_t pos = 0; pos < fullGenome.size(); pos += FASTA_LINE_WIDTH)
        out << fullGenome.substr(pos, FASTA_LINE_WIDTH) << "\n";

    out.close();
    std::cout << "Full pseudo-genome written to: " << outputPath << "\n";
    std::cout << "Gap resolution:  " << stats.gapsResolved << " resolved via local reassembly, "
              << stats.gapsPadded << " padded with N (" << stats.totalPaddedBases << " bases)\n";
}


// STAGE 6 - Write .visdata for visualizer

static void writeVisData(const std::vector<ContigAssembler::Contig>& contigs, const Scaffolder& scaffolder,
                         VisSession& session, const std::string& sourcePath, const std::string& strategyName,
                         size_t k, const GapFiller& filler) {

    session.k            = k;
    session.sourceFile   = sourcePath;
    session.strategyName = strategyName;

    GapFillStats stats;
    session.genomeSequence = buildGenomeSequence(scaffolder.getScaffolds(), contigs, filler, stats);

    scaffolder.toVisSession(session, k - 1);

    const std::string visPath =
        "../data/results/assembly_k" + std::to_string(k) + "_" + strategyName + ".visdata";

    VisExporter::write(session, visPath);
    std::cout << "Visualization data written: " << visPath << "\n";
}

// MAIN PIPELINE

static const std::vector<std::string> STRATEGY_NAMES = { "skip", "greedy", "scored" };

static ResolutionStrategy strategyByIndex(size_t index) {
    switch (index) {
        case 0:  return ResolutionStrategy::skip();
        case 1:  return ResolutionStrategy::greedy();
        default: return ResolutionStrategy::scored();
    }
}

int main() {
    try {
        const std::string path = ConsoleInput::promptFilePath(
            "FASTQ file path: ");

        const size_t totalBases = countFastqBases(path);
        std::cout << "Total bases:     " << totalBases << "\n";
        std::cout << "--------------------------------------\n";

        const size_t maxK = std::min<size_t>(KmerEncoding::MAX_K_128, totalBases);
        const std::string genomeName = deriveGenomeName(path);

        do {
            size_t k = ConsoleInput::promptSizeT("K-mer size", 2, maxK);

            size_t strategyIndex = ConsoleInput::promptChoice(
                "Resolution strategy:",
                { "skip   - never resolve ambiguous branches",
                  "greedy - take the highest length-scoring edge",
                  "scored - weighted length + k-mer frequency + overlap quality" });
            const std::string&  strategyName = STRATEGY_NAMES[strategyIndex];
            ResolutionStrategy  strategy     = strategyByIndex(strategyIndex);

            bool recordAnimation = ConsoleInput::promptYesNo(
                "Record a .visdata trace for the visualizer?");

            std::cout << "======================================\n";
            std::cout << "K = " << k << ", strategy = " << strategyName << "\n";
            std::cout << "======================================\n";

            // Set up session and recorder BEFORE traversal so steps are captured
            VisSession session;
            Recorder   recorder(&session);

            KmerTable     kTable = loadReads(path, k, totalBases);
            DeBruijnGraph graph  = buildGraph(kTable, k);

            std::vector<ContigAssembler::Contig> contigs = recordAnimation
                ? buildContigs(graph, recorder)
                : buildContigs(graph);

            const KmerTable* kTablePtr = (strategyName == "scored") ? &kTable : nullptr;

            Scaffolder scaffolder = buildScaffolds(contigs, graph, strategy, kTablePtr);

            GapEstimator gapEstimator(kTable, k);
            GapFiller    gapFiller(graph, gapEstimator);

            writeScaffoldFasta(scaffolder.getScaffolds(), contigs, k, strategyName);
            writeFullGenomeFasta(scaffolder.getScaffolds(), contigs, genomeName, k, strategyName, gapFiller);

            if (recordAnimation)
                writeVisData(contigs, scaffolder, session, path, strategyName, k, gapFiller);

        } while (ConsoleInput::promptYesNo("Run again with different settings?"));

    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}