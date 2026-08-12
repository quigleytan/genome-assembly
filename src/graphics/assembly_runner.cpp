/*
 * assembly_runner.cpp
 * Summary:
 * - Implements runAssembly() - the full pipeline from FASTQ to .visdata,
 *   running on a background thread launched by gui.cpp.
 * - Updates AssemblyProgress at each stage boundary so the UI can display
 *   a live progress bar.
 * - This file mirrors the pipeline logic in ScaffoldAssembly.cpp but is
 *   structured as a single callable function rather than a main().
 */

#include "assembler/graphics/assembly_runner.h"

#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <filesystem>

#include "assembler/io/sequence_reader.h"
#include "assembler/core/kmer_table.h"
#include "assembler/core/de_bruijn_graph.h"
#include "assembler/core/vis_data.h"
#include "assembler/core/recorder.h"
#include "assembler/construction/contig_assembler.h"
#include "assembler/construction/scaffolder.h"
#include "assembler/construction/gap_estimation.h"
#include "assembler/construction/gap_filling.h"
#include "assembler/graphics/vis_exporter.h"


// CONSTANTS
static constexpr size_t FASTA_LINE_WIDTH  = 60;


// PIPELINE STAGES
// Each stage updates progress before and after
// so the UI shows meaningful state transitions.

static KmerTable loadReads(const std::string& path,
                            size_t k,
                            size_t totalBases,
                            AssemblyProgress& progress) {
    progress.setMessage("Loading reads from " + path);
    progress.progress = 0.1f;

    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Could not open file: " + path);

    KmerTable kTable(totalBases, k);
    SequenceReader::encodeAllReads(file, k, kTable);
    file.close();

    progress.setMessage("Reads loaded: " + std::to_string(kTable.getNumItems())
                        + " unique k-mers");
    progress.progress = 1.0f;
    return kTable;
}

static DeBruijnGraph buildGraph(const KmerTable& kTable,
                                 size_t k,
                                 AssemblyProgress& progress) {
    progress.setMessage("Building De Bruijn graph...");
    progress.progress = 0.1f;

    DeBruijnGraph graph(k, kTable.getNumItems() * 2);
    for (const auto& entry : kTable)
        for (size_t i = 0; i < entry.value; ++i)
            graph.addKmer(entry.key);

    progress.setMessage("Graph built: " + std::to_string(graph.getNodeCount())
                        + " nodes, " + std::to_string(graph.getEdgeCount()) + " edges");
    progress.progress = 1.0f;
    return graph;
}

static std::vector<ContigAssembler::Contig> buildContigs(DeBruijnGraph& graph,
                                                          Recorder& recorder,
                                                          AssemblyProgress& progress) {
    progress.setMessage("Assembling contigs...");
    progress.progress = 0.1f;

    ContigAssembler ct(graph, &recorder);
    ct.computeContigs();

    const auto& contigs = ct.getContigs();
    progress.setMessage("Contigs assembled: " + std::to_string(contigs.size()));
    progress.progress = 1.0f;
    return contigs;
}

static Scaffolder buildScaffolds(
    const std::vector<ContigAssembler::Contig>& contigs,
    const DeBruijnGraph& graph,
    const std::string& strategyName,
    const KmerTable& kTable,
    AssemblyProgress& progress)
{
    progress.setMessage("Scaffolding contigs (" + strategyName + ")...");
    progress.progress = 0.1f;

    ResolutionStrategy strategy;
    if      (strategyName == "skip")   strategy = ResolutionStrategy::skip();
    else if (strategyName == "greedy") strategy = ResolutionStrategy::greedy();
    else                               strategy = ResolutionStrategy::scored();

    const KmerTable* kTablePtr = (strategyName == "scored") ? &kTable : nullptr;

    Scaffolder scaffolder(contigs, graph, strategy, kTablePtr);
    scaffolder.buildScaffolds();

    progress.setMessage("Scaffolding complete: "
                        + std::to_string(scaffolder.getScaffolds().size())
                        + " scaffolds");
    progress.progress = 1.0f;
    return scaffolder;
}

// STAGE 5 - Resolve gaps between scaffolds
//
// For each junction between consecutive scaffolds, estimates the gap length
// from the k-mer frequency drop at the flanking contig ends, then attempts a
// bounded local reassembly over the graph. A resolved junction contributes
// its real bridging sequence to genome; an unresolved one contributes an
// N-run sized to the estimate. Either way, a VisGapFill record is appended
// so the visualizer can show what happened at that junction.

struct GapFillStats {
    size_t gapsResolved     = 0;
    size_t gapsPadded       = 0;
    size_t totalPaddedBases = 0;
};

static void resolveGaps(
    const Scaffolder& scaffolder,
    const std::vector<ContigAssembler::Contig>& contigs,
    const GapFiller& filler,
    VisSession& session,
    AssemblyProgress& progress)
{
    progress.setMessage("Resolving gaps between scaffolds...");
    progress.progress = 0.1f;

    const auto& scaffolds = scaffolder.getScaffolds();

    std::string genome;
    session.gapFills.clear();
    if (!scaffolds.empty())
        session.gapFills.reserve(scaffolds.size() - 1);

    GapFillStats stats;

    for (size_t i = 0; i < scaffolds.size(); ++i) {
        for (const auto& entry : scaffolds[i].entries)
            genome += contigs[entry.contigIndex].sequence;

        if (i + 1 < scaffolds.size()) {
            const ContigAssembler::Contig& upstream   = contigs[scaffolds[i].entries.back().contigIndex];
            const ContigAssembler::Contig& downstream = contigs[scaffolds[i + 1].entries.front().contigIndex];

            GapFiller::Result result = filler.fillGap(
                upstream.sequence, upstream.endNode,
                downstream.sequence, downstream.startNode);

            VisGapFill gf;
            gf.resolved     = result.resolved;
            gf.estimatedGap = result.gapEstimate;

            if (result.resolved) {
                genome += result.sequence;
                gf.filledLength = result.sequence.length();
                ++stats.gapsResolved;
            } else {
                genome += std::string(result.gapEstimate, 'N');
                gf.filledLength = result.gapEstimate;
                ++stats.gapsPadded;
                stats.totalPaddedBases += result.gapEstimate;
            }

            session.gapFills.push_back(gf);
        }

        progress.progress = scaffolds.size() > 1
            ? static_cast<float>(i + 1) / static_cast<float>(scaffolds.size())
            : 1.0f;
    }

    session.genomeSequence = std::move(genome);

    progress.setMessage("Gaps resolved: " + std::to_string(stats.gapsResolved)
                        + " via local reassembly, " + std::to_string(stats.gapsPadded)
                        + " padded with N (" + std::to_string(stats.totalPaddedBases) + " bases)");
    progress.progress = 1.0f;
}

// STAGE 6 - Write outputs

static std::string writeOutputs(
    const Scaffolder& scaffolder,
    VisSession& session,
    const std::string& inputPath,
    const std::string& outputDir,
    const std::string& strategyName,
    size_t k,
    AssemblyProgress& progress)
{

    std::filesystem::path outDir(outputDir);
    std::error_code ec;
    std::filesystem::create_directories(outDir, ec);

    progress.progress = 0.2f;
    session.k            = k;
    session.sourceFile   = inputPath;
    session.strategyName = strategyName;

    scaffolder.toVisSession(session, k - 1);

    progress.progress = 0.5f;

    const std::string fastaPath = outputDir + "/full_reconstructed_genome.fna";
    {
        std::ofstream out(fastaPath);
        if (!out.is_open())
            throw std::runtime_error("Could not write FASTA: " + fastaPath);

        out << ">assembled_genome scaffolds=" << scaffolder.getScaffolds().size()
            << " k=" << k << " strategy=" << strategyName << "\n";

        const std::string& gs = session.genomeSequence;
        for (size_t pos = 0; pos < gs.size(); pos += FASTA_LINE_WIDTH)
            out << gs.substr(pos, FASTA_LINE_WIDTH) << "\n";
    }

    progress.progress = 0.8f;

    const std::string visPath = outputDir + "/assembly_k" + std::to_string(k)
                              + "_" + strategyName + ".visdata";

    VisExporter::write(session, visPath);

    progress.setMessage("Output written to " + visPath);
    progress.progress = 1.0f;

    return visPath;
}

// PUBLIC ENTRY POINT

void runAssembly(AssemblyConfig config, AssemblyProgress& progress) {
    try {
        // Stage 1 - Load reads
        progress.stage = AssemblyProgress::Stage::LoadingReads;
        KmerTable kTable = loadReads(config.inputPath, config.k,
                                     config.estimatedTotalBases, progress);

        // Stage 2 - Build graph
        progress.stage = AssemblyProgress::Stage::BuildingGraph;
        DeBruijnGraph graph = buildGraph(kTable, config.k, progress);

        // Stage 3 - Assemble contigs (with recorder for animation)
        progress.stage = AssemblyProgress::Stage::AssemblingContigs;
        VisSession session;
        Recorder   recorder(&session);
        auto contigs = buildContigs(graph, recorder, progress);

        // Stage 4 - Scaffold
        progress.stage = AssemblyProgress::Stage::Scaffolding;
        Scaffolder scaffolder = buildScaffolds(
            contigs, graph, config.strategy, kTable, progress);

        // Stage 5 - Resolve gaps between scaffolds
        progress.stage = AssemblyProgress::Stage::ResolvingGaps;
        GapEstimator gapEstimator(kTable, config.k);
        GapFiller    gapFiller(graph, gapEstimator);
        resolveGaps(scaffolder, contigs, gapFiller, session, progress);

        // Stage 6 - Write outputs
        progress.stage = AssemblyProgress::Stage::WritingOutput;
        std::string visPath = writeOutputs(
            scaffolder, session,
            config.inputPath, config.outputDir,
            config.strategy, config.k, progress);

        // Done
        progress.setOutputPath(visPath);
        progress.stage = AssemblyProgress::Stage::Complete;

    } catch (const std::exception& e) {
        progress.setError(e.what());
        progress.stage = AssemblyProgress::Stage::Failed;
    }
}