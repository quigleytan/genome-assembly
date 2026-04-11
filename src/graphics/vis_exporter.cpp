
#include "assembler/graphics/vis_exporter.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

// PRIVATE HELPERS

std::string VisExporter::encodeNodeId(NodeId id) {
    // Split __uint128_t into high and low uint64_t halves
    uint64_t hi = static_cast<uint64_t>(id >> 64);
    uint64_t lo = static_cast<uint64_t>(id);
    std::ostringstream oss;
    oss << std::hex << std::uppercase
        << std::setw(16) << std::setfill('0') << hi
        << ':'
        << std::setw(16) << std::setfill('0') << lo;
    return oss.str();
}

void VisExporter::writeHeader(std::ostream& out, const VisSession& session) {
    out << "VISDATA_VERSION " << FORMAT_VERSION                    << '\n';
    out << "K "               << session.k                         << '\n';
    out << "STRATEGY "        << session.strategyName              << '\n';
    out << "SOURCE "          << session.sourceFile                << '\n';
    out << "CONTIG_COUNT "    << session.contigs.size()            << '\n';
    out << "SCAFFOLD_COUNT "  << session.scaffolds.size()          << '\n';
    out << "CONTIG_STEPS "    << session.contigSteps.size()        << '\n';
    out << "GENOME_LENGTH "   << session.genomeSequence.length()   << '\n';
    out << "END_HEADER\n";
}

void VisExporter::writeGenome(std::ostream& out, const VisSession& session) {
    out << "BEGIN_GENOME\n";
    out << session.genomeSequence << '\n';
    out << "END_GENOME\n";
}

void VisExporter::writeContigs(std::ostream& out, const VisSession& session) {
    out << "BEGIN_CONTIGS\n";

    for (size_t i = 0; i < session.contigs.size(); ++i) {
        const VisContig& c = session.contigs[i];

        // Line 1: numeric fields and encoded node IDs
        out << "CONTIG "
            << i                        << ' '   // index
            << c.sequence.length()      << ' '   // length
            << (c.isCircular ? 1 : 0)   << ' '   // circular flag
            << std::fixed << std::setprecision(6)
            << c.score                  << ' '   // scaffolder score
            << c.scaffoldIndex          << ' '   // scaffold membership
            << encodeNodeId(c.startNode) << ' '
            << encodeNodeId(c.endNode)   << ' '
            // Labels last — may contain any characters except newline
            << c.startLabel             << ' '
            << c.endLabel               << '\n';

        // Line 2: full sequence on its own line, prefixed so loader knows what it is
        out << "SEQ " << c.sequence << '\n';
    }

    out << "END_CONTIGS\n";
}

void VisExporter::writeScaffolds(std::ostream& out, const VisSession& session) {
    out << "BEGIN_SCAFFOLDS\n";

    for (size_t i = 0; i < session.scaffolds.size(); ++i) {
        const VisScaffold& s = session.scaffolds[i];

        out << "SCAFFOLD "
            << i                         << ' '
            << (s.isCircular ? 1 : 0)    << ' '
            << s.contigIndices.size()     << '\n';

        // Contig index list
        out << "MEMBERS";
        for (size_t idx : s.contigIndices)
            out << ' ' << idx;
        out << '\n';

        // Gap list — one entry per contig (-1 = unknown, 0 = direct overlap)
        out << "GAPS";
        for (int gap : s.gaps)
            out << ' ' << gap;
        out << '\n';
    }

    out << "END_SCAFFOLDS\n";
}

void VisExporter::writeContigSteps(std::ostream& out, const VisSession& session) {
    out << "BEGIN_CONTIG_STEPS\n";

    // Run-length encode consecutive BaseAppended steps for the same contig.
    // Instead of one line per base (expensive for long sequences), we write:
    //   RUN <contigIndex> <count> <bases>
    // which the loader expands back into individual BaseAppended steps.

    size_t i = 0;
    while (i < session.contigSteps.size()) {
        const TraversalStep& step = session.contigSteps[i];

        if (step.type == TraversalStep::Type::BaseAppended) {
            // Collect consecutive BaseAppended steps for the same contig
            size_t j = i;
            std::string bases;
            while (j < session.contigSteps.size()
                   && session.contigSteps[j].type == TraversalStep::Type::BaseAppended
                   && session.contigSteps[j].contigIndex == step.contigIndex) {
                bases += session.contigSteps[j].base;
                ++j;
            }
            out << "RUN " << step.contigIndex << ' ' << bases.length() << ' ' << bases << '\n';
            i = j;

        } else if (step.type == TraversalStep::Type::ContigStarted) {
            // sequence field holds the startLabel (see Recorder.h)
            out << "STARTED " << step.contigIndex << ' ' << encodeNodeId(step.from)
                << ' ' << step.sequence << '\n';
            ++i;

        } else if (step.type == TraversalStep::Type::ContigFinished) {
            // base field encodes circularity: 'C' = circular, 'L' = linear
            out << "FINISHED " << step.contigIndex << ' '
                << (step.base == 'C' ? 1 : 0) << ' '
                << step.sequence << '\n';
            ++i;

        } else {
            // Unknown step type — skip with a comment so the file stays valid
            out << "# unknown step type " << static_cast<int>(step.type) << '\n';
            ++i;
        }
    }

    out << "END_CONTIG_STEPS\n";
}

// PUBLIC

void VisExporter::write(const VisSession& session, const std::string& filePath) {
    std::ofstream out(filePath);
    if (!out.is_open())
        throw std::runtime_error("data_exporter: could not open file for writing: " + filePath);

    writeHeader(out, session);
    writeGenome(out, session);
    writeContigs(out, session);
    writeScaffolds(out, session);
    writeContigSteps(out, session);

    if (out.fail())
        throw std::runtime_error("data_exporter: write error on file: " + filePath);
}