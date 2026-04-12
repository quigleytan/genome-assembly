/*
 * vis_loader.cpp
 */

#include "assembler/graphics/vis_loader.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>

// PRIVATE HELPERS

NodeId VisLoader::decodeNodeId(const std::string& encoded) {
    // Expected format: "XXXXXXXXXXXXXXXX:YYYYYYYYYYYYYYYY" (16 hex chars each)
    size_t colon = encoded.find(':');
    if (colon == std::string::npos)
        throw std::runtime_error("vis_loader: malformed NodeId (no colon): " + encoded);

    std::string hiStr = encoded.substr(0, colon);
    std::string loStr = encoded.substr(colon + 1);

    uint64_t hi, lo;
    try {
        hi = std::stoull(hiStr, nullptr, 16);
        lo = std::stoull(loStr, nullptr, 16);
    } catch (const std::exception&) {
        throw std::runtime_error("vis_loader: malformed NodeId hex value: " + encoded);
    }

    return (static_cast<NodeId>(hi) << 64) | static_cast<NodeId>(lo);
}


// HEADER

VisLoader::HeaderCounts VisLoader::parseHeader(std::istream& in, VisSession& session) {
    HeaderCounts counts;
    std::string line;

    while (std::getline(in, line)) {
        if (line == "END_HEADER") break;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "VISDATA_VERSION") {
            int version;
            iss >> version;
            if (version != FORMAT_VERSION)
                throw std::runtime_error(
                    "vis_loader: incompatible file version " +
                    std::to_string(version) +
                    " (expected " + std::to_string(FORMAT_VERSION) + ")");

        } else if (key == "K") {
            iss >> session.k;

        } else if (key == "STRATEGY") {
            iss >> session.strategyName;

        } else if (key == "SOURCE") {
            // Source path may contain spaces - read the rest of the line
            std::getline(iss, session.sourceFile);
            // Trim the leading space left by iss >> key
            if (!session.sourceFile.empty() && session.sourceFile[0] == ' ')
                session.sourceFile = session.sourceFile.substr(1);

        } else if (key == "GENOME_LENGTH") { iss >> counts.genomeLength; }
        else if (key == "CONTIG_COUNT")   { iss >> counts.contigCount;   }
        else if (key == "SCAFFOLD_COUNT")   { iss >> counts.scaffoldCount; }
        else if (key == "NODE_COUNT")       { iss >> counts.nodeCount;     }
        else if (key == "EDGE_COUNT")       { iss >> counts.edgeCount;     }
        else if (key == "CONTIG_STEPS")     { iss >> counts.contigSteps;   }
        else if (key == "EULER_STEPS")      { iss >> counts.eulerSteps;    }
        // Unknown header keys are silently ignored for forward compatibility
    }

    return counts;
}

// SEQUENCE

void VisLoader::parseGenome(std::istream& in, VisSession& session) {
    std::string line;
    std::getline(in, line);
    if (line != "BEGIN_GENOME")
        throw std::runtime_error("vis_loader: expected BEGIN_GENOME, got: " + line);

    // Read the genome sequence - single line
    std::getline(in, session.genomeSequence);

    // Consume END_GENOME
    std::getline(in, line);
    if (line != "END_GENOME")
        throw std::runtime_error("vis_loader: expected END_GENOME, got: " + line);
}

// CONTIGS

void VisLoader::parseContigs(std::istream& in, VisSession& session,
                              const HeaderCounts& counts) {
    session.contigs.reserve(counts.contigCount);

    std::string line;
    // Consume BEGIN_CONTIGS line
    std::getline(in, line);
    if (line != "BEGIN_CONTIGS")
        throw std::runtime_error("vis_loader: expected BEGIN_CONTIGS, got: " + line);

    while (std::getline(in, line)) {
        if (line == "END_CONTIGS") break;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "CONTIG") {
            VisContig c;
            size_t index;
            size_t length;
            int circular;
            std::string startEncoded, endEncoded;

            iss >> index >> length >> circular >> c.score
                >> c.scaffoldIndex >> startEncoded >> endEncoded
                >> c.startLabel   >> c.endLabel;

            c.isCircular = (circular == 1);
            c.startNode  = decodeNodeId(startEncoded);
            c.endNode    = decodeNodeId(endEncoded);

            // Next line must be the SEQ line
            std::string seqLine;
            std::getline(in, seqLine);
            if (seqLine.substr(0, 4) != "SEQ ")
                throw std::runtime_error(
                    "vis_loader: expected SEQ after CONTIG " +
                    std::to_string(index));
            c.sequence = seqLine.substr(4);

            session.contigs.push_back(std::move(c));
        }
    }
}


// SCAFFOLDS

void VisLoader::parseScaffolds(std::istream& in, VisSession& session,
                                const HeaderCounts& counts) {
    session.scaffolds.reserve(counts.scaffoldCount);

    std::string line;
    std::getline(in, line);
    if (line != "BEGIN_SCAFFOLDS")
        throw std::runtime_error("vis_loader: expected BEGIN_SCAFFOLDS, got: " + line);

    while (std::getline(in, line)) {
        if (line == "END_SCAFFOLDS") break;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "SCAFFOLD") {
            VisScaffold s;
            size_t index, contigCount;
            int circular;
            iss >> index >> circular >> contigCount;
            s.isCircular = (circular == 1);
            s.contigIndices.reserve(contigCount);
            s.gaps.reserve(contigCount);

            // MEMBERS line
            std::string membersLine;
            std::getline(in, membersLine);
            std::istringstream mss(membersLine);
            std::string mkey;
            mss >> mkey; // Consume "MEMBERS"
            size_t idx;
            while (mss >> idx)
                s.contigIndices.push_back(idx);

            // GAPS line
            std::string gapsLine;
            std::getline(in, gapsLine);
            std::istringstream gss(gapsLine);
            std::string gkey;
            gss >> gkey; // Consume "GAPS"
            int gap;
            while (gss >> gap)
                s.gaps.push_back(gap);

            session.scaffolds.push_back(std::move(s));
        }
    }
}

// CONTIG STEPS

void VisLoader::parseContigSteps(std::istream& in, VisSession& session,
                                  const HeaderCounts& counts) {
    // Reserve based on declared step count, not RUN-compressed line count
    session.contigSteps.reserve(counts.contigSteps);

    std::string line;
    std::getline(in, line);
    if (line != "BEGIN_CONTIG_STEPS")
        throw std::runtime_error("vis_loader: expected BEGIN_CONTIG_STEPS, got: " + line);

    while (std::getline(in, line)) {
        if (line == "END_CONTIG_STEPS") break;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        iss >> key;

        if (key == "STARTED") {
            TraversalStep step;
            step.type = TraversalStep::Type::ContigStarted;
            std::string fromEncoded;
            iss >> step.contigIndex >> fromEncoded;
            step.from = decodeNodeId(fromEncoded);
            // Label is the remainder of the line
            std::getline(iss, step.sequence);
            if (!step.sequence.empty() && step.sequence[0] == ' ')
                step.sequence = step.sequence.substr(1);
            session.contigSteps.push_back(std::move(step));

        } else if (key == "RUN") {
            // Expand run-length encoded base appends back into individual steps
            size_t contigIndex, count;
            std::string bases;
            iss >> contigIndex >> count >> bases;

            for (char base : bases) {
                TraversalStep step;
                step.type        = TraversalStep::Type::BaseAppended;
                step.contigIndex = contigIndex;
                step.base        = base;
                session.contigSteps.push_back(step);
            }

        } else if (key == "FINISHED") {
            TraversalStep step;
            step.type = TraversalStep::Type::ContigFinished;
            int circular;
            iss >> step.contigIndex >> circular;
            step.base = (circular == 1) ? 'C' : 'L';
            // Sequence is the remainder of the line
            std::getline(iss, step.sequence);
            if (!step.sequence.empty() && step.sequence[0] == ' ')
                step.sequence = step.sequence.substr(1);
            session.contigSteps.push_back(std::move(step));
        }
    }
}

// PUBLIC

VisSession VisLoader::load(const std::string& filePath) {
    std::ifstream in(filePath);
    if (!in.is_open())
        throw std::runtime_error("vis_loader: could not open file: " + filePath);

    VisSession session;

    HeaderCounts counts = parseHeader(in, session);
    parseGenome(in, session);
    parseContigs(in, session, counts);
    parseScaffolds(in, session, counts);
    parseContigSteps(in, session, counts);

    return session;
}