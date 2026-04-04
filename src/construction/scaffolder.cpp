#include "scaffolding/scaffolder.h"

#include <limits>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>

#include "encoding/kmer_encoding.h"


// PRIVATE HELPER FUNCTIONS

static double mean(const std::vector<double>& data) {
    if (data.empty()) return 0.0;
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

static double stddev(const std::vector<double>& data) {
    if (data.size() < 2) return 0.0;
    double m = mean(data);
    double sum = 0.0;
    for (double x : data) sum += (x - m) * (x - m);
    return std::sqrt(sum / data.size());
}

// PRIVATE

void Scaffolder::buildConnectionMap() {
    for (size_t i = 0; i < contigs_.size(); ++i) {
        auto [endVec,   _1] = endNodeMap_.insert(contigs_[i].endNode);
        endVec.push_back(i);

        auto [startVec, _2] = startNodeMap_.insert(contigs_[i].startNode);
        startVec.push_back(i);
    }
}

double Scaffolder::computeLengthScore(const ContigTraversal::Contig& contig) const {
    return (maxContigLength_ > 0)
        ? static_cast<double>(contig.sequence.length()) / maxContigLength_
        : 0.0;
}

double Scaffolder::computeFrequencyScore(const ContigTraversal::Contig& contig) const {
    const size_t k   = graph_.getK();
    const std::string& seq = contig.sequence;

    if (seq.length() < k) return 0.0;

    size_t totalCount = 0;
    size_t kmerCount  = 0;
    std::vector<double> counts;

    NodeId kmer = KmerEncoding::encode(seq.substr(0, k));
    const size_t* count = kmerTable_->find(kmer);
    if (count) {
        totalCount += *count;
        counts.push_back(static_cast<double>(*count));
    }
    ++kmerCount;

    for (size_t i = k; i < seq.length(); ++i) {
        kmer  = KmerEncoding::roll(kmer, seq[i], k);
        count = kmerTable_->find(kmer);
        if (count) {
            totalCount += *count;
            counts.push_back(static_cast<double>(*count));
        }
        ++kmerCount;
    }

    double avgFrequency = static_cast<double>(totalCount) / kmerCount;
    double freqCap      = mean(counts) + 2.0 * stddev(counts);

    return (freqCap > 0.0) ? std::min(avgFrequency / freqCap, 1.0) : 0.0;
}

double Scaffolder::computeOverlapScore(const ContigTraversal::Contig& contig) const {
    const size_t nodeLen = graph_.getK() - 1;

    if (contig.sequence.length() < nodeLen) return 0.0;

    std::string expectedSuffix = KmerEncoding::decode(contig.endNode, nodeLen);
    std::string actualSuffix   = contig.sequence.substr(contig.sequence.length() - nodeLen);

    size_t matches = 0;
    for (size_t i = 0; i < nodeLen; ++i)
        if (expectedSuffix[i] == actualSuffix[i]) ++matches;

    return static_cast<double>(matches) / nodeLen;
}

double Scaffolder::scoreContig(size_t contigIndex) const {
    const ContigTraversal::Contig& contig = contigs_[contigIndex];

    double effectiveLengthWeight    = strategy_.weights.lengthWeight;
    double effectiveFrequencyWeight = strategy_.weights.kmerFrequencyWeight;

    if (kmerTable_ == nullptr && effectiveFrequencyWeight > 0.0) {
        effectiveLengthWeight    += effectiveFrequencyWeight;
        effectiveFrequencyWeight  = 0.0;
    }

    return (effectiveLengthWeight * computeLengthScore(contig))
         + (effectiveFrequencyWeight > 0.0
                ? effectiveFrequencyWeight * computeFrequencyScore(contig) : 0.0)
         + (strategy_.weights.overlapQualityWeight > 0.0
                ? strategy_.weights.overlapQualityWeight * computeOverlapScore(contig) : 0.0);
}

size_t Scaffolder::resolveNext(NodeId boundaryNode,
                                     const std::vector<bool>& visited) const {
    const std::vector<size_t>* candidates = startNodeMap_.find(boundaryNode);

    if (!candidates || candidates->empty())
        return std::numeric_limits<size_t>::max();

    if (strategy_.skipAmbiguous && candidates->size() > 1)
        return std::numeric_limits<size_t>::max();

    size_t bestContig = std::numeric_limits<size_t>::max();
    double bestScore  = -1.0;

    for (size_t index : *candidates) {
        if (!visited[index]) {
            double score = scoreContig(index);
            if (score > bestScore) {
                bestScore  = score;
                bestContig = index;
            }
        }
    }
    return bestContig;
}

Scaffold Scaffolder::walkScaffold(size_t startIndex,
                                        std::vector<bool>& visited) const {
    Scaffold scaffold;
    size_t contigIndex = startIndex;

    while (true) {
        visited[contigIndex] = true;

        ScaffoldEntry entry;
        entry.contigIndex = contigIndex;
        entry.gapAfter    = ScaffoldEntry::DIRECT_OVERLAP;
        entry.score       = scoreContig(contigIndex); // ← now populated

        size_t next = resolveNext(contigs_[contigIndex].endNode, visited);

        if (next == std::numeric_limits<size_t>::max()) {
            entry.gapAfter = ScaffoldEntry::UNKNOWN_GAP;
            scaffold.entries.push_back(entry);
            break;
        }

        scaffold.entries.push_back(entry);
        contigIndex = next;
    }

    if (scaffold.entries.size() > 1) {
        NodeId lastEnd    = contigs_[scaffold.entries.back().contigIndex].endNode;
        NodeId firstStart = contigs_[scaffold.entries.front().contigIndex].startNode;
        scaffold.isCircular = (lastEnd == firstStart);
    }

    return scaffold;
}

bool Scaffolder::isScaffoldStart(size_t contigIndex) const {
    const auto* predecessors = endNodeMap_.find(contigs_[contigIndex].startNode);
    return predecessors == nullptr || predecessors->empty();
}

// PUBLIC

Scaffolder::Scaffolder(const std::vector<ContigTraversal::Contig>& contigs,
                                   const DeBruijnGraph& graph,
                                   ResolutionStrategy strategy,
                                   const KmerTable* kmerTable)
    : contigs_(contigs),
      graph_(graph),
      kmerTable_(kmerTable),
      strategy_(strategy),
      endNodeMap_(contigs.size() * 2),
      startNodeMap_(contigs.size() * 2) {}

void Scaffolder::buildScaffolds() {
    scaffolds_.clear();
    buildConnectionMap();

    maxContigLength_ = 0;
    for (const auto& c : contigs_)
        maxContigLength_ = std::max(maxContigLength_, c.sequence.length());

    std::vector<bool> visited(contigs_.size(), false);

    // Stage 1: linear entry points
    for (size_t i = 0; i < contigs_.size(); ++i) {
        if (!visited[i] && isScaffoldStart(i))
            scaffolds_.push_back(walkScaffold(i, visited));
    }

    // Stage 2: isolated cycles / circular scaffolds
    for (size_t i = 0; i < contigs_.size(); ++i) {
        if (!visited[i])
            scaffolds_.push_back(walkScaffold(i, visited));
    }
}

const std::vector<Scaffold>& Scaffolder::getScaffolds() const {
    return scaffolds_;
}

void Scaffolder::toVisSession(VisSession& session, size_t nodeLen) const {
    // Build a scaffold-index lookup so we can stamp each VisContig with which scaffold it belongs to in a single pass.
    std::vector<int> contigToScaffold(contigs_.size(), -1);
    for (size_t si = 0; si < scaffolds_.size(); ++si) {
        for (const auto& entry : scaffolds_[si].entries)
            contigToScaffold[entry.contigIndex] = static_cast<int>(si);
    }

    session.contigs.clear();
    session.contigs.reserve(contigs_.size());

    for (size_t i = 0; i < contigs_.size(); ++i) {
        const ContigTraversal::Contig& src = contigs_[i];

        VisContig vc;
        vc.sequence      = src.sequence;
        vc.startNode     = src.startNode;
        vc.endNode       = src.endNode;
        vc.startLabel    = KmerEncoding::decode(src.startNode, nodeLen);
        vc.endLabel      = KmerEncoding::decode(src.endNode,   nodeLen);
        vc.isCircular    = src.isCircular;
        vc.scaffoldIndex = contigToScaffold[i];

        // Pull the score from the corresponding ScaffoldEntry if available
        vc.score = 0.0;
        if (vc.scaffoldIndex >= 0) {
            const Scaffold& sc = scaffolds_[vc.scaffoldIndex];
            for (const auto& entry : sc.entries) {
                if (entry.contigIndex == i) {
                    vc.score = entry.score;
                    break;
                }
            }
        }

        session.contigs.push_back(std::move(vc));
    }

    // Populate VisScaffold list
    session.scaffolds.clear();
    session.scaffolds.reserve(scaffolds_.size());

    for (const auto& scaffold : scaffolds_) {
        VisScaffold vs;
        vs.isCircular = scaffold.isCircular;
        vs.contigIndices.reserve(scaffold.entries.size());
        vs.gaps.reserve(scaffold.entries.size());

        for (const auto& entry : scaffold.entries) {
            vs.contigIndices.push_back(entry.contigIndex);
            vs.gaps.push_back(entry.gapAfter);
        }

        session.scaffolds.push_back(std::move(vs));
    }
}

void Scaffolder::printStats() const {
    if (scaffolds_.empty()) {
        std::cout << "No scaffolds found\n";
        return;
    }

    size_t totalLength    = 0;
    size_t circularCount  = 0;
    size_t unresolvedGaps = 0;
    std::vector<size_t> lengths;
    lengths.reserve(scaffolds_.size());

    for (const auto& scaffold : scaffolds_) {
        size_t scaffoldLength = 0;
        for (const auto& entry : scaffold.entries) {
            scaffoldLength += contigs_[entry.contigIndex].sequence.length();
            if (entry.gapAfter == ScaffoldEntry::UNKNOWN_GAP)
                ++unresolvedGaps;
        }
        lengths.push_back(scaffoldLength);
        totalLength += scaffoldLength;
        if (scaffold.isCircular) ++circularCount;
    }

    std::sort(lengths.rbegin(), lengths.rend());
    size_t half        = (totalLength + 1) / 2;
    size_t accumulated = 0;
    size_t n50         = 0;
    for (size_t len : lengths) {
        accumulated += len;
        if (accumulated >= half) { n50 = len; break; }
    }

    std::string strategyName = strategy_.skipAmbiguous ? "skip"
                             : (strategy_.weights.kmerFrequencyWeight > 0 ||
                                strategy_.weights.overlapQualityWeight > 0)
                               ? "scored" : "greedy";

    std::cout << "--------------------------------------\n";
    std::cout << "Strategy:           " << strategyName      << "\n";
    std::cout << "Total scaffolds:    " << scaffolds_.size() << "\n";
    std::cout << "Circular scaffolds: " << circularCount     << "\n";
    std::cout << "Total bases:        " << totalLength       << "\n";
    std::cout << "Shortest scaffold:  " << lengths.back()    << " bases\n";
    std::cout << "Longest scaffold:   " << lengths.front()   << " bases\n";
    std::cout << "Unresolved gaps:    " << unresolvedGaps    << "\n";
    std::cout << "N50:                " << n50               << " bases\n";
    std::cout << "--------------------------------------\n";
}