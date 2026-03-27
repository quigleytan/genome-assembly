#include "ContigScaffolder.h"

#include <limits>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "DataProcessing/KmerEncoding.h"

// PRIVATE

double mean(const std::vector<double>& data) {
    if (data.empty()) return 0.0;
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

double stddev(const std::vector<double>& data) {
    if (data.size() < 2) return 0.0;
    double m = mean(data);
    double sum = 0.0;
    for (double x : data) sum += (x - m) * (x - m);
    return std::sqrt(sum / data.size());
}

void ContigScaffolder::buildConnectionMap() {
    for (size_t i = 0; i < contigs_.size(); ++i) {
        // Insert index i into endNodeMap_ under this contig's endNode
        auto [endVec, _1] = endNodeMap_.insert(contigs_[i].endNode);
        endVec.push_back(i);

        // Insert index i into startNodeMap_ under this contig's startNode
        auto [startVec, _2] = startNodeMap_.insert(contigs_[i].startNode);
        startVec.push_back(i);
    }
}

double ContigScaffolder::computeLengthScore(const ContigTraversal::Contig& contig) const {
    // Metric 1: Length
    return (maxContigLength_ > 0)
        ? static_cast<double>(contig.sequence.length()) / maxContigLength_: 0.0;
}


double ContigScaffolder::computeFrequencyScore(const ContigTraversal::Contig& contig) const {
    const size_t k = graph_.getK();
    const std::string& seq = contig.sequence;

    if (seq.length() >= k) {
        size_t totalCount = 0;
        size_t kmerCount  = 0;
        std::vector<double> counts;

        // encode each kmer in the contig and look up its count
        NodeId kmer = KmerEncoding::encode(seq.substr(0, k));
        const size_t* count = kmerTable_->find(kmer);
        if (count) {
            totalCount += *count;
            counts.push_back(static_cast<double>(*count)); // missing in your version
        }
        ++kmerCount;

        for (size_t i = k; i < seq.length(); ++i) {
            kmer = KmerEncoding::roll(kmer, seq[i], k);
            count = kmerTable_->find(kmer);
            if (count) {
                totalCount += *count;
                counts.push_back(static_cast<double>(*count));
            }
            ++kmerCount;
        }

        double avgFrequency = static_cast<double>(totalCount) / kmerCount;
        double freqCap = mean(counts) + 2.0 * stddev(counts);

        return (freqCap > 0.0) ? std::min(avgFrequency / freqCap, 1.0) : 0.0;
    }

    return 0.0;
}

double ContigScaffolder::computeOverlapScore(const ContigTraversal::Contig& contig) const {
    const size_t nodeLen = graph_.getK() - 1;

    if (contig.sequence.length() < nodeLen)
        return 0.0;

    std::string expectedSuffix = KmerEncoding::decode(contig.endNode, nodeLen);
    std::string actualSuffix   = contig.sequence.substr(contig.sequence.length() - nodeLen);

    size_t matches = 0;
    for (size_t i = 0; i < nodeLen; ++i)
        if (expectedSuffix[i] == actualSuffix[i]) ++matches;

    return static_cast<double>(matches) / nodeLen;
}

double ContigScaffolder::scoreContig(size_t contigIndex) const {
    const ContigTraversal::Contig& contig = contigs_[contigIndex];

    double effectiveLengthWeight    = strategy_.weights.lengthWeight;
    double effectiveFrequencyWeight = strategy_.weights.kmerFrequencyWeight;

    if (kmerTable_ == nullptr && effectiveFrequencyWeight > 0.0) {
        effectiveLengthWeight    += effectiveFrequencyWeight;
        effectiveFrequencyWeight  = 0.0;
    }

    return (effectiveLengthWeight * computeLengthScore(contig)) + (effectiveFrequencyWeight > 0.0
    ? effectiveFrequencyWeight * computeFrequencyScore(contig) : 0.0) + (strategy_.weights.overlapQualityWeight > 0.0
    ? strategy_.weights.overlapQualityWeight * computeOverlapScore(contig) : 0.0);
}

size_t ContigScaffolder::resolveNext(NodeId boundaryNode, const std::vector<bool>& visited) const {
    const std::vector<size_t>* candidates = startNodeMap_.find(boundaryNode);

    if (!candidates || candidates->empty())
        return std::numeric_limits<size_t>::max(); // Dead end reached

    if (strategy_.skipAmbiguous && candidates->size() > 1)
        return std::numeric_limits<size_t>::max(); // Path is ambiguous, strategy says stop

    size_t returnContig = std::numeric_limits<size_t>::max();
    double highestScore  = -1.0;

    for (size_t i = 0; i < candidates->size(); ++i) {
        size_t index = candidates->at(i);
        if (!visited[index]) {
            double score = scoreContig(index);
            if (score > highestScore) {
                highestScore = score;
                returnContig = index;
            }
        }
    }
    return returnContig;
}

Scaffold ContigScaffolder::walkScaffold(size_t startIndex, std::vector<bool>& visited) const {
    Scaffold scaffold;
    size_t contigIndex = startIndex;

    while (true) {

        visited[contigIndex] = true;

        ScaffoldEntry entry;
        entry.contigIndex = contigIndex;
        entry.gapAfter    = ScaffoldEntry::DIRECT_OVERLAP; // default, overridden on exit

        size_t next = resolveNext(contigs_[contigIndex].endNode, visited);

        if (next == std::numeric_limits<size_t>::max()) {
            // Dead end or ambiguous branch — terminate scaffold
            entry.gapAfter = ScaffoldEntry::UNKNOWN_GAP;
            scaffold.entries.push_back(entry);
            break;
        }

        scaffold.entries.push_back(entry);
        contigIndex = next;
    }

    if (scaffold.entries.size() > 1) {
        NodeId lastEnd   = contigs_[scaffold.entries.back().contigIndex].endNode;
        NodeId firstStart = contigs_[scaffold.entries.front().contigIndex].startNode;
        scaffold.isCircular = (lastEnd == firstStart);
    }

    return scaffold;
}

bool ContigScaffolder::isScaffoldStart(size_t contigIndex) const {
    const auto* predecessors = endNodeMap_.find(contigs_[contigIndex].startNode);
    return predecessors == nullptr || predecessors->empty();
}

// PUBLIC

ContigScaffolder::ContigScaffolder(const std::vector<ContigTraversal::Contig>& contigs,
                                   const DeBruijnGraph& graph,
                                   ResolutionStrategy strategy,
                                   const KmerTable* kmerTable)
    : contigs_(contigs),
      graph_(graph),
      kmerTable_(kmerTable),
      strategy_(strategy),
      endNodeMap_(contigs.size() * 2),
      startNodeMap_(contigs.size() * 2) {}

void ContigScaffolder::buildScaffolds() {
    scaffolds_.clear();
    buildConnectionMap();

    for (const auto& c : contigs_)
        maxContigLength_ = std::max(maxContigLength_, c.sequence.length());

    std::vector visited(contigs_.size(), false);

    // Phase 1: walk from all linear entry points (contigs with no predecessor)
    for (size_t i = 0; i < contigs_.size(); ++i) {
        if (!visited[i] && isScaffoldStart(i)) {
            scaffolds_.push_back(walkScaffold(i, visited));
        }
    }

    // Phase 2: handle remaining unvisited contigs (isolated cycles or
    // circular scaffolds with no external entry point)
    for (size_t i = 0; i < contigs_.size(); ++i) {
        if (!visited[i]) {
            scaffolds_.push_back(walkScaffold(i, visited));
        }
    }
}

const std::vector<Scaffold> &ContigScaffolder::getScaffolds() const {
    return scaffolds_;
}

void ContigScaffolder::printStats() const {
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
        // Scaffold length = sum of all contig lengths in this scaffold
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

    // N50 calculation
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
    std::cout << "Strategy:           " << strategyName       << "\n";
    std::cout << "Total scaffolds:    " << scaffolds_.size()  << "\n";
    std::cout << "Circular scaffolds: " << circularCount      << "\n";
    std::cout << "Total bases:        " << totalLength        << "\n";
    std::cout << "Shortest scaffold:  " << lengths.back()     << " bases\n";
    std::cout << "Longest scaffold:   " << lengths.front()    << " bases\n";
    std::cout << "Unresolved gaps:    " << unresolvedGaps     << "\n";
    std::cout << "N50:                " << n50                << " bases\n";
    std::cout << "--------------------------------------\n";
}
