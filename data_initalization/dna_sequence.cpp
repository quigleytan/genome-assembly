#include "dna_sequence.h"

#include "exceptions/DNASequenceException.h"

// PUBLIC

DNASequence::DNASequence(const std::string& name, const std::string& inputSequence)
    : sequence_(inputSequence),
    length_(inputSequence.length()) {
    this->name_ = name;
    AnalysisResult analysis = analyzeSequence(inputSequence);
    complementSequence_ = analysis.complementSequence_;
    gcTotal_ = analysis.gcTotal_;
    gcPercent_ = analysis.gcPercent_;
}

// GETTERS

const std::string& DNASequence::getSequence() const {
    return sequence_;
}

const std::string& DNASequence::getComplement() const {
    return complementSequence_;
}

const std::string& DNASequence::getName() const {
    return name_;
}

size_t DNASequence::getLength() const {
    return length_;
}

int DNASequence::getGCCount() const {
    return gcTotal_;
}

double DNASequence::getGCPercent() const {
    return gcPercent_;
}

// OVERLOADED OPERATORS

bool operator == (const DNASequence &sequenceOne, const DNASequence &sequenceTwo) {
    return sequenceOne.getSequence() == sequenceTwo.getSequence();
}

bool operator != (const DNASequence &sequenceOne, const DNASequence &sequenceTwo) {
    return sequenceOne.getSequence() != sequenceTwo.getSequence();
}

// PRIVATE

// USED ONLY IN CONSTRUCTOR
DNASequence::AnalysisResult DNASequence::analyzeSequence(const std::string& inputSequence) {
    // Ensures no empty sequences are added
    if (inputSequence.empty()) {
        throw DNASequenceException("Sequence is empty");
    }

    std::string complement;
    int gcCount = 0;

    for (char base : inputSequence) {
        switch(base) {
            case 'A': complement += 'T'; break;
            case 'T': complement += 'A'; break;
            case 'C': complement += 'G'; gcCount++; break;
            case 'G': complement += 'C'; gcCount++; break;
                // Throws an exception if inputSequence contains an invalid base
            default: throw DNASequenceException("Invalid base: " + std::string(1, base));
        }
    }
    // G and C base percentage calculation
    double gcPercent = double(gcCount) / double(inputSequence.length());

    return {complement, gcCount, gcPercent};
}