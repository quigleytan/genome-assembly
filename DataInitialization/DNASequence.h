/*
 * DNASequence.h
 * Summary:
 * - Represents a DNA sequence and provides analysis of its composition.
 * Features:
 * - Stores a DNA sequence (immutable once constructed).
 * - Computes the complement.
 * - Computes GC count and GC percent.
 * - Throws DNASequenceException for an invalid sequence input.
 * Additional:
 * - Internally uses a helper struct AnalysisResult for single-pass computation.
 * - All internal methods and structs are private.
 */

#ifndef DNA_SEQUENCE_H
#define DNA_SEQUENCE_H

#include <string>

class DNASequence {

public:

    /**
     * @brief Constructor for DNASequence.
     *
     * Takes the name and string of bases to create a DNASequence object.
     * The sequence is analyzed for its complement, GC count percentage.
     *
     * @param inputSequence String of char's to be analyzed and stored.
     * @param name          Name of the sequence derived from a data file.
     */
    DNASequence(const std::string& name, const std::string& inputSequence);

    /**
     * @brief Returns the DNA sequence string.
     * @return sequence The DNA sequence string (immutable after construction).
     */
    [[nodiscard]] const std::string& getSequence() const;

    /**
     * @brief Returns the complement of the sequence.
     * @return The complement of the sequence.
     */
    [[nodiscard]] const std::string& getComplement() const;

    /**
     * @brief Returns the associated name with the sequence.
     * @return Name of the sequence derived from a data file.
     */
    [[nodiscard]] const std::string& getName() const;

    /**
     * @brief Returns the length of the sequence.
     * @return size_t length of the sequence (number of bases).
     */
    [[nodiscard]] size_t getLength() const;

    /**
     * @brief Returns the number of instances of G's and C's in the string.
     * @return gcTotal The total count of guanine and cytosine in the sequence.
     */
    [[nodiscard]] int getGCCount() const;

    /**
     * @brief Returns the GC content of the string.
     * @return gcPercent The percentage of guanine and cytosine in the sequence.
     */
    [[nodiscard]] double getGCPercent() const;

    // Overloaded comparisons
    friend bool operator == (const DNASequence &sequenceOne, const DNASequence &sequenceTwo);

    friend bool operator != (const DNASequence &sequenceOne, const DNASequence &sequenceTwo);

private:

    std::string name_;
    std::string sequence_;
    std::string complementSequence_;
    int gcTotal_;
    double gcPercent_;
    size_t length_; // # of bases

    // Struct to allow for a single loop data analysis without modifying
    // members outside the constructor.
    struct AnalysisResult {
        std::string complementSequence_; // Complement sequence
        int gcTotal_;                    // Number of G/C bases
        double gcPercent_;               // Percentage of G/C
    };

    /**
     * @brief Computes basic sequence data
     *
     * Creates the complement, G/C counts, and percentage in a single pass through the input sequence.
     * Used only in the constructor.
     *
     * @param inputSequence Sequence passed in through the constructor.
     * @return AnalysisResult Struct containing the complement, gcTotal, and gcPercent.
     */
    static AnalysisResult analyzeSequence(const std::string& inputSequence);

};

#endif //DNA_SEQUENCE_H