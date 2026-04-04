/*
 * SequenceReader.h
 * Created by Tanner Quigley on 2/15/2026
 * Summary:
 * - Pulls information from FASTA and FASTQ files to create DNASequence objects
 * - Uses std::optional to handle empty files gracefully
 * - Throws std::runtime_error for various format errors in the input files
 * Important notes:
 * - Returns optionals, need to extract the object from the optional
 * TODO items for future work;
 * - Potentially add more support for other files
 */

#ifndef SEQUENCE_READER_H
#define SEQUENCE_READER_H

#include <optional>
#include <istream>
#include <string>

#include "dna_sequence.h"
#include "encoding/kmer_table.h"

class SequenceReader {

public:

    /**
     * @brief Reads a FASTA format sequence into optional DNASequence.
     * @param in Data file containing sequence information
     * @return Optional DNASequence object if successful, null optional if the file is empty
     */
    static std::optional<DNASequence> readFasta(std::istream& in);

    /**
     * @brief Reads a FASTQ format sequence into optional DNASequence.
     * @param in Data file containing sequence information
     * @return Optional DNASequence object if successful, null optional if the file is empty
     */
    static std::optional<DNASequence> readFastq(std::istream& in);

    /**
     * @breif Reads all sequences from a FASTQ stream and encodes k-mers into the provided table.
     * @param in    Input stream to be read from.
     * @param k     Length of k-mers to be encoded.
     * @param table Table used to store k-mers.
     */
    static void encodeAllReads(std::istream& in, size_t k, KmerTable& table);

};

#endif //SEQUENCE_READER_H