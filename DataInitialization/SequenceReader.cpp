#include "SequenceReader.h"

#include <stdexcept>

#include "DataProcessing/KmerEncoding.h"

std::optional<DNASequence> SequenceReader::readFasta(std::istream& in) {

    std::string header;
    std::string line;
    std::string sequence;

    // Read header
    if (!std::getline(in, header))
        return std::nullopt;   // Empty file

    if (header.empty() || header[0] != '>')
        throw std::runtime_error("Invalid FASTA header");

    header = header.substr(1);  // remove '>'

    // Read sequence lines
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>')
            throw std::runtime_error("Multiple sequences not supported");
        sequence += line;
    }

    return DNASequence(header, sequence);
}

std::optional<DNASequence> SequenceReader::readFastq(std::istream& in) {

    std::string header;
    std::string sequence;
    std::string plusLine;
    std::string quality;

    // Read header
    if (!std::getline(in, header))
        return std::nullopt;   // Empty file

    if (header.empty() || header[0] != '@')
        throw std::runtime_error("Invalid FASTQ header");

    header = header.substr(1);  // remove '@'

    // Read sequence
    if (!std::getline(in, sequence))
        throw std::runtime_error("Incomplete FASTQ sequence");

    // Read '+' line
    if (!std::getline(in, plusLine) || plusLine.empty() || plusLine[0] != '+')
        throw std::runtime_error("Invalid FASTQ '+' line");

    // Read quality line
    if (!std::getline(in, quality))
        throw std::runtime_error("Missing FASTQ quality line");

    if (quality.length() != sequence.length())
        throw std::runtime_error("Sequence and quality lengths differ");

    // We ignore quality for now
    return DNASequence(header, sequence);
}

void SequenceReader::encodeAllReads(std::istream& in, size_t k, KmerTable& table) {
    size_t readCount = 0;

    while (true) {
        auto read = readFastq(in);
        if (!read) break;  // EOF

        KmerEncoding::encodeSequence(read->getSequence(), k, table);
        ++readCount;
    }

    if (readCount == 0)
        throw std::runtime_error("No reads found in file");
}