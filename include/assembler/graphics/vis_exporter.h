/*
 * data_exporter.h
 * Summary:
 * - Serializes a VisSession to a plain-text .visdata file.
 * - Called at the end of a pipeline run when visualization is requested.
 * - The .visdata format is line-oriented and human-readable, making it
 *   easy to inspect and debug without a special tool.
 * Important notes:
 * - NodeId (__uint128_t) is written as two hex uint64_t values
 *   separated by ':' since __uint128_t has no standard stream support.
 * - Strings that may contain spaces (sequences, labels, file paths) are
 *   always written last on their line so the loader can use getline.
 * - The file is self-describing: a VERSION header lets the loader reject
 *   files written by an incompatible older exporter.
 */

#ifndef DATA_EXPORTER_H
#define DATA_EXPORTER_H

#include <string>

#include "assembler/core/vis_data.h"

class DataExporter {

public:

    /**
     * @brief Writes a VisSession to a .visdata file at the given path.
     *
     * Creates or overwrites the file. Throws std::runtime_error if the
     * file cannot be opened for writing.
     *
     * @param session  The fully populated VisSession to serialize.
     * @param filePath Output file path (e.g. "../Data/Results/run.visdata").
     */
    static void write(const VisSession& session, const std::string& filePath);

private:

    // Current file format version. Increment if the format changes
    // in a way that breaks backward compatibility with DataLoader.
    static constexpr int FORMAT_VERSION = 1;

    /**
     * @brief Encodes a NodeId as "hi:lo" hex string.
     * Splits __uint128_t into two uint64_t halves for portable serialization.
     */
    static std::string encodeNodeId(NodeId id);

    /**
     * @brief Writes the SESSION header block.
     * Contains version, k, sourceFile, and strategyName.
     */
    static void writeHeader(std::ostream& out, const VisSession& session);

    static void writeGenome(std::ostream& out, const VisSession& session);
    /**
     * @brief Writes all VisContig entries.
     * Format per contig:
     *   CONTIG <index> <length> <circular> <score> <scaffoldIndex>
     *           <startNodeHi:Lo> <endNodeHi:Lo> <startLabel> <endLabel>
     *   SEQ <sequence>
     */
    static void writeContigs(std::ostream& out, const VisSession& session);

    /**
     * @brief Writes all VisScaffold entries.
     * Format per scaffold:
     *   SCAFFOLD <index> <circular> <contigCount>
     *   MEMBERS <idx0> <idx1> ...
     *   GAPS    <gap0> <gap1> ...
     */
    static void writeScaffolds(std::ostream& out, const VisSession& session);

    /**
     * @brief Writes the contig animation step list.
     * Each step is one line. BaseAppended steps are compacted into
     * run-length encoded RUN blocks to keep file size manageable
     * for long sequences.
     */
    static void writeContigSteps(std::ostream& out, const VisSession& session);

};

#endif // DATA_EXPORTER_H