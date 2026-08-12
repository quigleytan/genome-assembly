/*
 * data_loader.h
 * Summary:
 * - Decodes a .visdata file written by DataExporter into a VisSession.
 * - Called by gui.cpp at startup when a file is selected.
 * Important notes:
 * - NodeId hi:lo hex pairs are decoded back into __uint128_t.
 * - RUN blocks in the contig step section are expanded into individual
 *   BaseAppended steps so contig_view's animator sees a flat step list.
 * - Section order in the file must match vis_exporter's write order.
 * - Unknown section headers are skipped for forward compatibility.
 */

#ifndef DATA_LOADER_H
#define DATA_LOADER_H

#include <string>
#include "assembler/core/vis_data.h"

class VisLoader {

public:

    /**
     * @brief Reads a .visdata file and returns a fully populated VisSession.
     *
     * Throws std::runtime_error if:
     *  - The file cannot be opened
     *  - The format version is incompatible
     *  - Any required section is missing or malformed
     *
     * @param filePath Path to the .visdata file to load.
     * @return Populated VisSession ready for the visualizer.
     */
    static VisSession load(const std::string& filePath);

private:

    static constexpr int FORMAT_VERSION = 2;

    /**
     * @brief Decodes a "hi:lo" hex string back into a NodeId (__uint128_t).
     * @throws std::runtime_error if the string is not valid hi:lo format.
     */
    static NodeId decodeNodeId(const std::string& encoded);

    /**
     * @brief Parses the header block, validates the version, fills session metadata.
     * Reads until END_HEADER. Returns a struct of the counts declared in the header
     * so the loader can pre-reserve vectors before parsing each section.
     */
    struct HeaderCounts {
        size_t contigCount   = 0;
        size_t scaffoldCount = 0;
        size_t nodeCount     = 0;
        size_t edgeCount     = 0;
        size_t contigSteps   = 0;
        size_t eulerSteps    = 0;
        size_t genomeLength  = 0;
        size_t gapFillCount  = 0;
    };

    /**
     *
     * @param in
     * @param session
     * @return
     */
    static HeaderCounts parseHeader(std::istream& in, VisSession& session);

    /**
     *
     * @param in
     * @param session
     */
    static void parseGenome(std::istream& in, VisSession& session);

    /**
     * @brief Parses the BEGIN_CONTIGS ... END_CONTIGS block.
     */
    static void parseContigs(std::istream& in, VisSession& session,
                             const HeaderCounts& counts);

    /**
     * @brief Parses the BEGIN_SCAFFOLDS ... END_SCAFFOLDS block.
     */
    static void parseScaffolds(std::istream& in, VisSession& session,
                               const HeaderCounts& counts);

    /**
     * @brief Parses the BEGIN_GAPFILLS ... END_GAPFILLS block.
     */
    static void parseGapFills(std::istream& in, VisSession& session,
                              const HeaderCounts& counts);

    /**
     * @brief Parses the BEGIN_CONTIG_STEPS ... END_CONTIG_STEPS block.
     * Expands RUN blocks into individual BaseAppended steps.
     */
    static void parseContigSteps(std::istream& in, VisSession& session,
                                 const HeaderCounts& counts);

};

#endif // DATA_LOADER_H