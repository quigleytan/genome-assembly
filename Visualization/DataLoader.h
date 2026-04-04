/*
 * DataLoader.h
 * Summary:
 * - Deserializes a .visdata file written by DataExporter into a VisSession.
 * - Called by VisualizerApp at startup when a file is selected.
 * - Validates the format version before parsing and throws descriptively
 *   on any malformed line so the user knows exactly what went wrong.
 * Important notes:
 * - NodeId hi:lo hex pairs are decoded back into __uint128_t.
 * - RUN blocks in the contig step section are expanded into individual
 *   BaseAppended steps so ContigView's animator sees a flat step list.
 * - Section order in the file must match DataExporter's write order.
 *   Unknown section headers are skipped for forward compatibility.
 */

#ifndef DATA_LOADER_H
#define DATA_LOADER_H

#include <string>
#include "VisData.h"

class DataLoader {

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

    static constexpr int FORMAT_VERSION = 1;

    /**
     * @brief Decodes a "hi:lo" hex string back into a NodeId (__uint128_t).
     * @throws std::runtime_error if the string is not valid hi:lo format.
     */
    static NodeId decodeNodeId(const std::string& encoded);

    /**
     * @brief Parses the header block, validates version, fills session metadata.
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
    };
    static HeaderCounts parseHeader(std::istream& in, VisSession& session);

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
     * @brief Parses the BEGIN_NODES ... END_NODES block (phase 2).
     * Produces an empty node list if the block is present but empty.
     */
    static void parseNodes(std::istream& in, VisSession& session,
                           const HeaderCounts& counts);

    /**
     * @brief Parses the BEGIN_EDGES ... END_EDGES block (phase 2).
     * Produces an empty edge list if the block is present but empty.
     */
    static void parseEdges(std::istream& in, VisSession& session,
                           const HeaderCounts& counts);

    /**
     * @brief Parses the BEGIN_CONTIG_STEPS ... END_CONTIG_STEPS block.
     * Expands RUN blocks into individual BaseAppended steps.
     */
    static void parseContigSteps(std::istream& in, VisSession& session,
                                 const HeaderCounts& counts);

    /**
     * @brief Parses the BEGIN_EULER_STEPS ... END_EULER_STEPS block.
     * Produces an empty step list if the block is present but empty.
     */
    static void parseEulerSteps(std::istream& in, VisSession& session,
                                const HeaderCounts& counts);
};

#endif // DATA_LOADER_H