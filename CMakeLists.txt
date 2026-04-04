cmake_minimum_required(VERSION 3.14)
project(M4OEP-tequigle)

set(CMAKE_CXX_STANDARD 17)

# DEPENDENCY VERSIONS

set(GLFW_VERSION  3.3.9)
set(GLM_VERSION   1.0.1)
set(IMGUI_VERSION 1.90.4)

# FETCH LIBRARIES

include(FetchContent)

FetchContent_Declare(
        glfw
        URL https://github.com/glfw/glfw/archive/refs/tags/${GLFW_VERSION}.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

FetchContent_Declare(
        glm
        URL https://github.com/g-truc/glm/archive/refs/tags/${GLM_VERSION}.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(glm)

FetchContent_Declare(
        imgui
        URL https://github.com/ocornut/imgui/archive/refs/tags/v${IMGUI_VERSION}.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(imgui)

# IMGUI STATIC LIBRARY

add_library(ImGui STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_demo.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(ImGui PUBLIC
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
        ${glfw_SOURCE_DIR}/include
)

target_link_libraries(ImGui PUBLIC glfw)
if (WIN32)
    target_link_libraries(ImGui PUBLIC opengl32)
endif()

# SHARED INCLUDE ROOT

set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR})

# ASSEMBLY EXECUTABLES

add_executable(FullReconstruction
        src/construction/gap_estimation.cpp
        src/construction/gap_filling.cpp
)

set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR})

add_executable(ScaffoldAssembly
        src/pipelines/scaffold_assembly.cpp
        src/core/data_list.cpp
        src/core/dna_sequence.cpp
        src/io/sequence_reader.cpp
        src/construction/kmer_encoding.cpp
        src/core/de_bruijn_graph.cpp
        src/construction/contig_assembler.cpp
        src/construction/scaffolder.cpp
        src/graphics/data_exporter.cpp
)
target_include_directories(ScaffoldAssembly PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(FastaEulerianAssembly
        src/pipelines/fasta_eulerian_assembly.cpp
        src/core/data_list.cpp
        src/core/dna_sequence.cpp
        src/io/sequence_reader.cpp
        src/construction/kmer_encoding.cpp
        src/core/de_bruijn_graph.cpp
        src/construction/eulerian_assembler.cpp
)
target_include_directories(FastaEulerianAssembly PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(FastaContigAssembly
        src/pipelines/fasta_contig_assembly.cpp
        src/core/data_list.cpp
        src/core/dna_sequence.cpp
        src/io/sequence_reader.cpp
        src/construction/kmer_encoding.cpp
        src/core/de_bruijn_graph.cpp
        src/construction/contig_assembler.cpp
        src/construction/scaffolder.cpp
)
target_include_directories(FastaContigAssembly PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(FastqContigAssembly
        src/pipelines/fastq_contig_assembly.cpp
        src/core/data_list.cpp
        src/core/dna_sequence.cpp
        src/io/sequence_reader.cpp
        src/construction/kmer_encoding.cpp
        src/core/de_bruijn_graph.cpp
        src/construction/contig_assembler.cpp
)
target_include_directories(FastqContigAssembly PRIVATE ${PROJECT_INCLUDE_DIR})

# TEST EXECUTABLES

add_executable(InitializationTests
        tests/InitializationTests.cpp
        src/core/dna_sequence.cpp
        src/io/sequence_reader.cpp
        src/construction/kmer_encoding.cpp
)
target_include_directories(InitializationTests PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(ProcessingTests
        tests/ProcessingTests.cpp
        src/construction/kmer_encoding.cpp
        src/core/dna_sequence.cpp
)
target_include_directories(ProcessingTests PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(ConstructionTests
        src/core/data_list.cpp
        tests/ConstructionTests.cpp
        src/core/de_bruijn_graph.cpp
        src/construction/eulerian_assembler.cpp
        src/construction/contig_assembler.cpp
        src/construction/scaffolder.cpp
        src/construction/kmer_encoding.cpp
        src/core/dna_sequence.cpp
)
target_include_directories(ConstructionTests PRIVATE ${PROJECT_INCLUDE_DIR})

# VISUALIZER EXECUTABLE

add_executable(Visualizer
        src/graphics/gui.cpp
        src/graphics/contig_view.cpp
        include/assembler/graphics/GraphView.cpp
        src/graphics/data_exporter.cpp
        src/graphics/data_loader.cpp
        src/core/dna_sequence.cpp
        src/io/sequence_reader.cpp
        src/construction/kmer_encoding.cpp
        src/core/de_bruijn_graph.cpp
        src/construction/eulerian_assembler.cpp
        src/construction/contig_assembler.cpp
        src/construction/scaffolder.cpp
)

target_include_directories(Visualizer PRIVATE
        ${PROJECT_INCLUDE_DIR}
        ${imgui_SOURCE_DIR}
        ${imgui_SOURCE_DIR}/backends
        ${glfw_SOURCE_DIR}/include
        ${glm_SOURCE_DIR}
)
target_link_libraries(Visualizer PRIVATE ImGui glfw)
if (WIN32)
    target_link_libraries(Visualizer PRIVATE opengl32)
endif()