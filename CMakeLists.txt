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
        gap_resolution/gap_estimation.cpp
        gap_resolution/gap_filling.cpp
)

set(PROJECT_INCLUDE_DIR ${CMAKE_SOURCE_DIR})

add_executable(ScaffoldAssembly
        pipelines/ScaffoldAssembly.cpp
        DataList.cpp
        data_initalization/dna_sequence.cpp
        data_initalization/sequence_reader.cpp
        encoding/kmer_encoding.cpp
        ScaffoldConstruction/de_bruijn_graph.cpp
        ScaffoldConstruction/contig_traversal.cpp
        ScaffoldConstruction/scaffolder.cpp
        graphics/data_exporter.cpp
)
target_include_directories(ScaffoldAssembly PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(FastaEulerianAssembly
        pipelines/FastaEulerianAssembly.cpp
        DataList.cpp
        data_initalization/dna_sequence.cpp
        data_initalization/sequence_reader.cpp
        encoding/kmer_encoding.cpp
        ScaffoldConstruction/de_bruijn_graph.cpp
        ScaffoldConstruction/eulerian_traversal.cpp
)
target_include_directories(FastaEulerianAssembly PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(FastaContigAssembly
        pipelines/FastaContigAssembly.cpp
        DataList.cpp
        data_initalization/dna_sequence.cpp
        data_initalization/sequence_reader.cpp
        encoding/kmer_encoding.cpp
        ScaffoldConstruction/de_bruijn_graph.cpp
        ScaffoldConstruction/contig_traversal.cpp
        ScaffoldConstruction/scaffolder.cpp
)
target_include_directories(FastaContigAssembly PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(FastqContigAssembly
        pipelines/FastqContigAssembly.cpp
        DataList.cpp
        data_initalization/dna_sequence.cpp
        data_initalization/sequence_reader.cpp
        encoding/kmer_encoding.cpp
        ScaffoldConstruction/de_bruijn_graph.cpp
        ScaffoldConstruction/contig_traversal.cpp
)
target_include_directories(FastqContigAssembly PRIVATE ${PROJECT_INCLUDE_DIR})

# TEST EXECUTABLES

add_executable(InitializationTests
        data_initalization/InitializationTests.cpp
        data_initalization/dna_sequence.cpp
        data_initalization/sequence_reader.cpp
        encoding/kmer_encoding.cpp
)
target_include_directories(InitializationTests PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(ProcessingTests
        encoding/ProcessingTests.cpp
        encoding/kmer_encoding.cpp
        data_initalization/dna_sequence.cpp
)
target_include_directories(ProcessingTests PRIVATE ${PROJECT_INCLUDE_DIR})

add_executable(ConstructionTests
        DataList.cpp
        ScaffoldConstruction/ConstructionTests.cpp
        ScaffoldConstruction/de_bruijn_graph.cpp
        ScaffoldConstruction/eulerian_traversal.cpp
        ScaffoldConstruction/contig_traversal.cpp
        ScaffoldConstruction/scaffolder.cpp
        encoding/kmer_encoding.cpp
        data_initalization/dna_sequence.cpp
)
target_include_directories(ConstructionTests PRIVATE ${PROJECT_INCLUDE_DIR})

# VISUALIZER EXECUTABLE

add_executable(Visualizer
        graphics/gui.cpp
        graphics/contig_view.cpp
        graphics/GraphView.cpp
        graphics/data_exporter.cpp
        graphics/data_loader.cpp
        data_initalization/dna_sequence.cpp
        data_initalization/sequence_reader.cpp
        encoding/kmer_encoding.cpp
        ScaffoldConstruction/de_bruijn_graph.cpp
        ScaffoldConstruction/eulerian_traversal.cpp
        ScaffoldConstruction/contig_traversal.cpp
        ScaffoldConstruction/scaffolder.cpp
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