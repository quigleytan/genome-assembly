cmake_minimum_required(VERSION 3.14)
project(M4OEP-tequigle)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

macro(fetch_dependency name url)
    FetchContent_Declare(
            ${name}
            URL ${url}
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(${name})
endmacro()

set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

fetch_dependency(glfw https://github.com/glfw/glfw/archive/refs/tags/3.3.9.tar.gz)
fetch_dependency(glm https://github.com/g-truc/glm/archive/refs/tags/1.0.1.tar.gz)
fetch_dependency(imgui https://github.com/ocornut/imgui/archive/refs/tags/v1.90.4.tar.gz)

#IMGUI
add_library(ImGui STATIC
        ${imgui_SOURCE_DIR}/imgui.cpp
        ${imgui_SOURCE_DIR}/imgui_draw.cpp
        ${imgui_SOURCE_DIR}/imgui_widgets.cpp
        ${imgui_SOURCE_DIR}/imgui_tables.cpp
        ${imgui_SOURCE_DIR}/imgui_demo.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
        ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(ImGui PUBLIC ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends ${glfw_SOURCE_DIR}/include)
target_link_libraries(ImGui PUBLIC glfw)

if(WIN32)
    target_link_libraries(ImGui PUBLIC opengl32)
endif()

macro(add_assembler_app target_name main_src)
    add_executable(${target_name} ${main_src} ${ARGN})
    target_include_directories(${target_name} PRIVATE "${CMAKE_SOURCE_DIR}/include")
endmacro()

set(CORE_SOURCES
        src/core/data_list.cpp
        src/core/dna_sequence.cpp
        src/core/de_bruijn_graph.cpp
        src/io/sequence_reader.cpp
)

set(CONSTRUCTION_SOURCES
        src/construction/kmer_encoding.cpp
        src/construction/contig_assembler.cpp
        src/construction/eulerian_assembler.cpp
        src/construction/scaffolder.cpp
        src/construction/gap_estimation.cpp
        src/construction/gap_filling.cpp
)

# Executables
add_assembler_app(ScaffoldAssembly      src/pipelines/scaffold_assembly.cpp ${CORE_SOURCES} ${CONSTRUCTION_SOURCES} src/graphics/vis_exporter.cpp)
add_assembler_app(FastaEulerianAssembly src/pipelines/fasta_eulerian_assembly.cpp ${CORE_SOURCES} ${CONSTRUCTION_SOURCES})
add_assembler_app(FastaContigAssembly   src/pipelines/fasta_contig_assembly.cpp ${CORE_SOURCES} ${CONSTRUCTION_SOURCES})
add_assembler_app(FastqContigAssembly   src/pipelines/fastq_contig_assembly.cpp ${CORE_SOURCES} ${CONSTRUCTION_SOURCES})

# Tests
add_assembler_app(InitializationTests   tests/InitializationTests.cpp src/core/dna_sequence.cpp src/io/sequence_reader.cpp src/construction/kmer_encoding.cpp)
add_assembler_app(ProcessingTests       tests/ProcessingTests.cpp src/construction/kmer_encoding.cpp src/core/dna_sequence.cpp)
add_assembler_app(ConstructionTests     tests/ConstructionTests.cpp ${CORE_SOURCES} ${CONSTRUCTION_SOURCES})

# Visualizer
add_assembler_app(Visualizer
        src/graphics/gui.cpp
        src/graphics/contig_view.cpp
        src/graphics/vis_exporter.cpp
        src/graphics/vis_loader.cpp
        ${CORE_SOURCES}
        ${CONSTRUCTION_SOURCES}
)

target_include_directories(Visualizer PRIVATE ${imgui_SOURCE_DIR} ${imgui_SOURCE_DIR}/backends ${glfw_SOURCE_DIR}/include ${glm_SOURCE_DIR})
target_link_libraries(Visualizer PRIVATE ImGui glfw)