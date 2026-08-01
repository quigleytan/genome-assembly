/*
 * gui.cpp
 * Summary:
 * - Main visualizer executable. Handles window creation, ImGui setup,
 *   file loading via path input or native Windows file picker,
 *   assembly configuration panel, background assembly thread,
 *   and the main render loop.
 */

#include <iostream>
#include <string>
#include <memory>
#include <chrono>
#include <thread>
#include <filesystem>
#include <functional>

// Windows native file dialog
#include <windows.h>

// GLFW + OpenGL
#include "GLFW/glfw3.h"

// Dear ImGui
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// Project
#include "assembler/graphics/contig_view.h"
#include "assembler/graphics/vis_loader.h"
#include "assembler/graphics/assembly_runner.h"
#include "assembler/core/vis_data.h"

// CONSTANTS

static constexpr int  WINDOW_WIDTH  = 1280;
static constexpr int  WINDOW_HEIGHT = 800;
static constexpr char WINDOW_TITLE[] = "Genome Assembler Visualizer";
static constexpr char GLSL_VERSION[] = "#version 330";


// GLFW ERROR


static void glfwErrorCallback(int error, const char* description) {
    std::cerr << "GLFW error " << error << ": " << description << "\n";
}

// IMGUI STYLE

static void applyStyle() {
    ImGuiStyle& style = ImGui::GetStyle();
    ImGui::StyleColorsDark();

    style.WindowRounding    = 4.0f;
    style.FrameRounding     = 3.0f;
    style.ScrollbarRounding = 3.0f;
    style.GrabRounding      = 3.0f;
    style.ItemSpacing       = ImVec2(8.0f, 5.0f);
    style.FramePadding      = ImVec2(6.0f, 4.0f);

    ImVec4* colors = style.Colors;
    colors[ImGuiCol_WindowBg]      = ImVec4(0.12f, 0.12f, 0.14f, 1.00f);
    colors[ImGuiCol_ChildBg]       = ImVec4(0.10f, 0.10f, 0.12f, 1.00f);
    colors[ImGuiCol_FrameBg]       = ImVec4(0.18f, 0.18f, 0.22f, 1.00f);
    colors[ImGuiCol_SliderGrab]    = ImVec4(0.30f, 0.60f, 0.90f, 1.00f);
    colors[ImGuiCol_Button]        = ImVec4(0.22f, 0.22f, 0.28f, 1.00f);
    colors[ImGuiCol_ButtonHovered] = ImVec4(0.30f, 0.30f, 0.40f, 1.00f);
}

// FILE DIALOGS
// Opens a native Windows file picker for .visdata files

static std::string openVisdataDialog() {
    char pathBuf[MAX_PATH] = "";
    OPENFILENAMEA ofn = {};
    ofn.lStructSize = sizeof(ofn);
    ofn.lpstrFilter = "Visdata Files\0*.visdata\0All Files\0*.*\0";
    ofn.lpstrFile   = pathBuf;
    ofn.nMaxFile    = MAX_PATH;
    ofn.lpstrTitle  = "Open .visdata file";
    ofn.Flags       = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (GetOpenFileNameA(&ofn))
        return std::string(pathBuf);
    return "";
}

// Opens a native Windows file picker for FASTQ files
static std::string openFastqDialog() {
    char pathBuf[MAX_PATH] = "";
    OPENFILENAMEA ofn = {};
    ofn.lStructSize = sizeof(ofn);
    ofn.lpstrFilter = "FASTQ Files\0*.fastq;*.fq\0All Files\0*.*\0";
    ofn.lpstrFile   = pathBuf;
    ofn.nMaxFile    = MAX_PATH;
    ofn.lpstrTitle  = "Select FASTQ file";
    ofn.Flags       = OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (GetOpenFileNameA(&ofn))
        return std::string(pathBuf);
    return "";
}

// Opens a file picker then copies the selected file into destDir.
// Returns the destination path, or "" if cancelled or copy failed.
static std::string uploadFastqToDir(const std::string& destDir) {
    std::string src = openFastqDialog();
    if (src.empty()) return "";

    try {
        std::filesystem::path srcPath(src);
        std::filesystem::path dstPath =
            std::filesystem::path(destDir) / srcPath.filename();
        std::filesystem::copy_file(
            srcPath, dstPath,
            std::filesystem::copy_options::overwrite_existing);
        return dstPath.string();
    } catch (const std::exception& e) {
        std::cerr << "File copy failed: " << e.what() << "\n";
        return "";
    }
}

// ERROR OVERLAY

static void renderErrorOverlay(const std::string& message) {
    ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(
        ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y - 60.0f),
        ImGuiCond_Always, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowBgAlpha(0.80f);
    ImGui::Begin("##error", nullptr,
                 ImGuiWindowFlags_NoDecoration  |
                 ImGuiWindowFlags_NoInputs       |
                 ImGuiWindowFlags_NoNav          |
                 ImGuiWindowFlags_NoMove         |
                 ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f),
                       "Error: %s", message.c_str());
    ImGui::End();
}

// MENU BAR

static void renderMenuBar(const std::string& currentFile,
                          bool& showLaunchPanel) {
    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("New Assembly / Load File..."))
                showLaunchPanel = true;
            ImGui::Separator();
            if (ImGui::MenuItem("Quit"))
                glfwSetWindowShouldClose(glfwGetCurrentContext(), GLFW_TRUE);
            ImGui::EndMenu();
        }

        if (!currentFile.empty()) {
            float textWidth = ImGui::CalcTextSize(currentFile.c_str()).x;
            ImGui::SetCursorPosX(
                ImGui::GetContentRegionMax().x - textWidth - 8.0f);
            ImGui::TextDisabled("%s", currentFile.c_str());
        }

        ImGui::EndMainMenuBar();
    }
}

// PROGRESS BAR LABEL
// Returns a human-readable label for each stage.

static const char* stageLabel(AssemblyProgress::Stage stage) {
    switch (stage) {
        case AssemblyProgress::Stage::Idle:              return "Idle";
        case AssemblyProgress::Stage::LoadingReads:      return "Loading reads";
        case AssemblyProgress::Stage::BuildingGraph:     return "Building graph";
        case AssemblyProgress::Stage::AssemblingContigs: return "Assembling contigs";
        case AssemblyProgress::Stage::Scaffolding:       return "Scaffolding";
        case AssemblyProgress::Stage::WritingOutput:     return "Writing output";
        case AssemblyProgress::Stage::Complete:          return "Complete";
        case AssemblyProgress::Stage::Failed:            return "Failed";
        default:                                         return "";
    }
}

// Maps stage to an overall pipeline fraction (0.0 – 1.0)
// so the progress bar moves smoothly across all stages.
static float stageBaseFraction(AssemblyProgress::Stage stage) {
    switch (stage) {
        case AssemblyProgress::Stage::Idle:              return 0.00f;
        case AssemblyProgress::Stage::LoadingReads:      return 0.00f;
        case AssemblyProgress::Stage::BuildingGraph:     return 0.20f;
        case AssemblyProgress::Stage::AssemblingContigs: return 0.40f;
        case AssemblyProgress::Stage::Scaffolding:       return 0.65f;
        case AssemblyProgress::Stage::WritingOutput:     return 0.85f;
        case AssemblyProgress::Stage::Complete:          return 1.00f;
        case AssemblyProgress::Stage::Failed:            return 1.00f;
        default:                                         return 0.00f;
    }
}

static float stageWidth(AssemblyProgress::Stage stage) {
    switch (stage) {
        case AssemblyProgress::Stage::LoadingReads:      return 0.20f;
        case AssemblyProgress::Stage::BuildingGraph:     return 0.20f;
        case AssemblyProgress::Stage::AssemblingContigs: return 0.25f;
        case AssemblyProgress::Stage::Scaffolding:       return 0.20f;
        case AssemblyProgress::Stage::WritingOutput:     return 0.15f;
        default:                                         return 0.00f;
    }
}

// LAUNCH PANEL
// Shown when no file is loaded or the user
// clicks File > New Assembly / Load File.
// Allows running a new assembly or loading
// an existing .visdata file.


static void renderLaunchPanel(
    bool&             showPanel,
    AssemblyConfig&   config,
    AssemblyProgress& progress,
    std::thread&      assemblyThread,
    bool&             threadRunning,
    const std::function<void(const std::string&)>& tryLoad)
{
    ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(
        ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f),
        ImGuiCond_Always, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(580, 420), ImGuiCond_Always);

    ImGui::Begin("Genome Assembler",
                 nullptr,
                 ImGuiWindowFlags_NoResize   |
                 ImGuiWindowFlags_NoMove     |
                 ImGuiWindowFlags_NoCollapse);

    ImGui::TextColored(ImVec4(0.9f, 0.9f, 0.5f, 1.0f), "Run Assembly");
    ImGui::Separator();
    ImGui::Spacing();

    // FASTQ input
    static char fastqBuf[512] = "";
    ImGui::Text("FASTQ Input File:");
    ImGui::SetNextItemWidth(-220.0f);
    ImGui::InputText("##fastqpath", fastqBuf, sizeof(fastqBuf));
    ImGui::SameLine();

    // Browse - pick an existing file
    if (ImGui::Button("Browse...", ImVec2(100, 0))) {
        std::string picked = openFastqDialog();
        if (!picked.empty()) {
            strncpy(fastqBuf, picked.c_str(), sizeof(fastqBuf) - 1);
            fastqBuf[sizeof(fastqBuf) - 1] = '\0';
        }
    }
    ImGui::SameLine();

    // Upload - copy a file into the Data directory
    if (ImGui::Button("Upload", ImVec2(80, 0))) {
        // Upload into the same directory as the current working directory
        std::string dest = uploadFastqToDir("../Data");
        if (!dest.empty()) {
            strncpy(fastqBuf, dest.c_str(), sizeof(fastqBuf) - 1);
            fastqBuf[sizeof(fastqBuf) - 1] = '\0';
        }
    }

    ImGui::Spacing();

    // Output directory
    static char outputBuf[512] = "../Data/Results";
    ImGui::Text("Output Directory:");
    ImGui::SetNextItemWidth(-8.0f);
    ImGui::InputText("##outputdir", outputBuf, sizeof(outputBuf));

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    // K value
    static int kVal = 9;
    ImGui::Text("K-mer size:");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(80.0f);
    ImGui::SliderInt("##kval", &kVal, 2, 63);
    ImGui::SameLine();
    ImGui::TextDisabled("(2 – 63)");

    ImGui::Spacing();

    // Strategy radio buttons
    static int strategyIdx = 2; // 0=skip, 1=greedy, 2=scored
    ImGui::Text("Resolution strategy:");
    ImGui::SameLine();
    ImGui::RadioButton("Skip",   &strategyIdx, 0); ImGui::SameLine();
    ImGui::RadioButton("Greedy", &strategyIdx, 1); ImGui::SameLine();
    ImGui::RadioButton("Scored", &strategyIdx, 2);

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    // Run button
    bool isRunning = threadRunning &&
                     progress.stage != AssemblyProgress::Stage::Complete &&
                     progress.stage != AssemblyProgress::Stage::Failed;

    if (isRunning) ImGui::BeginDisabled();

    if (ImGui::Button("Run Assembly", ImVec2(140, 32))) {
        if (fastqBuf[0] != '\0') {
            // Build config from UI state
            const char* strategies[] = { "skip", "greedy", "scored" };
            config.inputPath           = std::string(fastqBuf);
            config.outputDir           = std::string(outputBuf);
            config.k                   = static_cast<size_t>(kVal);
            config.strategy            = strategies[strategyIdx];
            config.estimatedTotalBases = 100000;

            // Reset progress
            progress.stage    = AssemblyProgress::Stage::Idle;
            progress.progress = 0.0f;
            progress.setMessage("");
            progress.setError("");
            progress.setOutputPath("");

            // Launch worker thread
            if (assemblyThread.joinable())
                assemblyThread.join();

            threadRunning  = true;
            assemblyThread = std::thread(runAssembly, config,
                                         std::ref(progress));
        }
    }

    if (isRunning) ImGui::EndDisabled();

    // Progress display bar
    auto stage = progress.stage.load();

    if (stage != AssemblyProgress::Stage::Idle) {
        ImGui::Spacing();

        // Overall progress bar combining stage fraction + within-stage progress
        float overall = stageBaseFraction(stage)
                      + stageWidth(stage) * progress.progress.load();
        overall = std::min(overall, 1.0f);

        ImGui::SetNextItemWidth(-8.0f);
        ImGui::ProgressBar(overall, ImVec2(-1, 20),
                           stageLabel(stage));

        // Status message
        std::string msg = progress.getMessage();
        if (!msg.empty()) {
            ImGui::TextDisabled("%s", msg.c_str());
        }

        // Error display
        if (stage == AssemblyProgress::Stage::Failed) {
            ImGui::TextColored(ImVec4(1.0f, 0.3f, 0.3f, 1.0f),
                               "Error: %s", progress.getError().c_str());
        }

        // Auto-load when complete
        if (stage == AssemblyProgress::Stage::Complete) {
            ImGui::TextColored(ImVec4(0.3f, 1.0f, 0.3f, 1.0f),
                               "Assembly complete!");
            std::string outPath = progress.getOutputPath();
            if (!outPath.empty()) {
                ImGui::TextDisabled("Output: %s", outPath.c_str());
                ImGui::Spacing();
                if (ImGui::Button("Open in Visualizer", ImVec2(160, 28))) {
                    tryLoad(outPath);
                    showPanel = false;
                }
            }
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::TextDisabled("─────────────────────────────────────────────");
    ImGui::Spacing();

    // Section 2: Load existing .visdata
    ImGui::TextColored(ImVec4(0.9f, 0.9f, 0.5f, 1.0f),
                       "Load Existing .visdata");
    ImGui::Separator();
    ImGui::Spacing();

    static char visdataBuf[512] = "";
    ImGui::SetNextItemWidth(-200.0f);
    bool pressedEnter = ImGui::InputText(
        "##visdatapath", visdataBuf, sizeof(visdataBuf),
        ImGuiInputTextFlags_EnterReturnsTrue);
    ImGui::SameLine();

    if (ImGui::Button("Browse##vd", ImVec2(90, 0))) {
        std::string picked = openVisdataDialog();
        if (!picked.empty()) {
            strncpy(visdataBuf, picked.c_str(), sizeof(visdataBuf) - 1);
            visdataBuf[sizeof(visdataBuf) - 1] = '\0';
        }
    }
    ImGui::SameLine();

    bool loadPressed = pressedEnter ||
                       ImGui::Button("Load##vd", ImVec2(80, 0));
    if (loadPressed && visdataBuf[0] != '\0') {
        tryLoad(std::string(visdataBuf));
        if (progress.stage != AssemblyProgress::Stage::Failed)
            showPanel = false;
        visdataBuf[0] = '\0';
    }

    ImGui::End();
}

// MAIN

int main(int argc, char** argv) {

    std::string initialPath;
    if (argc >= 2)
        initialPath = argv[1];

    glfwSetErrorCallback(glfwErrorCallback);
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return 1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(
        WINDOW_WIDTH, WINDOW_HEIGHT, WINDOW_TITLE, nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return 1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    applyStyle();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(GLSL_VERSION);

    // Application state
    std::unique_ptr<VisSession>   session;
    std::unique_ptr<ContigView>   contigView;
    std::string                   loadedFilePath;
    std::string                   errorMessage;
    float                         errorTimer    = 0.0f;
    bool                          showPanel     = initialPath.empty();

    AssemblyConfig   assemblyConfig;
    AssemblyProgress assemblyProgress;
    std::thread      assemblyThread;
    bool             threadRunning = false;

    // tryLoad - loads a .visdata file into session + contigView
    auto tryLoad = [&](const std::string& path) {
        if (path.empty()) return;
        try {
            auto s = std::make_unique<VisSession>(VisLoader::load(path));
            auto v = std::make_unique<ContigView>(*s);
            session        = std::move(s);
            contigView     = std::move(v);
            loadedFilePath = path;
            errorMessage.clear();
            showPanel = false;
        } catch (const std::exception& e) {
            errorMessage = e.what();
            errorTimer   = 4.0f;
        }
    };

    if (!initialPath.empty())
        tryLoad(initialPath);

    // Frame timing
    using Clock     = std::chrono::steady_clock;
    using TimePoint = Clock::time_point;
    TimePoint lastFrame = Clock::now();

    // Frame loop
    while (!glfwWindowShouldClose(window)) {

        glfwPollEvents();

        TimePoint now       = Clock::now();
        float     deltaTime = std::chrono::duration<float>(now - lastFrame).count();
        lastFrame  = now;
        deltaTime  = std::min(deltaTime, 0.1f);

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Menu bar
        renderMenuBar(loadedFilePath, showPanel);

        // Launch panel
        if (showPanel) {
            renderLaunchPanel(showPanel, assemblyConfig, assemblyProgress,
                              assemblyThread, threadRunning, tryLoad);
        }

        // Error overlay
        if (errorTimer > 0.0f) {
            renderErrorOverlay(errorMessage);
            errorTimer -= deltaTime;
        }

        // Main view
        if (contigView) {
            contigView->update(deltaTime);
            contigView->render();
        } else if (!showPanel) {
            ImGui::SetNextWindowPos(
                ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f),
                ImGuiCond_Always, ImVec2(0.5f, 0.5f));
            ImGui::Begin("##hint", nullptr,
                         ImGuiWindowFlags_NoDecoration     |
                         ImGuiWindowFlags_AlwaysAutoResize |
                         ImGuiWindowFlags_NoMove);
            ImGui::TextDisabled("No file loaded.");
            ImGui::TextDisabled("Use File > New Assembly / Load File...");
            ImGui::End();
        }

        ImGui::Render();

        int displayW, displayH;
        glfwGetFramebufferSize(window, &displayW, &displayH);
        glViewport(0, 0, displayW, displayH);
        glClearColor(0.10f, 0.10f, 0.12f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // Cleanup
    // Join the worker thread before destroying anything it might reference
    if (assemblyThread.joinable())
        assemblyThread.join();

    contigView.reset();
    session.reset();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}