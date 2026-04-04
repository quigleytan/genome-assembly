/*
 * VisualizerApp.cpp
 * Summary:
 * -
 */

#include <iostream>
#include <string>
#include <memory>
#include <stdexcept>
#include <chrono>

// GLFW + OpenGL
#include "GLFW/glfw3.h"

// Dear ImGui
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// Project
#include "graphics/contig_view.h"
#include "graphics/data_loader.h"
#include "../../include/core/vis_data.h"

// CONSTANTS

static constexpr int   WINDOW_WIDTH  = 1280;
static constexpr int   WINDOW_HEIGHT = 800;
static constexpr char  WINDOW_TITLE[] = "Genome Assembler Visualizer";
static constexpr char  GLSL_VERSION[] = "#version 330";


// GLFW ERROR

static void glfwErrorCallback(int error, const char* description) {
    std::cerr << "GLFW error " << error << ": " << description << "\n";
}

// IMGUI STYLE

static void applyStyle() {
    ImGuiStyle& style = ImGui::GetStyle();
    ImGui::StyleColorsDark();

    style.WindowRounding   = 4.0f;
    style.FrameRounding    = 3.0f;
    style.ScrollbarRounding = 3.0f;
    style.GrabRounding     = 3.0f;
    style.ItemSpacing      = ImVec2(8.0f, 5.0f);
    style.FramePadding     = ImVec2(6.0f, 4.0f);

    // Slightly warmer background so colored bars pop
    ImVec4* colors = style.Colors;
    colors[ImGuiCol_WindowBg]  = ImVec4(0.12f, 0.12f, 0.14f, 1.00f);
    colors[ImGuiCol_ChildBg]   = ImVec4(0.10f, 0.10f, 0.12f, 1.00f);
    colors[ImGuiCol_FrameBg]   = ImVec4(0.18f, 0.18f, 0.22f, 1.00f);
    colors[ImGuiCol_SliderGrab]= ImVec4(0.30f, 0.60f, 0.90f, 1.00f);
    colors[ImGuiCol_Button]    = ImVec4(0.22f, 0.22f, 0.28f, 1.00f);
    colors[ImGuiCol_ButtonHovered] = ImVec4(0.30f, 0.30f, 0.40f, 1.00f);
}

// LOAD PROMPT

static std::string renderLoadPrompt() {
    static char pathBuf[512] = "";
    std::string result;

    ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(
        ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f),
        ImGuiCond_Always, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(480, 140), ImGuiCond_Always);

    ImGui::Begin("Load .visdata file",
                 nullptr,
                 ImGuiWindowFlags_NoResize |
                 ImGuiWindowFlags_NoMove   |
                 ImGuiWindowFlags_NoCollapse);

    ImGui::TextWrapped("Enter the path to a .visdata file produced by the assembly pipeline:");
    ImGui::Spacing();

    // Focus the input field automatically on first open
    ImGui::SetNextItemWidth(-1);
    bool pressedEnter = ImGui::InputText(
        "##filepath", pathBuf, sizeof(pathBuf),
        ImGuiInputTextFlags_EnterReturnsTrue);

    ImGui::Spacing();

    bool load = pressedEnter || ImGui::Button("Load", ImVec2(80, 0));

    if (load && pathBuf[0] != '\0')
        result = std::string(pathBuf);

    ImGui::End();
    return result;
}

// ERROR OVERLAY

static void renderErrorOverlay(const std::string& message) {
    ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y - 60.0f),
                            ImGuiCond_Always, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowBgAlpha(0.80f);
    ImGui::Begin("##error",
                 nullptr,
                 ImGuiWindowFlags_NoDecoration   |
                 ImGuiWindowFlags_NoInputs        |
                 ImGuiWindowFlags_NoNav           |
                 ImGuiWindowFlags_NoMove          |
                 ImGuiWindowFlags_AlwaysAutoResize);
    ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "Error: %s", message.c_str());
    ImGui::End();
}

// MENU BAR

static std::string renderMenuBar(const std::string& currentFile) {
    std::string result;

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Load .visdata..."))
                result = "__SHOW_PROMPT__"; // sentinel to open the prompt
            ImGui::Separator();
            if (ImGui::MenuItem("Quit"))
                glfwSetWindowShouldClose(
                    glfwGetCurrentContext(), GLFW_TRUE);
            ImGui::EndMenu();
        }

        // Show the current filename on the right side of the menu bar
        if (!currentFile.empty()) {
            float textWidth = ImGui::CalcTextSize(currentFile.c_str()).x;
            ImGui::SetCursorPosX(
                ImGui::GetContentRegionMax().x - textWidth - 8.0f);
            ImGui::TextDisabled("%s", currentFile.c_str());
        }

        ImGui::EndMainMenuBar();
    }

    return result;
}

// MAIN

int main(int argc, char** argv) {

    // Initiation steps
    std::string initialPath;
    if (argc >= 2)
        initialPath = argv[1];

    glfwSetErrorCallback(glfwErrorCallback);
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return 1;
    }

    // OpenGL 3.3 core profile:  required by imgui_impl_opengl3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // required on macOS

    GLFWwindow* window = glfwCreateWindow(
        WINDOW_WIDTH, WINDOW_HEIGHT, WINDOW_TITLE, nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return 1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // VSync — caps frame rate to monitor refresh

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;

    applyStyle();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(GLSL_VERSION);

    std::unique_ptr<VisSession>   session;
    std::unique_ptr<ContigView>   contigView;
    std::string                   loadedFilePath;
    std::string                   errorMessage;
    float                         errorTimer    = 0.0f;
    bool                          showPrompt    = initialPath.empty();

    // Load initial file
    auto tryLoad = [&](const std::string& path) {
        try {
            auto s = std::make_unique<VisSession>(DataLoader::load(path));
            auto v = std::make_unique<ContigView>(*s);
            // Only commit once both succeed
            session        = std::move(s);
            contigView     = std::move(v);
            loadedFilePath = path;
            errorMessage.clear();
            showPrompt = false;
        } catch (const std::exception& e) {
            errorMessage = e.what();
            errorTimer   = 4.0f; // show error for 4 seconds
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

        // Delta time
        TimePoint now      = Clock::now();
        float     deltaTime = std::chrono::duration<float>(now - lastFrame).count();
        lastFrame = now;
        deltaTime = std::min(deltaTime, 0.1f); // clamp to avoid spiral on hang

        // ImGui new frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Menu bar
        std::string menuAction = renderMenuBar(loadedFilePath);
        if (menuAction == "__SHOW_PROMPT__")
            showPrompt = true;

        // ── Load prompt (shown when no file loaded or user requests load) ──
        if (showPrompt) {
            std::string promptResult = renderLoadPrompt();
            if (!promptResult.empty())
                tryLoad(promptResult);
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
        } else if (!showPrompt) {
            // No file loaded and prompt is dismissed
            ImGuiIO& imio = ImGui::GetIO();
            ImGui::SetNextWindowPos(
                ImVec2(imio.DisplaySize.x * 0.5f, imio.DisplaySize.y * 0.5f),
                ImGuiCond_Always, ImVec2(0.5f, 0.5f));
            ImGui::Begin("##hint",
                         nullptr,
                         ImGuiWindowFlags_NoDecoration |
                         ImGuiWindowFlags_AlwaysAutoResize |
                         ImGuiWindowFlags_NoMove);
            ImGui::TextDisabled("No file loaded.");
            ImGui::TextDisabled("Use File > Load .visdata... to open a session.");
            ImGui::End();
        }

        // Render step
        ImGui::Render();

        int displayW, displayH;
        glfwGetFramebufferSize(window, &displayW, &displayH);
        glViewport(0, 0, displayW, displayH);
        glClearColor(0.10f, 0.10f, 0.12f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // Deconstructors
    contigView.reset();
    session.reset();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}