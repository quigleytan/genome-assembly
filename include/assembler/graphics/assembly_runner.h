/*
 * assembly_runner.h
 * Summary:
 * - Defines AssemblyConfig (user parameters), AssemblyProgress (shared state
 *   between worker and render threads), and the runAssembly() entry point.
 * - runAssembly() is called on a std::thread by gui.cpp. It updates
 *   AssemblyProgress atomically so the render thread can display live progress
 *   without blocking.
 * Important notes:
 * - AssemblyProgress::stage is atomic - safe to read from the render thread
 *   at any time without a mutex.
 * - AssemblyProgress::message and errorMessage are guarded by messageMutex
 *   because std::string is not trivially copyable.
 * - outputPath is set once when stage reaches Complete. The render thread
 *   should call tryLoad(progress.outputPath) at that point.
 */

#ifndef ASSEMBLY_RUNNER_H
#define ASSEMBLY_RUNNER_H

#include <string>
#include <atomic>
#include <mutex>

// ASSEMBLY CONFIG
// Parameters set by the user in the UI.
// Passed by value to runAssembly().

struct AssemblyConfig {
    std::string inputPath;                       // Path to the FASTQ file
    std::string outputDir;                       // Directory to write output files
    size_t      k                    = 9;        // K-mer size
    std::string strategy             = "scored"; // "skip", "greedy", "scored"
    size_t      estimatedTotalBases  = 100000;
};

// ASSEMBLY PROGRESS
// Written by the worker thread, read by the
// render thread every frame.

struct AssemblyProgress {

    enum class Stage {
        Idle,
        LoadingReads,
        BuildingGraph,
        AssemblingContigs,
        Scaffolding,
        ResolvingGaps,
        WritingOutput,
        Complete,
        Failed
    };

    std::atomic<Stage> stage    { Stage::Idle };
    std::atomic<float> progress { 0.0f };  // 0.0 to 1.0 within current stage

    // Strings guarded by messageMutex - not safe to read without locking
    std::mutex  messageMutex;
    std::string message;       // Human-readable status text
    std::string errorMessage;  // Set when stage == Failed
    std::string outputPath;    // Set when stage == Complete

    void setMessage(const std::string& msg) {
        std::lock_guard<std::mutex> lock(messageMutex);
        message = msg;
    }

    std::string getMessage() {
        std::lock_guard<std::mutex> lock(messageMutex);
        return message;
    }

    void setError(const std::string& err) {
        std::lock_guard<std::mutex> lock(messageMutex);
        errorMessage = err;
    }

    std::string getError() {
        std::lock_guard<std::mutex> lock(messageMutex);
        return errorMessage;
    }

    void setOutputPath(const std::string& path) {
        std::lock_guard<std::mutex> lock(messageMutex);
        outputPath = path;
    }

    std::string getOutputPath() {
        std::lock_guard<std::mutex> lock(messageMutex);
        return outputPath;
    }
};

// ENTRY POINT

/**
 * @brief Runs the full assembly pipeline on the calling thread.
 *
 * On success: progress.stage = Complete, progress.outputPath set.
 * On failure: progress.stage = Failed,   progress.errorMessage set.
 *
 * @param config   Assembly parameters chosen by the user.
 * @param progress Shared progress struct updated throughout the run.
 */
void runAssembly(AssemblyConfig config, AssemblyProgress& progress);

#endif // ASSEMBLY_RUNNER_H