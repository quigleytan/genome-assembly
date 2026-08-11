#include "assembler/io/console_input.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

// PRIVATE HELPERS (translation-unit local)

namespace {

    std::string trim(const std::string& s) {
        size_t start = s.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        size_t end = s.find_last_not_of(" \t\r\n");
        return s.substr(start, end - start + 1);
    }

    // Reads a line from stdin, throwing if the stream is closed/exhausted.
    std::string readLine() {
        std::string line;
        if (!std::getline(std::cin, line))
            throw std::runtime_error("No input received - stdin closed unexpectedly");
        return trim(line);
    }

}

// PUBLIC

std::string ConsoleInput::promptFilePath(const std::string& prompt) {
    while (true) {
        std::cout << prompt;
        std::string path = readLine();

        if (path.empty()) {
            std::cout << "Path cannot be empty. Please try again.\n";
            continue;
        }

        std::ifstream file(path);
        if (file.is_open()) return path;

        std::cout << "Could not open '" << path << "'. Please try again.\n";
    }
}

size_t ConsoleInput::promptSizeT(const std::string& prompt, size_t min, size_t max) {
    while (true) {
        std::cout << prompt << " [" << min << "-" << max << "]: ";
        std::string input = readLine();

        if (input.empty() || input.find_first_not_of("0123456789") != std::string::npos) {
            std::cout << "Please enter a whole number.\n";
            continue;
        }

        try {
            size_t parsed = static_cast<size_t>(std::stoull(input));
            if (parsed < min || parsed > max) {
                std::cout << "Value must be between " << min << " and " << max << ".\n";
                continue;
            }
            return parsed;
        } catch (const std::out_of_range&) {
            std::cout << "Value is too large. Please try again.\n";
        }
    }
}

bool ConsoleInput::promptYesNo(const std::string& prompt, bool defaultValue) {
    while (true) {
        std::cout << prompt << (defaultValue ? " [Y/n]: " : " [y/N]: ");
        std::string input = readLine();

        if (input.empty()) return defaultValue;

        for (char& c : input) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

        if (input == "y" || input == "yes") return true;
        if (input == "n" || input == "no")  return false;

        std::cout << "Please answer 'y' or 'n'.\n";
    }
}

size_t ConsoleInput::promptChoice(const std::string& prompt, const std::vector<std::string>& options) {
    if (options.empty())
        throw std::runtime_error("promptChoice called with no options");

    std::cout << prompt << "\n";
    for (size_t i = 0; i < options.size(); ++i)
        std::cout << "  " << (i + 1) << ") " << options[i] << "\n";

    return promptSizeT("Enter a choice", 1, options.size()) - 1;
}
