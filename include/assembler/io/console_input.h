/*
 * console_input.h
 * Summary:
 * - Shared helpers for gathering pipeline configuration interactively from
 *   stdin, with retry-on-invalid-input loops so a bad entry never crashes
 *   or silently misconfigures a run.
 * Important notes:
 * - Every prompt* function blocks until it receives a valid answer, or
 *   throws std::runtime_error if stdin is closed/exhausted (e.g. piped
 *   input ran out) - callers don't need their own EOF handling.
 */

#ifndef CONSOLE_INPUT_H
#define CONSOLE_INPUT_H

#include <string>
#include <vector>

class ConsoleInput {

public:

    /**
     * @brief Prompts for a file path, retrying until one can be opened for reading.
     * @param prompt Message displayed before reading input.
     * @return Path to a file confirmed to be openable.
     * @throws std::runtime_error if stdin closes before a valid path is entered.
     */
    static std::string promptFilePath(const std::string& prompt);

    /**
     * @brief Prompts for an unsigned integer within [min, max], retrying on
     *        non-numeric input or values outside the range.
     * @param prompt Message displayed before reading input.
     * @param min    Minimum acceptable value (inclusive).
     * @param max    Maximum acceptable value (inclusive).
     * @return A validated value in [min, max].
     * @throws std::runtime_error if stdin closes before a valid value is entered.
     */
    static size_t promptSizeT(const std::string& prompt, size_t min, size_t max);

    /**
     * @brief Prompts for a yes/no answer.
     * @param prompt       Message displayed before reading input (a "[y/n]" hint is appended).
     * @param defaultValue Value returned if the user presses Enter with no input.
     * @return true for yes, false for no.
     * @throws std::runtime_error if stdin closes before a valid value is entered.
     */
    static bool promptYesNo(const std::string& prompt, bool defaultValue = false);

    /**
     * @brief Displays a numbered menu and prompts for a selection.
     * @param prompt  Message displayed above the menu.
     * @param options Menu entries, displayed 1-indexed.
     * @return Index into options (0-indexed) of the selected entry.
     * @throws std::runtime_error if stdin closes before a valid selection is entered,
     *         or if options is empty.
     */
    static size_t promptChoice(const std::string& prompt, const std::vector<std::string>& options);

};

#endif //CONSOLE_INPUT_H
