/*
 * DNASequenceException.h
 * Created by Tanner Quigley on 2/15/2026
 * Summary:
 *  - Custom exception used to flag composition problems of input sequences.
 */
#ifndef NODE_NOT_FOUND_EXCEPTION_H
#define NODE_NOT_FOUND_EXCEPTION_H

#include <exception>
#include <string>
#include <cstdint>

class NodeNotFoundException : public std::exception {

private:

    std::string message;

public:

    explicit NodeNotFoundException(uint64_t node)
        : message("Node not found: " + std::to_string(node)) {}

    // Overwriting of what() base method
    [[nodiscard]] const char* what() const noexcept override {
        return message.c_str();
    }
};
#endif //NODE_NOT_FOUND_EXCEPTION_H