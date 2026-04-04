/*
 * DNASequenceException.h
 * Created by Tanner Quigley on 1/29/2026
 * Summary:
 *  - Custom exception used to flag composition problems of input sequences.
 */

#ifndef DNA_SEQUENCE_ERROR_H
#define DNA_SEQUENCE_ERROR_H

#include <exception>
#include <string>
#include <utility>

class DNASequenceException : public std::exception {

private:

    std::string message;

public:
    // Constructor
    explicit DNASequenceException(std::string  msg) : message(std::move(msg)) {}

    // Overwriting of what() base method
    [[nodiscard]] const char* what() const noexcept override {
        return message.c_str();
    }
};

#endif //DNA_SEQUENCE_ERROR_H