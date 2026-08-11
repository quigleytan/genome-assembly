#ifndef NODE_NOT_FOUND_EXCEPTION_H
#define NODE_NOT_FOUND_EXCEPTION_H

#include <algorithm>
#include <exception>
#include <string>
#include <cstdint>

class NodeNotFoundException : public std::exception {

private:
    std::string message;

    static std::string uint128ToString(__uint128_t v) {
        if (v == 0) return "0";
        std::string s;
        while (v > 0) { s += ('0' + static_cast<int>(v % 10)); v /= 10; }
        std::reverse(s.begin(), s.end());
        return s;
    }

public:

    explicit NodeNotFoundException(__uint128_t node)
        : message("Node not found: " + uint128ToString(node)) {}

    [[nodiscard]] const char* what() const noexcept override {
        return message.c_str();
    }
};

#endif //NODE_NOT_FOUND_EXCEPTION_H