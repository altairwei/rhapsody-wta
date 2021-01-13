#ifndef SRC_UTILS_H
#define SRC_UTILS_H

#include <string>
#include <iostream>

int getRandomRange(const int& lowerBound, const int& upperBound);

class ProgressBar
{

public:
    ProgressBar(
        std::size_t total,
        std::string &prefix,
        std::string &suffix,
        int ncol = 60,
        std::ostream &file = std::cerr);
    ~ProgressBar();

    void update(const std::size_t amount);

private:
    std::size_t count;
    std::size_t done;
    std::string &prefix;
    std::string &suffix;
    int ncol;
    std::ostream &file;
};

#endif // SRC_UTILS_H