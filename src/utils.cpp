#include "utils.h"

#include <cstdlib>

using namespace std;

// Random select from a half open interval: [0, 10)
int getRandomRange(const int& lowerBound, const int& upperBound)
{
    const int range = upperBound - lowerBound;
    return (lowerBound + int(range * double(std::rand()) / double(RAND_MAX)));
}

ProgressBar::ProgressBar(
    std::size_t total,
    std::string &prefix,
    std::string &suffix,
    int ncol /*= 60*/,
    std::ostream &file /*= std::cerr*/)
    : count(total), prefix(prefix), suffix(suffix), ncol(ncol), file(file), done(0)
{
    update(0);
}

void ProgressBar::update(const size_t amount)
{
    done += amount;
    int x = int(ncol*double(done)/double(count));
    file << prefix << " [" << string(x, '#') << string(ncol - x, ' ') << "] "
        << done << '/' << count << " " << suffix << '\r' << flush;
}