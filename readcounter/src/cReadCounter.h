#ifndef READCOUNTER_H
#define READCOUNTER_H

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <functional>
#include <string>
#include <vector>

namespace std
{
template<> struct hash<std::pair<std::string, std::string>>
{
    size_t operator()(const std::pair<std::string, std::string> &p) const;
};
}

typedef std::unordered_map<std::string, std::string> SequenceSet;
typedef std::unordered_map<std::string, std::unordered_map<std::string, double>> SortedCellCounts;
typedef std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>> Counts;

enum NDSIS {noNDSI = 0, forward = 1, reverse = 2};
class Experiment
{
public:
    Experiment();
    Experiment(std::string);

    std::string name;
    std::string insert;

    SequenceSet fwBarcodeSet;
    SequenceSet revBarcodeSet;
    SequenceSet namedInserts;
    NDSIS ndsi;
    SortedCellCounts sortedCells;

    Counts counts;
};

class ReadCounter
{
public:
    ReadCounter(std::vector<std::shared_ptr<Experiment>>);
    void countReads(const std::string&, const std::string&, int threads=1);
    const std::unordered_map<std::string, std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>>>& getCounts() const;

    uint64_t read() const;
    uint64_t counted() const;
    uint64_t written() const;

private:
    struct ThreadSynchronization;
    std::unordered_map<std::string, std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>>> m_counts;

    uint64_t m_read;
    uint64_t m_counted;
    uint64_t m_written;

    std::vector<std::shared_ptr<Experiment>> m_experiments;


    void readFile(const std::string&, ThreadSynchronization*);
    void matchRead(ThreadSynchronization*);
    void writeReads(const std::string&, ThreadSynchronization*);
};

#endif
