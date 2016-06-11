#ifndef READCOUNTER_H
#define READCOUNTER_H

#include <unordered_map>
#include <unordered_set>
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
typedef std::unordered_map<std::string, std::vector<std::string>> UniqueBarcodes;

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

class InsertNode;
class Node;

class ReadCounter
{
public:
    ReadCounter(std::vector<Experiment*>, uint16_t insert_mismatches = 1);
    ~ReadCounter();

    void countReads(const std::string&, const std::string&, int threads=1);

    uint16_t allowedMismatches() const;

    uint64_t read() const;
    uint64_t counted() const;
    uint64_t unmatchedTotal() const;
    uint64_t unmatchedInsert() const;
    uint64_t unmatchedBarcodes() const;
    uint64_t unmatchedInsertSequence() const;
    uint64_t written() const;

    UniqueBarcodes uniqueForwardBarcodes() const;
    UniqueBarcodes uniqueReverseBarcodes() const;

private:
    struct ThreadSynchronization;

    uint16_t m_allowedMismatches;

    uint64_t m_read;
    uint64_t m_counted;
    uint64_t m_unmatchedTotal;
    uint64_t m_unmatchedInsert;
    uint64_t m_unmatchedBarcodes;
    uint64_t m_unmatchedInsertSequence;
    uint64_t m_written;

    std::vector<Experiment*> m_experiments;

    std::vector<InsertNode*> m_tree;
    std::vector<Node*> m_nodes;
    UniqueBarcodes m_uniqueFwCodes;
    UniqueBarcodes m_uniqueRevCodes;

    void readFile(const std::string&, ThreadSynchronization*);
    void matchRead(ThreadSynchronization*);
    void writeReads(const std::string&, ThreadSynchronization*);
};

#endif
