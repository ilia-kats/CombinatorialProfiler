#ifndef READCOUNTER_H
#define READCOUNTER_H
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <numeric>
#include <fstream>
#include <condition_variable>
#include <atomic>
#include <thread>
#include <cmath>

#include <cassert>

class Experiment;
class InsertNode;
class NodeBase;

class ReadCounter
{
public:
    virtual ~ReadCounter();

    void countReads(const std::string&, const std::string&, int threads=1);

    uint16_t allowedMismatches() const;

    uint64_t read() const;

    uint64_t counted() const;

    uint64_t unmatchedTotal() const;

    uint64_t unmatchedInsert() const;

    uint64_t unmatchedBarcodes() const;

    uint64_t unmatchedInsertSequence() const;

    uint64_t written() const;

protected:
    ReadCounter(std::vector<Experiment*> experiments, uint16_t insert_mismatches = 1);

    void init();

    const std::vector<Experiment*>& experiments() const;

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
    std::vector<NodeBase*> m_nodes;

    void readFile(const std::string&, ThreadSynchronization*);

    void matchRead(ThreadSynchronization*);

    void writeReads(const std::string&, ThreadSynchronization*);

    InsertNode* makeInsertNode(const Experiment*) const;
    virtual NodeBase* makeFwCodeNode(const Experiment*, const std::string&) = 0;
    virtual NodeBase* makeRevCodeNode(const Experiment*, const std::string&) = 0;
    virtual NodeBase* makeDummyCodeNode() = 0;
};

#endif
