#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "util.h"

#include<string>
#include <unordered_map>
#include <unordered_set>

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


#endif
