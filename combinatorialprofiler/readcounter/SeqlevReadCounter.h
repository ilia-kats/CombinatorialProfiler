#ifndef SEQLEV_READCOUNTER_H
#define SEQLEV_READCOUNTER_H

#include "ReadCounter.h"

class SeqlevReadCounter : public ReadCounter
{
public:
    static SeqlevReadCounter* getReadCounter(std::vector<Experiment*>, uint16_t insert_mismatches = 1, uint16_t allowed_barcode_mismatches = 0);

    uint16_t allowedBarcodeMismatches() const;

private:
    SeqlevReadCounter(std::vector<Experiment*>, uint16_t insert_mismatches, uint16_t allowed_barcode_mismatches);

    uint16_t m_allowedBarcodeMismatches;

    virtual NodeBase* makeFwCodeNode(const Experiment*, const std::string&);
    virtual NodeBase* makeRevCodeNode(const Experiment*, const std::string&);
    virtual NodeBase* makeDummyCodeNode();
};

#endif
