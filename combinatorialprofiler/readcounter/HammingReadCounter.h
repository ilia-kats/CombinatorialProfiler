#ifndef HAMMING_READCOUNTER_H
#define HAMMING_READCOUNTER_H

#include "ReadCounter.h"
#include "util.h"

class HammingReadCounter : public ReadCounter
{
public:
    static HammingReadCounter* getReadCounter(std::vector<Experiment*>, uint16_t insert_mismatches = 1, uint16_t unique_barcode_length = 0, float allowed_barcode_mismatches = 0);

    uint16_t minimumUniqueBarcodeLength() const;
    float allowedBarcodeMismatches() const;

    std::unordered_map<std::string, UniqueBarcodes> uniqueForwardBarcodes() const;
    std::unordered_map<std::string, UniqueBarcodes> uniqueReverseBarcodes() const;

private:
    HammingReadCounter(std::vector<Experiment*>, uint16_t insert_mismatches, uint16_t unique_barcode_length, float allowed_barcode_mismatches);

    uint16_t m_uniqueBarcodeLength;
    float m_allowedBarcodeMismatches;

    std::unordered_map<std::string, UniqueBarcodes> m_uniqueFwCodes;
    std::unordered_map<std::string, UniqueBarcodes> m_uniqueRevCodes;

    virtual NodeBase* makeFwCodeNode(const Experiment*, const std::string&);
    virtual NodeBase* makeRevCodeNode(const Experiment*, const std::string&);
    virtual NodeBase* makeDummyCodeNode();
};

#endif
