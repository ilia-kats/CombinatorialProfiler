#ifndef HAMMING_READCOUNTER_H
#define HAMMING_READCOUNTER_H

#include "ReadCounter.h"

class HammingReadCounter : public ReadCounter<HammingBarcodeNode>
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

    virtual HammingBarcodeNode* makeFwCodeNode(const Experiment*, const std::string&);
    virtual HammingBarcodeNode* makeRevCodeNode(const Experiment*, const std::string&);
    virtual HammingBarcodeNode* makeDummyCodeNode();
};

#endif
