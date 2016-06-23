#include "HammingReadCounter.h"

HammingReadCounter* HammingReadCounter::getReadCounter(std::vector<Experiment*> experiments, uint16_t insert_mismatches, uint16_t unique_barcode_length, float allowed_barcode_mismatches)
{
    HammingReadCounter *counter = new HammingReadCounter(experiments, insert_mismatches, unique_barcode_length, allowed_barcode_mismatches);
    counter->init();
    return counter;
}

HammingReadCounter::HammingReadCounter(std::vector<Experiment*> experiments, uint16_t insert_mismatches, uint16_t unique_barcode_length, float allowed_barcode_mismatches)
: ReadCounter(std::move(experiments), insert_mismatches), m_uniqueBarcodeLength(unique_barcode_length), m_allowedBarcodeMismatches(allowed_barcode_mismatches)
{
    std::unordered_map<std::string, std::unordered_set<std::string>> fwCodes;
    std::unordered_map<std::string, std::unordered_set<std::string>> revCodes;

    for (const auto &exp : this->experiments()) {
        for (const auto &c : exp->fwBarcodeSet)
            fwCodes[exp->insert].insert(c.first);
        for (const auto &c : exp->revBarcodeSet)
            revCodes[exp->insert].insert(c.first);
    }

    for (const auto &c : fwCodes) {
        m_uniqueFwCodes[c.first] = makeUnique(c.second, m_uniqueBarcodeLength);
    }
    for (const auto &c : revCodes) {
        m_uniqueRevCodes[c.first] = makeUnique(c.second, m_uniqueBarcodeLength);
    }
}

uint16_t HammingReadCounter::minimumUniqueBarcodeLength() const
{
    return m_uniqueBarcodeLength;
}

float HammingReadCounter::allowedBarcodeMismatches() const
{
    return m_allowedBarcodeMismatches;
}

std::unordered_map<std::string, UniqueBarcodes> HammingReadCounter::uniqueForwardBarcodes() const
{
    return m_uniqueFwCodes;
}

std::unordered_map<std::string, UniqueBarcodes> HammingReadCounter::uniqueReverseBarcodes() const
{
    return m_uniqueRevCodes;
}

HammingBarcodeNode* HammingReadCounter::makeFwCodeNode(const Experiment *e, const std::string &seq)
{
    return new FwHammingBarcodeNode(seq, m_allowedBarcodeMismatches, m_uniqueFwCodes[e->insert][seq]);
}

HammingBarcodeNode* HammingReadCounter::makeRevCodeNode(const Experiment *e, const std::string &seq)
{
    return new RevHammingBarcodeNode(seq, m_allowedBarcodeMismatches, m_uniqueRevCodes[e->insert][seq]);
}

HammingBarcodeNode* HammingReadCounter::makeDummyCodeNode()
{
    return new DummyHammingBarcodeNode();
}
