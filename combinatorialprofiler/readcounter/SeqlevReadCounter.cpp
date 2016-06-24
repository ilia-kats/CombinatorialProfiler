#include "SeqlevReadCounter.h"

SeqlevReadCounter* SeqlevReadCounter::getReadCounter(std::vector<Experiment*> experiments, uint16_t insert_mismatches, uint16_t allowed_barcode_mismatches)
{
    SeqlevReadCounter *counter = new SeqlevReadCounter(experiments, insert_mismatches, allowed_barcode_mismatches);
    counter->init();
    return counter;
}

SeqlevReadCounter::SeqlevReadCounter(std::vector<Experiment*> experiments, uint16_t insert_mismatches, uint16_t allowed_barcode_mismatches)
: ReadCounter(std::move(experiments), insert_mismatches),  m_allowedBarcodeMismatches(allowed_barcode_mismatches)
{}

uint16_t SeqlevReadCounter::allowedBarcodeMismatches() const
{
    return m_allowedBarcodeMismatches;
}

SeqlevBarcodeNode* SeqlevReadCounter::makeFwCodeNode(const Experiment *e, const std::string &seq)
{
    return new FwSeqlevBarcodeNode(seq, m_allowedBarcodeMismatches);
}

SeqlevBarcodeNode* SeqlevReadCounter::makeRevCodeNode(const Experiment *e, const std::string &seq)
{
    return new RevSeqlevBarcodeNode(seq, m_allowedBarcodeMismatches);
}

SeqlevBarcodeNode* SeqlevReadCounter::makeDummyCodeNode()
{
    return new DummySeqlevBarcodeNode();
}
