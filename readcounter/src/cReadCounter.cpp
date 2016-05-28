#include "cReadCounter.h"

BarcodeSet::BarcodeSet(const std::unordered_map<std::string, std::list<std::string>> &f, const std::unordered_map<std::string, std::list<std::string>> &r)
: fw(f), rev(r)
{}
