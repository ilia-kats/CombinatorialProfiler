#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

std::string revCompl(const std::string &);

typedef std::unordered_map<std::string, std::vector<std::string>> UniqueBarcodes;
UniqueBarcodes makeUnique(const std::unordered_set<std::string>&, uint16_t);

namespace std
{
    template<> struct hash<std::pair<std::string, std::string>>
    {
        size_t operator()(const std::pair<std::string, std::string> &p) const;
    };
}

#endif
