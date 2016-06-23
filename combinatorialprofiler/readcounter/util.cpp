#include "util.h"

#include <algorithm>

std::string revCompl(const std::string &seq)
{
    std::string seq2;
    seq2.reserve(seq.size());
    std::transform(seq.crbegin(), seq.crend(), std::inserter(seq2, seq2.begin()), [](const std::remove_reference<decltype(seq)>::type::value_type &c) -> std::remove_reference<decltype(seq)>::type::value_type {
        switch (c) {
            case 'A':
                return 'T';
            case 'a':
                return 't';
            case 'T':
                return 'A';
            case 't':
                return 'a';
            case 'G':
                return 'C';
            case 'g':
                return 'c';
            case 'C':
                return 'G';
            case 'c':
                return 'g';
            default:
                return 'N';
        }
    });
    return seq2;
}

UniqueBarcodes makeUnique(const std::unordered_set<std::string> &codes, uint16_t minlength)
{
    std::unordered_map<std::string, std::vector<std::string>> uniqueCodes;
    for (const auto &c : codes) {
        std::vector<std::string> u;
        u.push_back(c);
        if (minlength) {
            for (std::remove_reference<decltype(c)>::type::size_type i = 1; i < c.size() && c.size() - i >= minlength; ++i) {
                bool found = false;
                for (const auto &cc : codes) {
                    if (cc != c && cc.find(c.c_str() + i) != std::remove_reference<decltype(cc)>::type::npos) {
                        found = true;
                        break;
                    }
                }
                if (!found)
                    u.push_back(c.substr(i));
            }
        }
        uniqueCodes[c] = std::move(u);
    }
    return uniqueCodes;
}

namespace std
{
    size_t hash<std::pair<std::string, std::string>>::operator()(const std::pair<std::string, std::string> &p) const
    {
        auto h1 = std::hash<std::string>()(p.first);
        auto h2 = std::hash<std::string>()(p.second);
        return h1 ^ (h2 << 1);
    }
}
