#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <numeric>
#include <valarray>

#include <iostream>

#ifndef NDEBUG
template<typename T>
void test_print(std::valarray<T>& v, int rows, int cols)
{
    for(int r=0; r<rows; ++r) {
        for(int c=0; c<cols; ++c) {
            std::cout << v[rows * c + r] << ' ';
        }
        std::cout << '\n';
    }
}
#endif

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

template<class NeedleIt, class HaystackIt>
std::string::size_type hamming_distance(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin)
{
    return std::inner_product(nbegin,
                              nend,
                              hbegin,
                              static_cast<std::string::size_type>(0),
                              std::plus<std::string::size_type>(),
                              [](const typename decltype(nbegin)::value_type &n, const typename decltype(hbegin)::value_type &h) -> bool {return static_cast<bool>(n ^ h);});
}

template<class NeedleIt, class HaystackIt>
std::string::size_type seqlev_distance(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin)
{
    return seqlev_distance_begin(nbegin, nend, hbegin, hbegin + std::distance(nbegin, nend));
}

// Buschmann, Tilo, Bystrykh, Leonid V: Levenshtein error-correcting barcodes for multiplexed DNA sequencing, BMC bioinformatics 14(1), BioMed Central Ltd, 272, 2013
template<class NeedleIt, class HaystackIt>
std::string::size_type seqlev_distance(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin, HaystackIt hend)
{
    auto needleLength = std::distance(nbegin, nend) + 1;
    auto haystackLength = std::distance(hbegin, hend) + 1;
    std::valarray<std::string::size_type> distances(needleLength * haystackLength);
    {
        auto it = std::begin(distances);
        std::iota(it, it + needleLength, 0);
        for (std::size_t i = 0; i < haystackLength; ++i)
        {
            distances[i * needleLength] = i;
        }
    }

    for (std::size_t i = 1; i < needleLength; ++i, ++nbegin) {
        auto hit = hbegin;
        for (std::size_t j = 1; j < haystackLength; ++j, ++hit) {
            uint8_t cost = 0;
            if (*nbegin != *hit)
                cost = 1;
            auto cindex = (i - 1) + (j - 1) * needleLength;
            distances[i + j * needleLength] = std::min({
                distances[cindex] + cost,
                distances[cindex + 1] + 1,
                distances[cindex + needleLength] + 1
            });
        }
    }
    auto min_distance = *(std::end(distances) - 1);

#ifndef NDEBUG
    test_print(distances, needleLength, haystackLength);
    std::cout << min_distance << std::endl;
#endif

    for (std::size_t i = needleLength * (haystackLength - 1); i < distances.size(); ++i)
        min_distance = std::min(min_distance, distances[i]);
    for (std::size_t i = needleLength - 1; i < distances.size(); i += needleLength)
        min_distance = std::min(min_distance, distances[i]);
    return min_distance;
}

template<class NeedleIt, class HaystackIt>
std::string::size_type seqlev_distance_end(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin, HaystackIt hend)
{
    auto needleLength = std::distance(nbegin, nend) + 1;
    auto haystackLength = std::distance(hbegin, hend) + 1;
    std::valarray<std::string::size_type> distances(needleLength * haystackLength);
    {
        auto it = std::begin(distances);
        std::iota(it, it + needleLength, 0);
        for (std::size_t i = 0; i < haystackLength; ++i)
        {
            distances[i * needleLength] = i;
        }
    }

    --nend;
    --hend;
    for (std::size_t i = 1; i < needleLength; ++i, --nend) {
        auto hit = hend;
        for (std::size_t j = 1; j < haystackLength; ++j, --hit) {
            uint8_t cost = 0;
            if (*nend != *hit)
                cost = 1;
            auto cindex = (i - 1) + (j - 1) * needleLength;
            distances[i + j * needleLength] = std::min({
                distances[cindex] + cost,
                distances[cindex + 1] + 1,
                distances[cindex + needleLength] + 1
            });
        }
    }
    auto min_distance = *(std::end(distances) - 1);

#ifndef NDEBUG
    test_print(distances, needleLength, haystackLength);
    std::cout << min_distance << std::endl;
#endif

    for (std::size_t i = needleLength * (haystackLength - 1); i < distances.size(); ++i)
        min_distance = std::min(min_distance, distances[i]);
    for (std::size_t i = needleLength - 1; i < distances.size(); i += needleLength)
        min_distance = std::min(min_distance, distances[i]);
    return min_distance;
}

std::string::size_type seqlev_distance(const std::string&, const std::string&, bool at_end=false);

#endif
