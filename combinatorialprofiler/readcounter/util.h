#ifndef UTIL_H
#define UTIL_H
#include <string>
std::string revCompl(const std::string &);

namespace std
{
    template<> struct hash<std::pair<std::string, std::string>>
    {
        size_t operator()(const std::pair<std::string, std::string> &p) const;
    };
}

#endif
