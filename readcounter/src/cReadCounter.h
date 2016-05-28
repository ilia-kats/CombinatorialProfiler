#ifndef READCOUNTER_H
#define READCOUNTER_H

#include <unordered_map>
#include <string>
#include <list>

class BarcodeSet
{
public:
    BarcodeSet(const std::unordered_map<std::string, std::list<std::string>>&, const std::unordered_map<std::string, std::list<std::string>>&);

    const std::unordered_map<std::string, std::list<std::string>> fw;
    const std::unordered_map<std::string, std::list<std::string>> rev;

};


#endif
