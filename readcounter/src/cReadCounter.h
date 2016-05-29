#ifndef READCOUNTER_H
#define READCOUNTER_H

#include <unordered_map>
#include <functional>
#include <string>
#include <vector>

namespace std
{
template<> struct hash<std::pair<std::string, std::string>>
{
    size_t operator()(const std::pair<std::string, std::string> &p) const;
};
}

class BarcodeSet
{
public:
    BarcodeSet(const std::unordered_map<std::string, std::vector<std::string>>&, const std::unordered_map<std::string, std::vector<std::string>>&);

    const std::unordered_map<std::string, std::vector<std::string>> fw;
    const std::unordered_map<std::string, std::vector<std::string>> rev;

};

class Read
{
public:
    Read();
    Read(const std::string &name);
    Read(const std::string &name, const std::string &sequence);
    Read(const std::string &name, const std::string &sequence, const std::string &description);
    Read(const std::string &name, const std::string &sequence, const std::string &description, const std::string &quality);

    void setName(const std::string&);
    const std::string& getName() const;

    void setSequence(const std::string&);
    const std::string& getSequence() const;

    void setDescription(const std::string&);
    const std::string& getDescription() const;

    void setQuality(const std::string&);
    const std::string& getQuality() const;

    Read reverseComplement() const;

private:
    std::string m_name;
    std::string m_sequence;
    std::string m_description;
    std::string m_quality;
};

class SequenceMatcher
{
public:
    SequenceMatcher(const std::string&);

    std::pair<std::string::size_type, std::string::size_type> match(Read&, uint16_t) const;

private:
    const std::string m_fullseq;
    std::string m_upstreamseq;
    std::string m_downstreamseq;
    uint16_t m_insertlength;

    static std::pair<std::string::size_type, uint16_t> fuzzy_find(const std::string&, const std::string&);
};

class ReadCounter
{
public:
    ReadCounter(const std::string&, BarcodeSet *fw=nullptr, BarcodeSet *rev=nullptr);
    void countReads(const std::string&, int threads=1);
    const std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>>& getCounts() const;

private:
    class ThreadSynchronization;
    std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>> m_counts;

    SequenceMatcher m_matcher;
    BarcodeSet *m_barcodes_fw;
    BarcodeSet *m_barcodes_rev;

    void readFile(const std::string&, ThreadSynchronization*);
    void matchRead(ThreadSynchronization*);

    static const std::string dummykey;
};

#endif
