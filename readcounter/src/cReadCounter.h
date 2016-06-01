#ifndef READCOUNTER_H
#define READCOUNTER_H

#include <unordered_map>
#include <unordered_set>
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

typedef std::unordered_map<std::string, std::unordered_set<std::string>> InsertSet;

typedef std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> BarcodeSet;

class Read
{
public:
    Read();
    Read(std::string name);
    Read(std::string name, std::string sequence);
    Read(std::string name, std::string sequence, std::string description);
    Read(std::string name, std::string sequence, std::string description, std::string quality);

    void setName(std::string);
    const std::string& getName() const;

    void setSequence(std::string);
    const std::string& getSequence() const;

    void setDescription(std::string);
    const std::string& getDescription() const;

    void setQuality(std::string);
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
    ReadCounter(const std::unordered_map<std::string, std::string>&, BarcodeSet *fw=nullptr, BarcodeSet *rev=nullptr, InsertSet *ins = nullptr);
    void countReads(const std::string&, const std::string&, int threads=1);
    const std::unordered_map<std::string, std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>>>& getCounts() const;

    uint64_t read() const;
    uint64_t counted() const;
    uint64_t unmatchedInsert() const;
    uint64_t unmatchedBarcodeFw() const;
    uint64_t unmatchedBarcodeRev() const;

private:
    struct ThreadSynchronization;
    std::unordered_map<std::string, std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>>> m_counts;

    uint64_t m_read;
    uint64_t m_counted;
    uint64_t m_unmatched_insert;
    uint64_t m_unmatched_fw;
    uint16_t m_unmatched_rev;

    BarcodeSet *m_barcodes_fw;
    BarcodeSet *m_barcodes_rev;
    InsertSet *m_inserts;

    std::unordered_map<std::string, SequenceMatcher> m_matcher;

    void readFile(const std::string&, ThreadSynchronization*);
    void matchRead(ThreadSynchronization*);
    void writeReads(const std::string&, ThreadSynchronization*);

    static const std::string dummykey;
};

#endif
