#include "cReadCounter.h"

#include <algorithm>
#include <iterator>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <numeric>
#include <fstream>
#include <condition_variable>
#include <atomic>
#include <thread>
#include <cmath>

#include <cassert>

static std::string dummykey = std::string();

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

class Read
{
public:
    Read() {}
    Read(std::string name)
    : m_name(std::move(name))
    {}
    Read(std::string name, std::string sequence)
    : m_name(std::move(name)), m_sequence(std::move(sequence))
    {}
    Read(std::string name, std::string sequence, std::string description)
    : m_name(std::move(name)), m_sequence(std::move(sequence)), m_description(std::move(description))
    {}
    Read(std::string name, std::string sequence, std::string description, std::string quality)
    : m_name(std::move(name)), m_sequence(std::move(sequence)), m_description(std::move(description)), m_quality(std::move(quality))
    {}

    void setName(std::string name)
    {
        m_name = std::move(name);
    }
    const std::string& getName() const
    {
        return m_name;
    }

    void setSequence(std::string sequence)
    {
        m_sequence = std::move(sequence);
    }
    const std::string& getSequence() const
    {
        return m_sequence;
    }

    void setDescription(std::string description)
    {
        m_description = std::move(description);
    }
    const std::string& getDescription() const
    {
        return m_description;
    }

    void setQuality(std::string quality)
    {
        m_quality = std::move(quality);
    }
    const std::string& getQuality() const
    {
        return m_quality;
    }

    Read reverseComplement() const
    {
        std::string qual;
        qual.reserve(m_quality.size());
        std::copy(m_quality.crbegin(), m_quality.crend(), std::inserter(qual, qual.begin()));
        return Read(m_name, revCompl(m_sequence), m_description, qual);
    }

private:
    std::string m_name;
    std::string m_sequence;
    std::string m_description;
    std::string m_quality;
};

Experiment::Experiment()
: ndsi(NDSIS::noNDSI)
{}

Experiment::Experiment(std::string n)
: name(std::move(n)), ndsi(NDSIS::noNDSI)
{}

class Node
{
public:
class MatchBase
{
public:
    MatchBase(const Node *n)
    : m_node(n), m_match(true) {}
    MatchBase(const Node *n, bool m)
    : m_node(n), m_match(m) {}
    virtual ~MatchBase() {}

    operator bool() const
    {
        return m_match;
    }

    const Node* node()
    {
        return m_node;
    }

protected:
    const Node *m_node;
    bool m_match;
};

template<typename T>
class Match : public MatchBase
{
public:
    virtual ~Match() {}

    friend bool operator<(const Match<T> &l, const Match<T> &r)
    {
        if (l and !r)
            return true;
        else if (r and !l)
            return false;
        else
            return l.m_mismatches < r.m_mismatches;
    }
    friend bool operator>(const Match &l, const Match &r)
    {
        return r < l;
    }
    friend bool operator<=(const Match &l, const Match &r)
    {
        return !(l > r);
    }
    friend bool operator>=(const Match &l, const Match &r)
    {
        return !(l < r);
    }
    friend bool operator==(const Match &l, const Match &r)
    {
        return l.match && r.match && l.m_calcBackMismatches == r.m_calcBackMismatches;
    }
    friend bool operator!=(const Match &l, const Match &r)
    {
        return !(l == r);
    }

protected:
    Match(const Node *n, bool m)
    : MatchBase(n, m), m_mismatches(0) {}

    Match(const Node *n, bool m, T mismatches)
    : MatchBase(n, m), m_mismatches(mismatches) {}

    T m_mismatches;
};

    Node() : experiment(nullptr) {}
    Node(std::string seq, std::string::size_type mismatches = 0)
    : experiment(nullptr), sequence(std::move(seq)), m_allowedMismatches(mismatches) {}

    virtual ~Node() {}

    virtual MatchBase* match(Read&) const = 0;

    Experiment *experiment;

    std::vector<Node*> children;
    std::string sequence;

protected:
    std::string::size_type m_allowedMismatches;

    template<class NeedleIt, class HaystackIt>
    static std::pair<std::string::size_type, std::string::size_type> fuzzy_find(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin, HaystackIt hend)
    {
        auto hsize = std::distance(hbegin, hend);
        auto nsize = std::distance(nbegin, nend);
        auto totest = hsize - nsize + 1;
        if (totest <= 0)
            return std::make_pair(0, UINT16_MAX);
        std::pair<std::string::size_type, std::string::size_type> bestMatch(0, SIZE_MAX);
        for (std::string::size_type i = 0; i < totest; ++i) {
            auto mm = std::inner_product(nbegin,
                                         nend,
                                         hbegin + i,
                                         static_cast<std::string::size_type>(0),
                                         std::plus<std::string::size_type>(),
                                         [](const typename decltype(nbegin)::value_type &n, const typename decltype(hbegin)::value_type &h) -> bool {return static_cast<bool>(n ^ h);});
            if (mm < bestMatch.second) {
                bestMatch.first = i;
                bestMatch.second = mm;
            }
            if (!mm)
                break;
        }
        return bestMatch;
    }
};

class BarcodeNode : public Node
{
public:
class BarcodeMatch : public Match<float>
{
public:
    BarcodeMatch(const BarcodeNode *n, std::string::size_type length, std::string::size_type mismatches)
    : Match(n, false, (float)mismatches * (float)n->sequence.size() / (float)length), m_matchedLength(length), m_actualMismatches(mismatches)
    {
        if (m_mismatches > n->m_allowedMismatches)
            m_match = false;
        else
            m_match = true;
    }

private:
    std::string::size_type m_matchedLength;
    std::string::size_type m_actualMismatches;
};

    BarcodeNode() : Node(dummykey, 0) {}
    virtual ~BarcodeNode() {}

    virtual BarcodeMatch* match(Read &rd) const
    {
        return new BarcodeMatch(this, SIZE_MAX, 0);
    }

protected:
    BarcodeNode(std::string seq, std::vector<std::string> uniqueseq, std::string::size_type mismatches=0)
    : Node(std::move(seq), mismatches), m_uniqueSequences(std::move(uniqueseq))
    {}

    std::vector<std::string> m_uniqueSequences;
};

class RevBarcodeNode : public BarcodeNode
{
public:
    RevBarcodeNode() : BarcodeNode() {}
    RevBarcodeNode(std::string fullseq, std::vector<std::string> uniqueSeqs, std::string::size_type mismatches=0)
    : BarcodeNode(std::move(fullseq), std::move(uniqueSeqs), mismatches)
    {
        for (auto &s : m_uniqueSequences) {
            s = revCompl(s);
        }
    }

    virtual BarcodeMatch* match(Read &rd) const
    {
        auto rdseq = rd.getSequence();
        std::pair<decltype(m_uniqueSequences)::size_type, std::string::size_type> bestFind(SIZE_MAX, SIZE_MAX);
        if (!m_allowedMismatches) {
            for (const auto &seq : m_uniqueSequences) {
                if (!rdseq.compare(rdseq.size() - seq.size(), seq.size(), seq)) {
                    return new BarcodeMatch(this, seq.size(), 0);
                }
            }
        } else {
            auto it = m_uniqueSequences.cbegin();
            for (decltype(m_uniqueSequences)::size_type i = 0; i < m_uniqueSequences.size(); ++i, ++it) {
                auto found = fuzzy_find(it->cbegin(), it->cend(), rdseq.cend() - it->size(), rdseq.cend());
                if (found.second < bestFind.second) {
                    bestFind.first = i;
                    bestFind.second = found.second;
                }
                if (!found.second)
                    break;
            }
        }
        return new BarcodeMatch(this, m_uniqueSequences[bestFind.first].size(), bestFind.second);
    }
};

class FwBarcodeNode : public BarcodeNode
{
public:
    FwBarcodeNode() : BarcodeNode() {}
    FwBarcodeNode(std::string fullseq, std::vector<std::string> uniqueSeqs, std::string::size_type mismatches=0)
    : BarcodeNode(std::move(fullseq), std::move(uniqueSeqs), mismatches)
    {}

    virtual BarcodeMatch* match(Read &rd) const
    {
        auto rdseq = rd.getSequence();
        std::pair<decltype(m_uniqueSequences)::size_type, std::string::size_type> bestFind(SIZE_MAX, SIZE_MAX);
        if (!m_allowedMismatches){
            for (const auto &seq : m_uniqueSequences) {
                if (!rdseq.compare(0, seq.size(), seq)) {
                    return new BarcodeMatch(this, seq.size(), 0);
                }
            }
        } else {
            auto it = m_uniqueSequences.cbegin();
            for (decltype(m_uniqueSequences)::size_type i = 0; i < m_uniqueSequences.size(); ++i, ++it) {
                auto found = fuzzy_find(it->cbegin(), it->cend(), rdseq.cbegin(), rdseq.cbegin() + it->size());
                if (found.second < bestFind.second) {
                    bestFind.first = i;
                    bestFind.second = found.second;
                }
                if (!found.second)
                    break;
            }
        }
        return new BarcodeMatch(this, m_uniqueSequences[bestFind.first].size(), bestFind.second);
    }
};

class InsertNode : public Node
{
public:
class InsertMatch : public Match<std::string::size_type>
{
public:
    InsertMatch(const InsertNode *n, std::string::size_type mismatches, std::string insert = std::string())
    : Match(n, mismatches <= n->m_allowedMismatches, mismatches), m_insert(std::move(insert)) {}

    InsertMatch(const InsertNode *n, bool match)
    : Match(n, match) {}

    const std::string& insertSequence() const
    {
        return m_insert;
    }

private:
    std::string m_insert;
};

    InsertNode(std::string seq, std::string::size_type mismatches = 1)
    : Node(std::move(seq), mismatches)
    {
        // Currently assuming there is only one insert sequence
        auto insert_start = sequence.find_first_of("Nn");
        auto insert_end = sequence.find_last_of("Nn");
        if (insert_start == std::string::npos || insert_end == std::string::npos)
            throw std::invalid_argument("No dynamic insert sequence found");
        m_insertlength = insert_end - insert_start + 1;
        m_upstreamseq = sequence.substr(0, insert_start);
        m_downstreamseq = sequence.substr(insert_end + 1, std::string::npos);
    }

    virtual InsertMatch* match(Read &read) const
    {
        decltype(m_upstreamseq)::size_type start, end;
        std::string::size_type mismatches_cnt = 0;
        if (m_allowedMismatches > 0) {
            auto upstream = fuzzy_find(m_upstreamseq.cbegin(), m_upstreamseq.cend(), read.getSequence().cbegin(), read.getSequence().cend());
            if (upstream.second > m_allowedMismatches) {
                read = read.reverseComplement();
                upstream = fuzzy_find(m_upstreamseq.cbegin(), m_upstreamseq.cend(), read.getSequence().cbegin(), read.getSequence().cend());
                if (upstream.second > m_allowedMismatches)
                    return new InsertMatch(this, upstream.second);
            }
            start = upstream.first + m_upstreamseq.size();
            auto downstream = fuzzy_find(m_downstreamseq.cbegin(), m_downstreamseq.cend(), read.getSequence().cbegin() + start, read.getSequence().cend());
            downstream.first += start;
            mismatches_cnt = downstream.second + upstream.second;
            if (mismatches_cnt > m_allowedMismatches)
                return new InsertMatch(this, mismatches_cnt);
            end = downstream.first;
        } else {
            start = read.getSequence().find(m_upstreamseq);
            if (start == decltype(m_upstreamseq)::npos) {
                read = read.reverseComplement();
                start = read.getSequence().find(m_upstreamseq);
                if (start == decltype(m_upstreamseq)::npos)
                    return new InsertMatch(this, false);
            }
            start += m_upstreamseq.size();
            end = read.getSequence().find(m_downstreamseq, start);
            if (end == decltype(m_downstreamseq)::npos)
                return new InsertMatch(this, false);
        }

        if (end - start != m_insertlength)
            return new InsertMatch(this, false);

        std::string insert;
        auto it = read.getSequence().cbegin();
        std::copy(it + start, it + end, std::inserter(insert, insert.begin()));
        return new InsertMatch(this, mismatches_cnt, std::move(insert));
    }

private:
    std::string m_upstreamseq;
    std::string m_downstreamseq;
    std::string::size_type m_insertlength;
};

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

ReadCounter::ReadCounter(std::vector<Experiment*> experiments, uint16_t insert_mismatches, uint16_t unique_barcode_length, uint16_t allowed_barcode_mismatches)
: m_allowedMismatches(insert_mismatches), m_uniqueBarcodeLength(unique_barcode_length), m_allowedBarcodeMismatches(allowed_barcode_mismatches), m_read(0), m_counted(0), m_unmatchedTotal(0), m_unmatchedInsert(0), m_unmatchedBarcodes(0), m_unmatchedInsertSequence(0), m_written(0), m_experiments(experiments)
{
    std::unordered_map<std::string, std::unordered_set<std::string>> fwCodes;
    std::unordered_map<std::string, std::unordered_set<std::string>> revCodes;

    std::unordered_map<std::string, InsertNode*> inserts;
    for (const auto &exp : m_experiments) {
        if (!inserts.count(exp->insert)) {
            InsertNode *n = new InsertNode(exp->insert, m_allowedMismatches);
            inserts[exp->insert] = n;
            m_nodes.push_back(n);
            m_tree.push_back(n);
        }
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

    std::unordered_map<std::string, std::unordered_map<std::string, BarcodeNode*>> fwNodes;
    std::unordered_map<std::string,BarcodeNode*> dummyNodes;

    for (const auto &exp : m_experiments) {
        std::vector<BarcodeNode*> revNodes;
        for (const auto &c : exp->revBarcodeSet) {
            BarcodeNode *n = new RevBarcodeNode(c.first, m_uniqueRevCodes[exp->insert][c.first], m_allowedBarcodeMismatches);
            n->experiment = exp;
            revNodes.push_back(n);
            m_nodes.push_back(n);
        }

        for (const auto &c : exp->fwBarcodeSet) {
            if (!fwNodes[exp->insert].count(c.first)) {
                BarcodeNode *n = new FwBarcodeNode(c.first, m_uniqueFwCodes[exp->insert][c.first], m_allowedBarcodeMismatches);
                fwNodes[exp->insert][c.first] = n;
                m_nodes.push_back(n);
                inserts[exp->insert]->children.push_back(n);
            }
            std::copy(revNodes.cbegin(), revNodes.cend(), std::inserter(fwNodes[exp->insert][c.first]->children, fwNodes[exp->insert][c.first]->children.end()));
        }
        if (!exp->fwBarcodeSet.size()) {
            if (!dummyNodes.count(exp->insert)) {
                BarcodeNode *n = new BarcodeNode();
                m_nodes.push_back(n);
                dummyNodes[exp->insert] = n;
                std::copy(revNodes.cbegin(), revNodes.cend(), std::inserter(n->children, n->children.end()));
            }
        }
    }

    // have to make sure that dummy nodes are always last in the list
    for (const auto &d : dummyNodes) {
        inserts[d.first]->children.push_back(d.second);
    }
    for (const auto &exp : m_experiments) {
        if (!exp->revBarcodeSet.size()) {
            BarcodeNode *n = new BarcodeNode();
            n->experiment = exp;
            m_nodes.push_back(n);
            if (exp->fwBarcodeSet.size()) {
                for (const auto &c : exp->fwBarcodeSet) {
                    fwNodes[exp->insert][c.first]->children.push_back(n);
                }
            } else {
                dummyNodes[exp->insert]->children.push_back(n);
            }
        }
    }
}

ReadCounter::~ReadCounter()
{
    for (auto &n : m_nodes)
        delete n;
}

struct ReadCounter::ThreadSynchronization
{
    struct FailedMatch
    {
        enum class Fail {noFail, insertFailed, barcodeFailed, insertMatchFailed};
        Read read;
        Fail fail;
        std::vector<Node::MatchBase*> nodes;

        FailedMatch(Read r, Fail f, std::vector<Node::MatchBase*> ns)
        : read(std::move(r)), fail(f), nodes(std::move(ns))
        {}
    };
    static constexpr int maxsize = 500;
    std::queue<Read> inqueue;
    std::mutex inqueuemutex;
    std::atomic<bool> eof;
    std::condition_variable inqueuefull;
    std::condition_variable inqueueempty;

    std::atomic<bool> finishedMatching;

    std::queue<FailedMatch> outqueue;
    std::mutex outqueuemutex;
    std::condition_variable outqueuefull;
    std::condition_variable outqueueempty;

    std::mutex countsmutex;
};

void ReadCounter::countReads(const std::string &file, const std::string &outprefix, int threads)
{
    ThreadSynchronization sync;
    sync.finishedMatching = false;
    sync.eof = false; // need this to prevent race condition: worker threads start and exit before reader thread
    std::thread reader(&ReadCounter::readFile, this, file, &sync);
    std::thread writer(&ReadCounter::writeReads, this, outprefix, &sync);
    std::vector<std::thread> workers;
    for (int i = 0; i < threads; ++i)
        workers.emplace_back(&ReadCounter::matchRead, this, &sync);
    reader.join();
    for (auto &worker : workers)
        worker.join();
    sync.finishedMatching = true;
    sync.outqueueempty.notify_one();
    writer.join();
    assert(m_read == m_counted + m_unmatchedTotal);
    assert(m_unmatchedTotal = m_unmatchedInsert + m_unmatchedBarcodes + m_unmatchedInsertSequence);
    assert(m_unmatchedTotal == m_written);
}

uint16_t ReadCounter::allowedMismatches() const
{
    return m_allowedMismatches;
}

uint16_t ReadCounter::minimumUniqueBarcodeLength() const
{
    return m_uniqueBarcodeLength;
}

uint16_t ReadCounter::allowedBarcodeMismatches() const
{
    return m_allowedBarcodeMismatches;
}

uint64_t ReadCounter::read() const
{
    return m_read;
}

uint64_t ReadCounter::counted() const
{
    return m_counted;
}

uint64_t ReadCounter::unmatchedTotal() const
{
    return m_unmatchedTotal;
}

uint64_t ReadCounter::unmatchedInsert() const
{
    return m_unmatchedInsert;
}

uint64_t ReadCounter::unmatchedBarcodes() const
{
    return m_unmatchedBarcodes;
}

uint64_t ReadCounter::unmatchedInsertSequence() const
{
    return m_unmatchedInsertSequence;
}

uint64_t ReadCounter::written() const
{
    return m_written;
}

std::unordered_map<std::string, UniqueBarcodes> ReadCounter::uniqueForwardBarcodes() const
{
    return m_uniqueFwCodes;
}

std::unordered_map<std::string, UniqueBarcodes> ReadCounter::uniqueReverseBarcodes() const
{
    return m_uniqueRevCodes;
}

void ReadCounter::readFile(const std::string &file, ThreadSynchronization *sync)
{
    std::ifstream in(file);
    std::array<std::string, 4> read;
    std::unique_lock<std::mutex> qlock(sync->inqueuemutex, std::defer_lock);
    sync->eof = false;
    for (int i = 0; in; ++i) {
        auto ri = i % 4;
        std::getline(in, read[ri]);
        if (ri == 3) {
            qlock.lock();
            if (sync->inqueue.size() >= sync->maxsize)
                sync->inqueuefull.wait(qlock, [&sync]{return sync->inqueue.size() < sync->maxsize;});
            sync->inqueue.emplace(std::move(read[0]), std::move(read[1]), std::move(read[2]), std::move(read[3]));
            ++m_read;
            qlock.unlock();
            sync->inqueueempty.notify_one();
        }
    }
    sync->eof = true;
    sync->inqueueempty.notify_all();
}

void ReadCounter::writeReads(const std::string &prefix, ThreadSynchronization *sync)
{
    std::unordered_map<std::string, std::ofstream> files;
    std::unique_lock<std::mutex> qlock(sync->outqueuemutex, std::defer_lock);
    while (true) {
        qlock.lock();
        if (!sync->outqueue.size())
            sync->outqueueempty.wait(qlock, [&sync]{return sync->outqueue.size() || sync->finishedMatching;});
        if (!sync->outqueue.size() && sync->finishedMatching) {
            qlock.unlock();
            return;
        }
        auto fmatch = sync->outqueue.front();
        sync->outqueue.pop();
        qlock.unlock();
        sync->outqueuefull.notify_one();
        std::vector<std::string> fname;
        fname.push_back(prefix);
        if (fmatch.fail == ThreadSynchronization::FailedMatch::Fail::insertFailed)
            fname.push_back("noinsert");
        else if (fmatch.fail == ThreadSynchronization::FailedMatch::Fail::barcodeFailed) {
            fname.push_back("noBarcodeForInsert");
            for (const auto &n : fmatch.nodes)
                fname.push_back(n->node()->sequence);
        }
        if (fmatch.fail == ThreadSynchronization::FailedMatch::Fail::insertMatchFailed) {
            fname.push_back("noNamedInsertForExperiment");
            fname.push_back(fmatch.nodes.back()->node()->experiment->name);
        }
        auto it = fname.cbegin();
        std::string outfname = *it++;
        for (; it != fname.cend(); ++it)
            outfname.append("_").append(*it);

        outfname.append(".fastq");
        if (!files[outfname].is_open())
            files[outfname].open(outfname, std::ios_base::ate);
        files[outfname] << fmatch.read.getName() << std::endl << fmatch.read.getSequence() << std::endl << fmatch.read.getDescription() << std::endl << fmatch.read.getQuality() << std::endl;
        ++m_written;
        for (auto &n : fmatch.nodes)
            delete n;
    }
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

void ReadCounter::matchRead(ThreadSynchronization *sync)
{
    std::unique_lock<std::mutex> qlock(sync->inqueuemutex, std::defer_lock);
    std::unique_lock<std::mutex> oqlock(sync->outqueuemutex, std::defer_lock);

    std::unordered_map<Experiment*, Counts> localcounts;
    uint64_t unmatchedInsert = 0;
    uint64_t unmatchedBarcodes = 0;
    uint64_t unmatchedInsertSequence = 0;
    uint64_t unmatchedTotal = 0;
    while (true) {
        qlock.lock();
        if (!sync->inqueue.size() && !sync->eof)
            sync->inqueueempty.wait(qlock, [&sync]{return sync->inqueue.size() || sync->eof;});
        if (sync->eof && !sync->inqueue.size()) {
            qlock.unlock();
            break;
        }
        Read rd(sync->inqueue.front());
        sync->inqueue.pop();
        qlock.unlock();
        sync->inqueuefull.notify_one();

        ThreadSynchronization::FailedMatch::Fail fail = ThreadSynchronization::FailedMatch::Fail::noFail;
        std::vector<Node::MatchBase*> nodes;

        std::string ins;
        InsertNode::InsertMatch *in = nullptr;
        for (const auto &n : m_tree) {
            InsertNode::InsertMatch *cn = n->match(rd);
            if (*cn) {
                nodes.push_back(cn);
                in = cn;
                break;
            }
            else
                delete cn;
        }
        if (in && *in) {
            Node::MatchBase *cn = in;
            while (!cn->node()->experiment && cn->node()->children.size()) {
                Node::MatchBase *best = nullptr;
                for (const auto &n : cn->node()->children) {
                    Node::MatchBase *bn = n->match(rd);
                    if (*bn && best && *bn < *best) {
                        delete best;
                        best = bn;
                        break;
                    } else if (*bn && !best) {
                        best = bn;
                    } else
                        delete bn;
                }
                if (best && *best) {
                    nodes.push_back(best);
                    cn = best;
                } else {
                    if (best)
                        delete best;
                    break;
                }
            }
            if (cn->node()->experiment) {
                if (nodes.size() != 3) {
                    for (auto &n : nodes)
                        delete n;
                    throw std::logic_error("Unexpected number of tree levels: " + std::to_string(nodes.size()));
                }

                // using exact sequence matching at the moment, maybe add option for mismatches?
                if (!cn->node()->experiment->namedInserts.size() || cn->node()->experiment->namedInserts.count(in->insertSequence()))
                    ++localcounts[cn->node()->experiment][std::make_pair(cn->node()->experiment->fwBarcodeSet[nodes[1]->node()->sequence], cn->node()->experiment->revBarcodeSet[nodes[2]->node()->sequence])][in->insertSequence()];
                else {
                    fail = ThreadSynchronization::FailedMatch::Fail::insertMatchFailed;
                    ++unmatchedInsertSequence;
                }
            } else {
                fail = ThreadSynchronization::FailedMatch::Fail::barcodeFailed;
                ++unmatchedBarcodes;
            }
        } else {
            fail = ThreadSynchronization::FailedMatch::Fail::insertFailed;
            ++unmatchedInsert;
        }

        if (fail != ThreadSynchronization::FailedMatch::Fail::noFail) {
            ++unmatchedTotal;
            oqlock.lock();
            if (sync->outqueue.size() >= sync->maxsize)
                sync->outqueuefull.wait(oqlock, [&sync]{return sync->outqueue.size() < sync->maxsize;});
            sync->outqueue.emplace(std::move(rd), fail, std::move(nodes));
            oqlock.unlock();
            sync->outqueueempty.notify_one();
        } else {
            for (auto &n : nodes)
                delete n;
        }
    }

    std::unique_lock<std::mutex> clock(sync->countsmutex);
    for (const auto &exp : localcounts) {
        for (const auto &bcodes : exp.second) {
            for (const auto &seq : bcodes.second) {
                exp.first->counts[bcodes.first][seq.first] += seq.second;
                m_counted += seq.second;
            }
        }
    }
    m_unmatchedTotal += unmatchedTotal;
    m_unmatchedInsert += unmatchedInsert;
    m_unmatchedBarcodes += unmatchedBarcodes;
    m_unmatchedInsertSequence += unmatchedInsertSequence;
    clock.unlock();
}
