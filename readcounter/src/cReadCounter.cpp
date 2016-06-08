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
    Node() : experiment(nullptr) {}
    Node(std::string seq)
    : experiment(nullptr), sequence(std::move(seq)) {}

    virtual ~Node() {}

    virtual bool match(Read&, std::string *s=nullptr) const = 0;

    Experiment *experiment;

    std::vector<Node*> children;
    std::string sequence;

protected:
    static std::pair<std::string::size_type, uint16_t> fuzzy_find(const std::string &needle, const std::string &haystack)
    {
        auto totest = haystack.size() - needle.size();
        std::vector<uint16_t> mismatches(totest);
        auto mit = mismatches.begin();
        for (size_t i = 0; i < totest; ++i, ++mit) {
            *mit = std::inner_product(needle.cbegin(),
                                      needle.cend(),
                                      haystack.cbegin() + i,
                                      static_cast<decltype(mismatches)::value_type>(0),
                                      std::plus<decltype(mismatches)::value_type>(),
                                      [](const std::remove_reference<decltype(needle)>::type::value_type &n, const std::remove_reference<decltype(haystack)>::type::value_type &h) -> bool {return static_cast<bool>(n ^ h);});
        }
        auto bestmatch = std::min_element(mismatches.cbegin(), mismatches.cend());
        return std::make_pair(std::distance(mismatches.cbegin(), bestmatch), *bestmatch);
    }
};


class BarcodeNode : public Node
{
public:
    BarcodeNode() : Node(dummykey) {}
    virtual ~BarcodeNode() {}

    virtual bool match(Read &rd, std::string *s=nullptr) const
    {
        return true;
    }

protected:
    BarcodeNode(std::string seq, std::vector<std::string> uniqueseq)
    : Node(std::move(seq)), m_uniqueSequences(std::move(uniqueseq))
    {}

    std::vector<std::string> m_uniqueSequences;
};

class RevBarcodeNode : public BarcodeNode
{
public:
    RevBarcodeNode() : BarcodeNode() {}
    RevBarcodeNode(std::string fullseq, std::vector<std::string> uniqueSeqs)
    : BarcodeNode(std::move(fullseq), std::move(uniqueSeqs))
    {
        for (auto &s : m_uniqueSequences) {
            s = revCompl(s);
        }
    }

    virtual bool match(Read &rd, std::string *s=nullptr) const
    {
        auto rdseq = rd.getSequence();
        for (const auto &seq : m_uniqueSequences) {
            if (!rdseq.compare(rdseq.size() - seq.size(), seq.size(), seq)) {
                return true;
            }
        }
        return false;
    }
};

class FwBarcodeNode : public BarcodeNode
{
public:
    FwBarcodeNode() : BarcodeNode() {}
    FwBarcodeNode(std::string fullseq, std::vector<std::string> uniqueSeqs)
    : BarcodeNode(std::move(fullseq), std::move(uniqueSeqs))
    {}

    virtual bool match(Read &rd, std::string *s=nullptr) const
    {
        for (const auto &seq : m_uniqueSequences) {
            auto rdseq = rd.getSequence();
            if (!rdseq.compare(0, seq.size(), seq)) {
                return true;
            }
        }
        return false;
    }
};

class InsertNode : public Node
{
public:
    InsertNode(std::string seq, uint16_t mismatches = 1)
    : Node(std::move(seq)), m_mismatches(mismatches)
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

    bool match(Read &read, std::string *insert=nullptr) const
    {
        auto upstream = fuzzy_find(m_upstreamseq, read.getSequence());
        if (upstream.second > m_mismatches) {
            read = read.reverseComplement();
            upstream = fuzzy_find(m_upstreamseq, read.getSequence());
            if (upstream.second > m_mismatches)
                return false;
        }
        auto downstream = fuzzy_find(m_downstreamseq, read.getSequence());
        if (downstream.second > m_mismatches)
            return false;
        if (downstream.first - (upstream.first + m_upstreamseq.size()) != m_insertlength)
            return false;
        auto start = upstream.first + m_upstreamseq.size();

        if (insert) {
            insert->clear();
            auto it = read.getSequence().cbegin();
            std::copy(it + start, it + downstream.first, std::inserter(*insert, insert->begin()));
        }
        return true;
    }
private:
    uint16_t m_mismatches;

    std::string m_upstreamseq;
    std::string m_downstreamseq;
    uint16_t m_insertlength;
};

std::unordered_map<std::string, std::vector<std::string>> makeUnique(const std::unordered_set<std::string> &codes)
{
    std::unordered_map<std::string, std::vector<std::string>> uniqueCodes;
    for (const auto &c : codes) {
        std::vector<std::string> u;
        for (std::remove_reference<decltype(c)>::type::size_type i = 0; i < c.size(); ++i) {
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
        uniqueCodes[c] = std::move(u);
    }
    return uniqueCodes;
}

ReadCounter::ReadCounter(std::vector<Experiment*> experiments)
: m_read(0), m_counted(0), m_unmatchedTotal(0), m_unmatchedInsert(0), m_unmatchedBarcodes(0), m_unmatchedInsertSequence(0), m_written(0), m_experiments(experiments)
{
    std::unordered_set<std::string> fwCodes;
    std::unordered_set<std::string> revCodes;

    std::unordered_map<std::string, InsertNode*> inserts;
    for (const auto &exp : m_experiments) {
        if (!inserts.count(exp->insert)) {
            InsertNode *n = new InsertNode(exp->insert);
            inserts[exp->insert] = n;
            m_nodes.push_back(n);
            m_tree.push_back(n);
        }
        for (const auto &c : exp->fwBarcodeSet)
            fwCodes.insert(c.first);
        for (const auto &c : exp->revBarcodeSet)
            revCodes.insert(c.first);
    }

    std::unordered_map<std::string, std::vector<std::string>> uniqueFwCodes = makeUnique(fwCodes);
    std::unordered_map<std::string, std::vector<std::string>> uniqueRevCodes = makeUnique(revCodes);

    std::unordered_map<std::string, std::unordered_map<std::string, BarcodeNode*>> fwNodes;
    std::unordered_map<std::string,BarcodeNode*> dummyNodes;

    for (const auto &exp : m_experiments) {
        std::vector<BarcodeNode*> revNodes;
        for (const auto &c : exp->revBarcodeSet) {
            BarcodeNode *n = new RevBarcodeNode(c.first, uniqueRevCodes[c.first]);
            n->experiment = exp;
            revNodes.push_back(n);
            m_nodes.push_back(n);
        }

        for (const auto &c : exp->fwBarcodeSet) {
            if (!fwNodes[exp->insert].count(c.first)) {
                BarcodeNode *n = new FwBarcodeNode(c.first, uniqueFwCodes[c.first]);
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
        std::vector<const Node*> nodes;

        FailedMatch(Read r, Fail f, std::vector<const Node*> ns)
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

void ReadCounter::readFile(const std::string &file, ThreadSynchronization *sync)
{
    auto in = std::ifstream(file);
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
                fname.push_back(n->sequence);
        }
        if (fmatch.fail == ThreadSynchronization::FailedMatch::Fail::insertMatchFailed) {
            fname.push_back("noNamedInsertForExperiment");
            fname.push_back(fmatch.nodes.back()->experiment->name);
        }
        auto it = fname.cbegin();
        std::string outfname = *it++;
        for (; it != fname.cend(); ++it)
            outfname.append("_").append(*it);

        outfname.append(".fastq");
        if (!files[outfname].is_open())
            files[outfname].open(outfname);
        files[outfname] << fmatch.read.getName() << std::endl << fmatch.read.getSequence() << std::endl << fmatch.read.getDescription() << std::endl << fmatch.read.getQuality() << std::endl;
        ++m_written;
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
        std::vector<const Node*> nodes;

        std::string ins;
        const Node *cn = nullptr;
        for (const auto &n : m_tree) {
            if (n->match(rd, &ins)) {
                cn = n;
            }
        }
        if (cn) {
            nodes.push_back(cn);
            while (!cn->experiment && cn->children.size()) {
                bool found = false;
                for (const auto &n : cn->children) {
                    if (n->match(rd)) {
                        cn = n;
                        nodes.push_back(n);
                        found = true;
                        break;
                    }
                }
                if (!found)
                    break;
            }
            if (cn->experiment) {
                if (nodes.size() != 3)
                    throw std::logic_error("Unexpected number of tree levels: " + std::to_string(nodes.size()));

                // using exact sequence matching at the moment, maybe add option for mismatches?
                if (!cn->experiment->namedInserts.size() || cn->experiment->namedInserts.count(ins))
                    ++localcounts[cn->experiment][std::make_pair(cn->experiment->fwBarcodeSet[nodes[1]->sequence], cn->experiment->revBarcodeSet[nodes[2]->sequence])][ins];
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
