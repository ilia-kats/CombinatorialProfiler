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

static std::string dummykey = "";

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

Experiment::Experiment()
: ndsi(NDSIS::noNDSI)
{}

Experiment::Experiment(std::string n)
: name(std::move(n)), ndsi(NDSIS::noNDSI)
{}

class Node
{
public:
    Node() {}
    virtual bool match(Read&, std::string*) const = 0;

    std::shared_ptr<Experiment> experiment;

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
    BarcodeNode() : Node() {};
    BarcodeNode(std::string fullseq, std::vector<std::string> uniqueSeqs)
    : Node(), fullSequence(std::move(fullseq)), m_uniqueSequences(std::move(uniqueSeqs))
    {}
    std::string fullSequence;

protected:
    std::vector<std::string> m_uniqueSequences;
};

class RevBarcodeNode : public BarcodeNode
{
public:
    RevBarcodeNode() : BarcodeNode() {}
    RevBarcodeNode(std::string fullseq, std::vector<std::string> uniqueSeqs)
    : BarcodeNode(std::move(fullseq), std::move(uniqueSeqs))
    {
        // make reverse complement of unique sequences
    }

    virtual bool match(Read &rd, std::string*)
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

    virtual bool match(Read &rd, std::string*)
    {
        for (const auto &seq : m_uniqueSequences) {
            auto rdseq = rd.getSequence();
            if (!rdseq.compare(0, seq.size(), seq)) {
                return true;
            }
        }
        return false;
    }

    std::vector<RevBarcodeNode*> revNodes;
};

class InsertNode : public Node
{
public:
    InsertNode(const std::string &sequence, uint16_t mismatches = 1)
    : Node(), m_fullseq(sequence), m_mismatches(mismatches)
    {
        // Currently assuming there is only one insert sequence
        auto insert_start = m_fullseq.find_first_of("Nn");
        auto insert_end = m_fullseq.find_last_of("Nn");
        if (insert_start == std::string::npos || insert_end == std::string::npos)
            throw std::invalid_argument("No dynamic insert sequence found");
        m_insertlength = insert_end - insert_start + 1;
        m_upstreamseq = m_fullseq.substr(0, insert_start);
        m_downstreamseq = m_fullseq.substr(insert_end + 1, std::string::npos);
    }

    bool match(Read &read, std::string *insert) const
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

    std::vector<FwBarcodeNode*> fwNodes;
    std::vector<RevBarcodeNode*> revNodes;

private:
    const std::string m_fullseq;
    uint16_t m_mismatches;

    std::string m_upstreamseq;
    std::string m_downstreamseq;
    uint16_t m_insertlength;
};

Read::Read()
{}

Read::Read(std::string name)
: m_name(std::move(name))
{}

Read::Read(std::string name, std::string sequence)
: m_name(std::move(name)), m_sequence(std::move(sequence))
{}

Read::Read(std::string name, std::string sequence, std::string description)
: m_name(std::move(name)), m_sequence(std::move(sequence)), m_description(std::move(description))
{}

Read::Read(std::string name, std::string sequence, std::string description, std::string quality)
: m_name(std::move(name)), m_sequence(std::move(sequence)), m_description(std::move(description)), m_quality(std::move(quality))
{}

void Read::setName(std::string name)
{
    m_name = std::move(name);
}

const std::string& Read::getName() const
{
    return m_name;
}

void Read::setSequence(std::string sequence)
{
    m_sequence = std::move(sequence);
}

const std::string& Read::getSequence() const
{
    return m_sequence;
}

void Read::setDescription(std::string description)
{
    m_description = std::move(description);
}

const std::string& Read::getDescription() const
{
    return m_description;
}

void Read::setQuality(std::string quality)
{
    m_quality = std::move(quality);
}

const std::string& Read::getQuality() const
{
    return m_quality;
}

Read Read::reverseComplement() const
{
    std::string seq2;
    std::string qual;
    seq2.reserve(m_sequence.size());
    qual.reserve(m_quality.size());
    std::transform(m_sequence.crbegin(), m_sequence.crend(), std::inserter(seq2, seq2.begin()), [](const decltype(m_sequence)::value_type &c) -> decltype(m_sequence)::value_type {
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
    std::copy(m_quality.crbegin(), m_quality.crend(), std::inserter(qual, qual.begin()));
    return Read(m_name, seq2, m_description, qual);
}

ReadCounter::ReadCounter(std::vector<std::shared_ptr<Experiment>> experiments)
: m_read(0), m_counted(0), m_written(0), m_experiments(experiments)
{
}

struct ReadCounter::ThreadSynchronization
{
    struct FailedMatch
    {
        enum Fail {noFail = 0x0, insertFailed = 0x1, fwBarcodeFailed = 0x2, revBarcodeFailed = 0x4, insertMatchFailed = 0x8};
        Read read;
        int fail;
        std::string insert_name;
        std::string barcode_fw;
        std::string barcode_rev;

        FailedMatch(Read r, int f, std::string fw, std::string rev)
        : read(std::move(r)), fail(f), barcode_fw(std::move(fw)), barcode_rev(std::move(rev))
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
    assert(m_read == m_counted + m_unmatched_total);
    assert(m_unmatched_total == m_written);
}

uint64_t ReadCounter::read() const
{
    return m_read;
}

uint64_t ReadCounter::counted() const
{
    return m_counted;
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
        if (fmatch.fail == ThreadSynchronization::FailedMatch::Fail::insertMatchFailed || fmatch.fail & ThreadSynchronization::FailedMatch::Fail::revBarcodeFailed && !(fmatch.fail & ThreadSynchronization::FailedMatch::Fail::fwBarcodeFailed)) {
            fname.push_back("fw");
            fname.push_back(fmatch.barcode_fw);
        }
        if (fmatch.fail == ThreadSynchronization::FailedMatch::Fail::insertMatchFailed || fmatch.fail & ThreadSynchronization::FailedMatch::Fail::fwBarcodeFailed && !(fmatch.fail & ThreadSynchronization::FailedMatch::Fail::revBarcodeFailed)) {
            fname.push_back("rev");
            fname.push_back(fmatch.barcode_rev);
        }
        if (fmatch.fail == ThreadSynchronization::FailedMatch::Fail::insertMatchFailed)
            fname.push_back("unmatched");
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
    decltype(m_counts) localcounts;

//     while (true) {
//         qlock.lock();
//         if (!sync->inqueue.size() && !sync->eof)
//             sync->inqueueempty.wait(qlock, [&sync]{return sync->inqueue.size() || sync->eof;});
//         if (sync->eof && !sync->inqueue.size()) {
//             qlock.unlock();
//             break;
//         }
//         Read rd(sync->inqueue.front());
//         sync->inqueue.pop();
//         qlock.unlock();
//         sync->inqueuefull.notify_one();
//
//         int fail = ThreadSynchronization::FailedMatch::Fail::noFail;
//         std::string fwkey;
//         std::string revkey;
//         std::string insert;
//         for (const auto &i : m_matcher) {
//             try {
//                 auto insertmatch = i.second.match(rd, 1);
//                 bool fwfound = true;
//                 bool revfound = true;
//                 if (!m_barcodes_fw || !m_barcodes_fw->count(i.first))
//                     fwkey = dummykey;
//                 else {
//                     fwfound = false;
//                     for (const auto &code : m_barcodes_fw->at(i.first)) {
//                         for (const auto &seq : code.second) {
//                             if (!rd.getSequence().compare(0, seq.size(), seq)) {
//                                 fwfound = true;
//                                 fwkey = code.first;
//                                 break;
//                             }
//                         }
//                         if (fwfound)
//                             break;
//                     }
//                 }
//                 if (!m_barcodes_rev || !m_barcodes_rev->count(i.first))
//                     revkey = dummykey;
//                 else {
//                     revfound = false;
//                     for (const auto &code : m_barcodes_rev->at(i.first)) {
//                         for (const auto &seq : code.second) {
//                             auto rdseq = rd.getSequence();
//                             if (!rdseq.compare(rdseq.size() - seq.size(), seq.size(), seq)) {
//                                 revfound = true;
//                                 revkey = code.first;
//                                 break;
//                             }
//                         }
//                         if (revfound)
//                             break;
//                     }
//                 }
//
//                 if (fwfound && revfound) {
//                     auto seq = rd.getSequence().substr(insertmatch.first, insertmatch.second);
//                     if (!m_inserts || !m_inserts->count(i.first) || m_inserts->at(i.first).count(seq)) {
//                         ++localcounts[i.first][std::make_pair(fwkey, revkey)][seq];
//                         fail = ThreadSynchronization::FailedMatch::Fail::noFail;
//                     } else {
//                         fail = ThreadSynchronization::FailedMatch::Fail::insertMatchFailed;
//                     }
//                 } else {
//                     if (!fwfound) {
//                         ++unmatched_fw;
//                         fail |= ThreadSynchronization::FailedMatch::Fail::fwBarcodeFailed;
//                     }
//                     if (!revfound) {
//                         ++unmatched_rev;
//                         fail |= ThreadSynchronization::FailedMatch::Fail::revBarcodeFailed;
//                     }
//                 }
//                 break;
//             } catch (std::exception &e) {
//                 fail = ThreadSynchronization::FailedMatch::Fail::insertFailed;
//             }
//         }
//
//         if (fail == ThreadSynchronization::FailedMatch::Fail::insertFailed)
//             ++unmatched_insert;
//         if (fail != ThreadSynchronization::FailedMatch::Fail::noFail) {
//             ++unmatched_total;
//             oqlock.lock();
//             if (sync->outqueue.size() >= sync->maxsize)
//                 sync->outqueuefull.wait(oqlock, [&sync]{return sync->outqueue.size() < sync->maxsize;});
//             sync->outqueue.emplace(std::move(rd), fail, std::move(fwkey), std::move(revkey));
//             oqlock.unlock();
//             sync->outqueueempty.notify_one();
//         }
//     }
//
//     std::unique_lock<std::mutex> clock(sync->countsmutex);
//     for (const auto &insname : localcounts) {
//         for (const auto &codes : insname.second) {
//             for (const auto &insert : codes.second) {
//                 m_counts[insname.first][codes.first][insert.first] += insert.second;
//                 m_counted += insert.second;
//             }
//         }
//     }
//     clock.unlock();
}

const std::unordered_map<std::string, std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>>>& ReadCounter::getCounts() const
{
    return m_counts;
}
