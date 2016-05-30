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

BarcodeSet::BarcodeSet(const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> &f, const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> &r)
: fw(f), rev(r)
{}

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

const std::string ReadCounter::dummykey = std::string("");

SequenceMatcher::SequenceMatcher(const std::string &sequence)
: m_fullseq(sequence)
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

std::pair<std::string::size_type, uint16_t> SequenceMatcher::fuzzy_find(const std::string &needle, const std::string &haystack)
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

std::pair<std::string::size_type, std::string::size_type> SequenceMatcher::match(Read &read, uint16_t mismatches) const
{
    auto upstream = fuzzy_find(m_upstreamseq, read.getSequence());
    if (upstream.second > mismatches) {
        read = read.reverseComplement();
        upstream = fuzzy_find(m_upstreamseq, read.getSequence());
        if (upstream.second > mismatches)
            throw std::runtime_error("Upstream sequence not found in read");
    }
    auto downstream = fuzzy_find(m_downstreamseq, read.getSequence());
    if (downstream.second > mismatches)
        throw std::runtime_error("Downstream sequence not found in read");
    if (downstream.first - (upstream.first + m_upstreamseq.size()) != m_insertlength)
        throw std::runtime_error("Insert length doesn't match");
    auto start = upstream.first + m_upstreamseq.size();
    return std::make_pair(start, downstream.first - start);
}

ReadCounter::ReadCounter(const std::unordered_map<std::string, std::string> &insertseqs, BarcodeSet *fw, BarcodeSet *rev)
: m_read(0), m_counted(0), m_unmatched_insert(0), m_unmatched_nobarcode(0), m_unmatched_fw(0), m_unmatched_rev(0), m_barcodes_fw(fw), m_barcodes_rev(rev)
{
    for (const auto &ins : insertseqs) {
        m_matcher.emplace(ins);
    }
}

struct ReadCounter::ThreadSynchronization
{
    struct FailedMatch
    {
        enum class Fail {noFail, insertFailed, noBarcodeForInsert, fwBarcodeFailed, revBarcodeFailed};
        Read read;
        Fail fail;
        std::string insert_name;
        std::string barcode_fw;
        std::string barcode_rev;

        FailedMatch(Read r, Fail f, std::string fw, std::string rev)
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
}

uint64_t ReadCounter::read() const
{
    return m_read;
}

uint64_t ReadCounter::counted() const
{
    return m_counted;
}

uint64_t ReadCounter::unmatchedInsert() const
{
    return m_unmatched_insert;
}

uint16_t ReadCounter::insertsWithoutBarcodes() const
{
    return m_unmatched_nobarcode;
}

uint64_t ReadCounter::unmatchedBarcodeFw() const
{
    return m_unmatched_fw;
}

uint64_t ReadCounter::unmatchedBarcodeRev() const
{
    return m_unmatched_rev;
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
        std::string fname = prefix;
        switch (fmatch.fail) {
            case ThreadSynchronization::FailedMatch::Fail::insertFailed:
                fname.append("noinsert");
                break;
            case ThreadSynchronization::FailedMatch::Fail::noBarcodeForInsert:
                fname.append(fmatch.insert_name);
                break;
            case ThreadSynchronization::FailedMatch::Fail::fwBarcodeFailed:
            case ThreadSynchronization::FailedMatch::Fail::revBarcodeFailed:
                fname.append(fmatch.barcode_fw).append(fmatch.barcode_rev);
                break;
            default:
                break;
        }
        fname.append(".fastq");
        if (!files[fname].is_open())
            files[fname].open(fname);
        files[fname] << fmatch.read.getName() << std::endl << fmatch.read.getSequence() << std::endl << fmatch.read.getDescription() << std::endl << fmatch.read.getQuality() << std::endl;
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

    uint64_t unmatched_insert = 0;
    uint64_t unmatched_nobarcode = 0;
    uint64_t unmatched_fw = 0;
    uint64_t unmatched_rev = 0;

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
        std::string fwkey;
        std::string revkey;
        std::string insert;
        for (const auto &i : m_matcher) {
            try {
                auto insertmatch = i.second.match(rd, 1);
                bool fwfound = true;
                bool revfound = true;
                if (m_barcodes_fw == nullptr)
                    fwkey = dummykey;
                else {
                    fwfound = false;
                    for (const auto &code : m_barcodes_fw->fw.at(i.first)) {
                        for (const auto &seq : code.second) {
                            if (!rd.getSequence().compare(0, seq.size(), seq)) {
                                fwfound = true;
                                fwkey = code.first;
                                break;
                            }
                        }
                        if (fwfound)
                            break;
                    }
                }
                if (m_barcodes_rev == nullptr)
                    revkey = dummykey;
                else {
                    revfound = false;
                    for (const auto &code : m_barcodes_rev->rev.at(i.first)) {
                        for (const auto &seq : code.second) {
                            auto rdseq = rd.getSequence();
                            if (!rdseq.compare(rdseq.size() - seq.size(), seq.size(), seq)) {
                                revfound = true;
                                revkey = code.first;
                                break;
                            }
                        }
                        if (revfound)
                            break;
                    }
                }

                if (fwfound && revfound) {
                    ++localcounts[i.first][std::make_pair(fwkey, revkey)][rd.getSequence().substr(insertmatch.first, insertmatch.second)];
                    fail = ThreadSynchronization::FailedMatch::Fail::noFail;
                }
                else if (!fwfound) {
                    ++unmatched_fw;
                    fail = ThreadSynchronization::FailedMatch::Fail::fwBarcodeFailed;
                } else {
                    ++unmatched_rev;
                    fail = ThreadSynchronization::FailedMatch::Fail::revBarcodeFailed;
                }
                break;
            } catch (std::out_of_range &e) {
                fail = ThreadSynchronization::FailedMatch::Fail::noBarcodeForInsert;
                ++unmatched_nobarcode;
                break;
            } catch (std::exception &e) {
                fail = ThreadSynchronization::FailedMatch::Fail::insertFailed;
            }
        }

        if (fail == ThreadSynchronization::FailedMatch::Fail::insertFailed)
            ++unmatched_insert;
        if (fail != ThreadSynchronization::FailedMatch::Fail::noFail) {
            oqlock.lock();
            if (sync->outqueue.size() >= sync->maxsize)
                sync->outqueuefull.wait(oqlock, [&sync]{return sync->outqueue.size() < sync->maxsize;});
            sync->outqueue.emplace(std::move(rd), fail, std::move(fwkey), std::move(revkey));
            oqlock.unlock();
            sync->outqueueempty.notify_one();
        }
    }

    std::unique_lock<std::mutex> clock(sync->countsmutex);
    for (const auto &insname : localcounts) {
        for (const auto &codes : insname.second) {
            for (const auto &insert : codes.second) {
                m_counts[insname.first][codes.first][insert.first] += insert.second;
                m_counted += insert.second;
            }
        }
    }
    m_unmatched_insert += unmatched_insert;
    m_unmatched_nobarcode += unmatched_nobarcode;
    m_unmatched_fw += unmatched_fw;
    m_unmatched_rev += unmatched_rev;
    clock.unlock();
}

const std::unordered_map<std::string, std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>>>& ReadCounter::getCounts() const
{
    return m_counts;
}
