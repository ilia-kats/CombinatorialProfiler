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

BarcodeSet::BarcodeSet(const std::unordered_map<std::string, std::vector<std::string>> &f, const std::unordered_map<std::string, std::vector<std::string>> &r)
: fw(f), rev(r)
{}

Read::Read()
{}

Read::Read(const std::string &name)
: m_name(name)
{}

Read::Read(const std::string &name, const std::string &sequence)
: m_name(name), m_sequence(sequence)
{}

Read::Read(const std::string &name, const std::string &sequence, const std::string &description)
: m_name(name), m_sequence(sequence), m_description(description)
{}

Read::Read(const std::string &name, const std::string &sequence, const std::string &description, const std::string &quality)
: m_name(name), m_sequence(sequence), m_description(description), m_quality(quality)
{}

void Read::setName(const std::string &name)
{
    m_name = name;
}

const std::string& Read::getName() const
{
    return m_name;
}

void Read::setSequence(const std::string &sequence)
{
    m_sequence = sequence;
}

const std::string& Read::getSequence() const
{
    return m_sequence;
}

void Read::setDescription(const std::string &description)
{
    m_description = description;
}

const std::string& Read::getDescription() const
{
    return m_description;
}

void Read::setQuality(const std::string &quality)
{
    m_quality = quality;
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

ReadCounter::ReadCounter(const std::string &insertseq, BarcodeSet *fw, BarcodeSet *rev)
: m_matcher(insertseq), m_barcodes_fw(fw), m_barcodes_rev(rev)
{}

class ReadCounter::ThreadSynchronization
{
public:
    static constexpr int maxsize = 500;
    std::queue<Read> queue;
    std::mutex queuemutex;
    std::mutex countsmutex;
    std::atomic<bool> eof;
    std::condition_variable queuefull;
    std::condition_variable queueempty;
};

void ReadCounter::countReads(const std::string &file, int threads)
{
    ThreadSynchronization sync;
    std::thread reader(&ReadCounter::readFile, this, file, &sync);
    std::vector<std::thread> workers;
    for (int i = 0; i < threads; ++i)
        workers.emplace_back(&ReadCounter::matchRead, this, &sync);
    reader.join();
    for (auto &worker : workers)
        worker.join();
}

void ReadCounter::readFile(const std::string &file, ThreadSynchronization *sync)
{
    auto in = std::ifstream(file);
    std::array<std::string, 4> read;
    std::unique_lock<std::mutex> qlock(sync->queuemutex, std::defer_lock);
    sync->eof = false;
    for (int i = 0; in; ++i) {
        auto ri = i % 4;
        std::getline(in, read[ri]);
        if (ri == 3) {
            qlock.lock();
            sync->queuefull.wait(qlock, [&sync]{return sync->queue.size() < sync->maxsize;});
            sync->queue.emplace(read[0], read[1], read[2], read[3]);
            qlock.unlock();
            sync->queueempty.notify_one();
        }
    }
    sync->queueempty.notify_all();
    sync->eof = true;
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
    std::unique_lock<std::mutex> qlock(sync->queuemutex, std::defer_lock);
    decltype(m_counts) localcounts;

    while (true) {
        qlock.lock();
        if (!sync->queue.size())
            sync->queueempty.wait(qlock, [&sync]{return sync->queue.size() || sync->eof;});
        if (sync->eof && !sync->queue.size()) {
            qlock.unlock();
            break;
        }
        Read rd(sync->queue.front());
        sync->queue.pop();
        qlock.unlock();
        sync->queuefull.notify_one();

        try {
            auto insertmatch = m_matcher.match(rd, 1);
            std::string fwkey;
            std::string revkey;
            bool fwfound = true;
            bool revfound = true;
            if (m_barcodes_fw == nullptr)
                fwkey = dummykey;
            else {
                fwfound = false;
                for (const auto &code : m_barcodes_fw->fw) {
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
                for (const auto &code : m_barcodes_rev->rev) {
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

            if (fwfound && revfound)
               ++localcounts[std::make_pair(fwkey, revkey)][rd.getSequence().substr(insertmatch.first, insertmatch.second)];
            //TODO: write to mismatches file
        } catch (std::exception &e) {
            //TODO: write to mismatches file
        }
    }
    sync->countsmutex.lock();
    for (const auto &codes : localcounts) {
        for (const auto &insert : codes.second) {
            m_counts[codes.first][insert.first] += insert.second;
        }
    }
    sync->countsmutex.unlock();
}

const std::unordered_map<std::pair<std::string, std::string>, std::unordered_map<std::string, uint64_t>>& ReadCounter::getCounts() const
{
    return m_counts;
}
