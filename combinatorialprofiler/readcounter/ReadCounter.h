#ifndef READCOUNTER_H
#define READCOUNTER_H

#include "Experiment.h"
#include "Read.h"
#include "Node.h"
#include "Match.h"
#include "util.h"

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <string>
#include <vector>
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

class InsertNode;
class NodeBase;

template<class BarcodeNodeBase>
class ReadCounter
{
public:
    virtual ~ReadCounter()
    {
        for (auto &n : m_nodes)
            delete n;
    }

    void countReads(const std::string &file, const std::string &outprefix, int threads=1)
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

    uint16_t allowedMismatches() const
    {
        return m_allowedMismatches;
    }

    uint64_t read() const
    {
        return m_read;
    }

    uint64_t counted() const
    {
        return m_counted;
    }

    uint64_t unmatchedTotal() const
    {
        return m_unmatchedTotal;
    }

    uint64_t unmatchedInsert() const
    {
        return m_unmatchedInsert;
    }

    uint64_t unmatchedBarcodes() const
    {
        return m_unmatchedBarcodes;
    }

    uint64_t unmatchedInsertSequence() const
    {
        return m_unmatchedInsertSequence;
    }

    uint64_t written()const
    {
        return m_written;
    }

protected:
    ReadCounter(std::vector<Experiment*> experiments, uint16_t insert_mismatches = 1)
    : m_allowedMismatches(insert_mismatches), m_read(0), m_counted(0), m_unmatchedTotal(0), m_unmatchedInsert(0), m_unmatchedBarcodes(0), m_unmatchedInsertSequence(0), m_written(0), m_experiments(std::move(experiments)) {}

    void init()
    {
        std::unordered_map<std::string, InsertNode*> inserts;
        for (const auto &exp : m_experiments) {
            if (!inserts.count(exp->insert)) {
                InsertNode *n = new InsertNode(exp->insert, m_allowedMismatches);
                inserts[exp->insert] = n;
                m_nodes.push_back(n);
                m_tree.push_back(n);
            }
        }

        std::unordered_map<std::string, std::unordered_map<std::string, BarcodeNodeBase*>> fwNodes;
        std::unordered_map<std::string,BarcodeNodeBase*> dummyNodes;

        for (const auto &exp : m_experiments) {
            std::vector<BarcodeNodeBase*> revNodes;
            for (const auto &c : exp->revBarcodeSet) {
                BarcodeNodeBase *n = makeRevCodeNode(exp, c.first);
                n->experiment = exp;
                revNodes.push_back(n);
                m_nodes.push_back(n);
            }

            for (const auto &c : exp->fwBarcodeSet) {
                if (!fwNodes[exp->insert].count(c.first)) {
                    BarcodeNodeBase *n = makeFwCodeNode(exp, c.first);
                    fwNodes[exp->insert][c.first] = n;
                    m_nodes.push_back(n);
                    inserts[exp->insert]->children.push_back(n);
                }
                std::copy(revNodes.cbegin(), revNodes.cend(), std::inserter(fwNodes[exp->insert][c.first]->children, fwNodes[exp->insert][c.first]->children.end()));
            }
            if (!exp->fwBarcodeSet.size()) {
                if (!dummyNodes.count(exp->insert)) {
                    BarcodeNodeBase *n = makeDummyCodeNode();
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
                BarcodeNodeBase *n = makeDummyCodeNode();
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

    const std::vector<Experiment*>& experiments() const
    {
        return m_experiments;
    }

private:
    struct ThreadSynchronization
    {
        struct FailedMatch
        {
            enum class Fail {noFail, insertFailed, barcodeFailed, insertMatchFailed};
            Read read;
            Fail fail;
            std::vector<MatchBase*> nodes;

            FailedMatch(Read r, Fail f, std::vector<MatchBase*> ns)
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

    uint16_t m_allowedMismatches;

    uint64_t m_read;
    uint64_t m_counted;
    uint64_t m_unmatchedTotal;
    uint64_t m_unmatchedInsert;
    uint64_t m_unmatchedBarcodes;
    uint64_t m_unmatchedInsertSequence;
    uint64_t m_written;

    std::vector<Experiment*> m_experiments;

    std::vector<InsertNode*> m_tree;
    std::vector<NodeBase*> m_nodes;

    void readFile(const std::string &file, ThreadSynchronization *sync)
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

    void matchRead(ThreadSynchronization *sync)
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

            typename ThreadSynchronization::FailedMatch::Fail fail = ThreadSynchronization::FailedMatch::Fail::noFail;
            std::vector<MatchBase*> nodes;

            std::string ins;
            InsertMatch *in = nullptr;
            for (const auto &n : m_tree) {
                InsertMatch *cn = n->match(rd);
                if (*cn) {
                    nodes.push_back(cn);
                    in = cn;
                    break;
                }
                else
                    delete cn;
            }
            if (in && *in) {
                MatchBase *cn = in;
                while (!cn->node()->experiment && cn->node()->children.size()) {
                    MatchBase *best = nullptr;
                    for (const auto &n : cn->node()->children) {
                        MatchBase *bn = n->match(rd);
                        if (*bn && best && *bn < *best) {
                            delete best;
                            best = bn;
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

    void writeReads(const std::string &prefix, ThreadSynchronization *sync)
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

    virtual BarcodeNodeBase* makeFwCodeNode(const Experiment*, const std::string&) = 0;
    virtual BarcodeNodeBase* makeRevCodeNode(const Experiment*, const std::string&) = 0;
    virtual BarcodeNodeBase* makeDummyCodeNode() = 0;
};

#endif
