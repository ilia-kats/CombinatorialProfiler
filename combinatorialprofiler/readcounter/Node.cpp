#include "Node.h"
#include "util.h"

#include <stdexcept>
#include <functional>

std::string node_dummykey = std::string();

NodeBase::NodeBase() : experiment(nullptr) {}

NodeBase::NodeBase(std::string seq)
: experiment(nullptr), sequence(std::move(seq)) {}

NodeBase::~NodeBase() {}

template<class NeedleIt, class HaystackIt>
std::pair<std::string::size_type, std::string::size_type> NodeBase::fuzzy_find(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin, HaystackIt hend)
{
    auto hsize = std::distance(hbegin, hend);
    auto nsize = std::distance(nbegin, nend);
    auto totest = hsize - nsize + 1;
    if (totest <= 0)
        return std::make_pair(0, UINT16_MAX);
    std::pair<std::string::size_type, std::string::size_type> bestMatch(0, SIZE_MAX);
    for (std::string::size_type i = 0; i < totest; ++i) {
        auto mm = hamming_distance(nbegin, nend, hbegin + i);
        if (mm < bestMatch.second) {
            bestMatch.first = i;
            bestMatch.second = mm;
        }
        if (!mm)
            break;
    }
    return bestMatch;
}

std::array<std::string::const_iterator, 4> FwBarcodeNode::getIterators(const std::string &seq, const Read &rd, std::string::size_type l) const
{
    return {seq.cbegin(), seq.cend(), rd.getSequence().cbegin(), rd.getSequence().cbegin() + l};
}

std::array<std::string::const_reverse_iterator, 4> RevBarcodeNode::getIterators(const std::string &seq, const Read &rd, std::string::size_type l) const
{
    return {seq.crbegin(), seq.crend(), rd.getSequence().crbegin(), rd.getSequence().crbegin() + l};
}

FwHammingBarcodeNode::FwHammingBarcodeNode(std::string fullseq, float mismatches, std::vector<std::string> uniqueSeqs)
: HammingBarcodeNode(std::move(fullseq), mismatches, std::move(uniqueSeqs))
{}

RevHammingBarcodeNode::RevHammingBarcodeNode(std::string fullseq, float mismatches, std::vector<std::string> uniqueSeqs)
: HammingBarcodeNode(std::move(fullseq), mismatches, std::move(uniqueSeqs))
{
    for (auto &s : m_uniqueSequences) {
        s = revCompl(s);
    }
}

DummyHammingBarcodeNode::DummyHammingBarcodeNode()
{}

HammingBarcodeNode<DummyHammingBarcodeNode>::match_type* DummyHammingBarcodeNode::match(Read &rd) const
{
    return new HammingBarcodeMatch<DummyHammingBarcodeNode>(this);
}

FwSeqlevBarcodeNode::FwSeqlevBarcodeNode(std::string fullseq, std::string::size_type mismatches)
: SeqlevBarcodeNode(std::move(fullseq), mismatches)
{}

RevSeqlevBarcodeNode::RevSeqlevBarcodeNode(std::string fullseq, std::string::size_type mismatches)
: SeqlevBarcodeNode(std::move(fullseq), mismatches)
{
    m_sequence = revCompl(m_sequence);
}

DummySeqlevBarcodeNode::DummySeqlevBarcodeNode()
{}

SeqlevBarcodeNode<DummySeqlevBarcodeNode>::match_type* DummySeqlevBarcodeNode::match(Read &rd) const
{
    return new SeqlevBarcodeMatch<DummySeqlevBarcodeNode>(this);
}

InsertNode::InsertNode(std::string seq, std::string::size_type mismatches)
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
    m_minReadLength = m_insertlength + m_upstreamseq.size() + m_downstreamseq.size();
}

InsertNode::match_type* InsertNode::getMatch(const Read &rd, std::string::size_type start, std::string::size_type end, std::string::size_type mismatches) const
{
    std::string insert;
    auto it = rd.getSequence().cbegin();
    std::copy(it + start, it + end, std::inserter(insert, insert.begin()));
    if (insert.find('N') != std::string::npos || insert.find('n') != std::string::npos)
        return new InsertMatch(this, false);
    return new InsertMatch(this, mismatches, std::move(insert));
}

ExactMatchingInsertNode::ExactMatchingInsertNode(std::string seq)
: InsertNode(std::move(seq), 0)
{}

InsertNode::match_type* ExactMatchingInsertNode::match(Read &read) const
{
    if (read.length() < m_minReadLength)
        return new InsertMatch(this, false);
    auto firstmatch = std::ref(m_upstreamseq);
    bool matchDownstream = true;
    auto toadd = firstmatch.get().size();
    if (firstmatch.get().empty()) {
        firstmatch = std::ref(m_downstreamseq);
        matchDownstream = false;
        toadd = -m_insertlength;
    }
    auto start = read.getSequence().find(firstmatch.get());
    if (start == decltype(firstmatch)::type::npos) {
        read = read.reverseComplement();
        start = read.getSequence().find(firstmatch);
        if (start == decltype(firstmatch)::type::npos)
            return new InsertMatch(this, false);
    }
    start += toadd;
    decltype(start) end;
    if (matchDownstream && !m_downstreamseq.empty()) {
        end = read.getSequence().find(m_downstreamseq, start);
        if (end == decltype(m_downstreamseq)::npos || end - start != m_insertlength)
            return new InsertMatch(this, false);
    } else {
        end = start + m_insertlength;
    }
    return getMatch(read, start, end, 0);
}

InsertNode::match_type* HammingMatchingInsertNode::match(Read &read) const
{
    if (read.length() < m_minReadLength)
        return new InsertMatch(this, false);
    decltype(m_upstreamseq)::size_type start, end;
    std::string::size_type mismatches_cnt = 0;

    auto firstmatch = std::ref(m_upstreamseq);
    bool matchDownstream = true;
    auto toadd = firstmatch.get().size();
    if (firstmatch.get().empty()) {
        firstmatch = std::ref(m_downstreamseq);
        matchDownstream = false;
        toadd = -m_insertlength;
    }

    auto upstream = fuzzy_find(firstmatch.get().cbegin(), firstmatch.get().cend(), read.getSequence().cbegin(), read.getSequence().cend());
    if (upstream.second > m_allowedMismatches) {
        read = read.reverseComplement();
        upstream = fuzzy_find(firstmatch.get().cbegin(), firstmatch.get().cend(), read.getSequence().cbegin(), read.getSequence().cend());
        if (upstream.second > m_allowedMismatches)
            return new InsertMatch(this, upstream.second);
    }
    start = upstream.first + toadd;
    if (matchDownstream && !m_downstreamseq.empty()) {
        auto downstream = fuzzy_find(m_downstreamseq.cbegin(), m_downstreamseq.cend(), read.getSequence().cbegin() + start, read.getSequence().cend());
        downstream.first += start;
        mismatches_cnt = downstream.second + upstream.second;
        if (mismatches_cnt > m_allowedMismatches)
            return new InsertMatch(this, mismatches_cnt);
        end = downstream.first;
        if (end - start != m_insertlength)
        return new InsertMatch(this, false);
    } else {
        end = start + m_insertlength;
    }

    return getMatch(read, start, end, mismatches_cnt);
}

