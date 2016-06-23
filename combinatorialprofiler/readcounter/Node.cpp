#include "Node.h"
#include "util.h"

#include <stdexcept>

static std::string dummykey = std::string();

NodeBase::NodeBase() : experiment(nullptr) {}

NodeBase::NodeBase(std::string seq)
: experiment(nullptr), sequence(std::move(seq)) {}

NodeBase::~NodeBase() {}

BarcodeNode::BarcodeNode() : Node(dummykey, 0) {}

BarcodeNode::~BarcodeNode() {}

BarcodeMatch* BarcodeNode::match(Read &rd) const
{
    return new BarcodeMatch(this, SIZE_MAX, 0);
}

BarcodeNode::BarcodeNode(std::string seq, std::vector<std::string> uniqueseq, float mismatches)
: Node(std::move(seq), mismatches), m_uniqueSequences(std::move(uniqueseq))
{}

RevBarcodeNode::RevBarcodeNode() : BarcodeNode() {}

RevBarcodeNode::RevBarcodeNode(std::string fullseq, std::vector<std::string> uniqueSeqs, float mismatches)
: BarcodeNode(std::move(fullseq), std::move(uniqueSeqs), mismatches)
{
    for (auto &s : m_uniqueSequences) {
        s = revCompl(s);
    }
}

BarcodeMatch* RevBarcodeNode::match(Read &rd) const
{
    auto rdseq = rd.getSequence();
    std::pair<decltype(m_uniqueSequences)::size_type, std::string::size_type> bestFind(0, SIZE_MAX);
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

FwBarcodeNode::FwBarcodeNode() : BarcodeNode() {}

FwBarcodeNode::FwBarcodeNode(std::string fullseq, std::vector<std::string> uniqueSeqs, float mismatches)
: BarcodeNode(std::move(fullseq), std::move(uniqueSeqs), mismatches)
{}

BarcodeMatch* FwBarcodeNode::match(Read &rd) const
{
    auto rdseq = rd.getSequence();
    std::pair<decltype(m_uniqueSequences)::size_type, std::string::size_type> bestFind(0, SIZE_MAX);
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
}

InsertMatch* InsertNode::match(Read &read) const
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
