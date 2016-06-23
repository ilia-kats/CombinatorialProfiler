#include "Node.h"
#include "util.h"

#include <stdexcept>

static std::string dummykey = std::string();

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

template<class NeedleIt, class HaystackIt>
std::string::size_type NodeBase::hamming_distance(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin)
{
    return std::inner_product(nbegin,
                              nend,
                              hbegin,
                              static_cast<std::string::size_type>(0),
                              std::plus<std::string::size_type>(),
                               [](const typename decltype(nbegin)::value_type &n, const typename decltype(hbegin)::value_type &h) -> bool {return static_cast<bool>(n ^ h);});
}

std::string::const_iterator FwBarcodeNode::getReadPart(const Read &rd, std::string::size_type l) const
{
    return rd.getSequence().cbegin();
}

std::string::const_iterator RevBarcodeNode::getReadPart(const Read &rd, std::string::size_type l) const
{
    return rd.getSequence().cend() - l;
}

HammingBarcodeNode::HammingBarcodeNode(std::string fullseq, float mismatches, std::vector<std::string> uniqueSeqs)
: Node(std::move(fullseq), mismatches), m_uniqueSequences(std::move(uniqueSeqs))
{}

HammingBarcodeNode::HammingBarcodeNode()
: Node(dummykey, 0)
{}

HammingBarcodeMatch* HammingBarcodeNode::match(Read &rd) const
{
    std::pair<decltype(m_uniqueSequences)::size_type, std::string::size_type> bestFind(0, SIZE_MAX);
    if (!m_allowedMismatches){
        for (const auto &seq : m_uniqueSequences) {
            auto tomatch = getReadPart(rd, seq.size());
            if (std::equal(seq.cbegin(), seq.cend(), tomatch))
                return new HammingBarcodeMatch(this, seq.size(), 0);
        }
    } else {
        auto it = m_uniqueSequences.cbegin();
        for (decltype(m_uniqueSequences)::size_type i = 0; i < m_uniqueSequences.size(); ++i, ++it) {
            auto tomatch = getReadPart(rd, it->size());
            auto found = hamming_distance(it->cbegin(), it->cend(), tomatch);
            if (found < bestFind.second) {
                bestFind.first = i;
                bestFind.second = found;
            }
            if (!found)
                break;
        }

    }
    return new HammingBarcodeMatch(this, m_uniqueSequences[bestFind.first].size(), bestFind.second);
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

HammingBarcodeMatch* DummyHammingBarcodeNode::match(Read &rd) const
{
    return new HammingBarcodeMatch(this, SIZE_MAX, 0);
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
