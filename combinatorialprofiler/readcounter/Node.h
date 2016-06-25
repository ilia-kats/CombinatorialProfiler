#ifndef NODE_H
#define NODE_H

#include "Read.h"
#include "Match.h"
#include "util.h"

#include <string>
#include <vector>
#include <utility>
#include <array>

class Experiment;
class MatchBase;
class InsertMatch;

extern std::string node_dummykey;

class NodeBase
{
public:
    NodeBase();
    NodeBase(std::string);

    virtual ~NodeBase();

    virtual MatchBase* match(Read&) const = 0;

    Experiment *experiment;

    std::vector<NodeBase*> children;
    std::string sequence;

protected:
    template<class NeedleIt, class HaystackIt>
    static std::pair<std::string::size_type, std::string::size_type> fuzzy_find(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin, HaystackIt hend);
};

template<typename T>
class Node : public NodeBase
{
public:
    typedef T mismatch_type;

    Node(std::string seq, T mismatches)
    : NodeBase(seq), m_allowedMismatches(mismatches) {}
    virtual ~Node() {}

    T allowedMismatches() const
    {
        return m_allowedMismatches;
    }

protected:
    T m_allowedMismatches;
};


class FwBarcodeNode
{
public:
    virtual std::array<std::string::const_iterator, 4> getIterators(const std::string&, const Read&, std::string::size_type) const;
};

class RevBarcodeNode
{
public:
    virtual std::array<std::string::const_reverse_iterator, 4> getIterators(const std::string&, const Read&, std::string::size_type) const;
};

template<class T>
class HammingBarcodeNode : public Node<float>
{
public:
    typedef HammingBarcodeMatch<T> match_type;

    HammingBarcodeNode(std::string fullseq, float mismatches, std::vector<std::string> uniqueSeqs);

    virtual match_type* match(Read &rd) const;

protected:
    HammingBarcodeNode();
    std::vector<std::string> m_uniqueSequences;
};

class FwHammingBarcodeNode : public HammingBarcodeNode<FwHammingBarcodeNode>, public FwBarcodeNode
{
public:
    FwHammingBarcodeNode(std::string, float mismatches, std::vector<std::string>);
};

class RevHammingBarcodeNode: public HammingBarcodeNode<RevHammingBarcodeNode>, public RevBarcodeNode
{
public:
    RevHammingBarcodeNode(std::string, float mismatches, std::vector<std::string>);
};

class DummyHammingBarcodeNode : public HammingBarcodeNode<DummyHammingBarcodeNode>, public FwBarcodeNode
{
public:
    DummyHammingBarcodeNode();
    virtual HammingBarcodeNode<DummyHammingBarcodeNode>::match_type* match(Read&) const;
};

template<class T>
class SeqlevBarcodeNode : public Node<std::string::size_type>
{
public:
    typedef SeqlevBarcodeMatch<T> match_type;

    SeqlevBarcodeNode(std::string fullseq, std::string::size_type mismatches);

    virtual match_type* match(Read &rd) const;

protected:
    SeqlevBarcodeNode();
    std::string m_sequence;
};

class FwSeqlevBarcodeNode : public SeqlevBarcodeNode<FwSeqlevBarcodeNode>, public FwBarcodeNode
{
public:
    FwSeqlevBarcodeNode(std::string, std::string::size_type);
};

class RevSeqlevBarcodeNode : public SeqlevBarcodeNode<RevSeqlevBarcodeNode>, public RevBarcodeNode
{
public:
    RevSeqlevBarcodeNode(std::string, std::string::size_type);
};

class DummySeqlevBarcodeNode : public SeqlevBarcodeNode<DummySeqlevBarcodeNode>, public FwBarcodeNode
{
public:
    DummySeqlevBarcodeNode();
    virtual SeqlevBarcodeNode<DummySeqlevBarcodeNode>::match_type* match(Read&) const;
};

class InsertNode : public Node<std::string::size_type>
{
public:
    typedef InsertMatch match_type;

    InsertNode(std::string, std::string::size_type mismatches);
    virtual match_type* match(Read&) const = 0;

protected:
    std::string m_upstreamseq;
    std::string m_downstreamseq;
    std::string::size_type m_insertlength;
    std::string::size_type m_minReadLength;

    match_type* getMatch(const Read&, std::string::size_type, std::string::size_type, std::string::size_type) const;
};

class ExactMatchingInsertNode : public InsertNode
{
public:
    ExactMatchingInsertNode(std::string);
    virtual InsertNode::match_type* match(Read&) const;
};

class HammingMatchingInsertNode : public InsertNode
{
public:
    using InsertNode::InsertNode;
    virtual InsertNode::match_type* match(Read&) const;
};

template<class T>
HammingBarcodeNode<T>::HammingBarcodeNode(std::string fullseq, float mismatches, std::vector<std::string> uniqueSeqs)
: Node(std::move(fullseq), mismatches), m_uniqueSequences(std::move(uniqueSeqs))
{}

template<class T>
typename HammingBarcodeNode<T>::match_type* HammingBarcodeNode<T>::match(Read &rd) const
{
    const T* cthis = static_cast<const T*>(this);
    std::pair<typename decltype(m_uniqueSequences)::size_type, std::string::size_type> bestFind(0, SIZE_MAX);
    if (rd.length() > sequence.size()) {
        if (!m_allowedMismatches){
            for (const auto &seq : m_uniqueSequences) {
                auto tomatch = cthis->getIterators(seq, rd, seq.size());
                if (std::equal(std::get<0>(tomatch), std::get<1>(tomatch), std::get<2>(tomatch)))
                    return new HammingBarcodeMatch<T>(this, seq.size(), 0);
            }
        } else {
            auto it = m_uniqueSequences.cbegin();
            for (typename decltype(m_uniqueSequences)::size_type i = 0; i < m_uniqueSequences.size(); ++i, ++it) {
                auto tomatch = cthis->getIterators(*it, rd, it->size());
                auto found = hamming_distance(tomatch[0], tomatch[1], tomatch[2]);
                if (found < bestFind.second) {
                    bestFind.first = i;
                    bestFind.second = found;
                }
                if (!found)
                    break;
            }

        }
    }
    return new HammingBarcodeMatch<T>(this, m_uniqueSequences[bestFind.first].size(), bestFind.second);
}

template<class T>
HammingBarcodeNode<T>::HammingBarcodeNode()
: Node(node_dummykey, 0)
{}

template<class T>
SeqlevBarcodeNode<T>::SeqlevBarcodeNode(std::string fullseq, std::string::size_type mismatches)
: Node(std::move(fullseq), mismatches), m_sequence(sequence)
{}

template<class T>
SeqlevBarcodeNode<T>::SeqlevBarcodeNode()
: Node(node_dummykey, 0)
{}

template<class T>
typename SeqlevBarcodeNode<T>::match_type* SeqlevBarcodeNode<T>::match(Read &rd) const
{
    if (rd.length() < m_sequence.size())
        return new SeqlevBarcodeMatch<T>(this, SIZE_MAX);
    auto length = std::min(rd.length(), m_sequence.size() + allowedMismatches()); // account for insertions
    auto tomatch = static_cast<const T*>(this)->getIterators(m_sequence, rd, length);
    return new SeqlevBarcodeMatch<T>(this, seqlev_distance(std::get<0>(tomatch), std::get<1>(tomatch), std::get<2>(tomatch), std::get<3>(tomatch), m_sequence.size(), length));
}

#endif
