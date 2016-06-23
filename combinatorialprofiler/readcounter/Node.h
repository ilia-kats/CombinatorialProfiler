#ifndef NODE_H
#define NODE_H

#include "Read.h"
#include "Match.h"

#include <string>
#include <vector>
#include <utility>
#include <numeric>

class Experiment;

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

    template<class NeedleIt, class HaystackIt>
    static std::string::size_type hamming_distance(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin);
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

class BarcodeNode
{
protected:
    virtual std::string::const_iterator getReadPart(const Read&, std::string::size_type) const = 0;
};

class FwBarcodeNode : virtual public BarcodeNode
{
protected:
    virtual std::string::const_iterator getReadPart(const Read&, std::string::size_type) const;
};

class RevBarcodeNode : virtual public BarcodeNode
{
protected:
    virtual std::string::const_iterator getReadPart(const Read&, std::string::size_type) const;
};

class HammingBarcodeNode : public Node<float>, virtual public BarcodeNode
{
public:
    typedef HammingBarcodeMatch match_type;

    HammingBarcodeNode(std::string, float mismatches, std::vector<std::string>);

    virtual HammingBarcodeMatch* match(Read&) const;

protected:
    HammingBarcodeNode();
    std::vector<std::string> m_uniqueSequences;
};

class FwHammingBarcodeNode : virtual public HammingBarcodeNode, virtual public FwBarcodeNode
{
public:
    FwHammingBarcodeNode(std::string, float mismatches, std::vector<std::string>);
};

class RevHammingBarcodeNode: virtual public HammingBarcodeNode, virtual public RevBarcodeNode
{
public:
    RevHammingBarcodeNode(std::string, float mismatches, std::vector<std::string>);
};

class DummyHammingBarcodeNode : virtual public HammingBarcodeNode, virtual public FwBarcodeNode
{
public:
    DummyHammingBarcodeNode();
    virtual HammingBarcodeMatch *match(Read&) const;
};

class InsertNode : public Node<std::string::size_type>
{
public:
    typedef InsertMatch match_type;

    InsertNode(std::string, std::string::size_type mismatches);
    virtual InsertMatch* match(Read&) const;

private:
    std::string m_upstreamseq;
    std::string m_downstreamseq;
    std::string::size_type m_insertlength;
};


#endif
