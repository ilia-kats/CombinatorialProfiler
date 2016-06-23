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


class BarcodeNode : public Node<float>
{
public:
    BarcodeNode();
    virtual ~BarcodeNode();

    virtual BarcodeMatch* match(Read&) const;

protected:
    BarcodeNode(std::string, std::vector<std::string>, float mismatches=0);

    std::vector<std::string> m_uniqueSequences;
};

class MatchingBarcodeNode : public BarcodeNode
{
public:
    virtual BarcodeMatch *match(Read&) const;

protected:
    MatchingBarcodeNode(std::string, std::vector<std::string>, float mismatches=0);

private:
    virtual std::string::const_iterator getReadPart(const Read&, std::string::size_type) const = 0;
};

class RevBarcodeNode : public MatchingBarcodeNode
{
public:
    RevBarcodeNode(std::string, std::vector<std::string>, float mismatches=0);

private:
    virtual std::string::const_iterator getReadPart(const Read&, std::string::size_type) const;
};

class FwBarcodeNode : public MatchingBarcodeNode
{
public:
    FwBarcodeNode(std::string, std::vector<std::string>, float mismatches=0);

private:
    virtual std::string::const_iterator getReadPart(const Read&, std::string::size_type) const;
};

class InsertNode : public Node<std::string::size_type>
{
public:
    InsertNode(std::string, std::string::size_type mismatches = 1);
    virtual InsertMatch* match(Read&) const;

private:
    std::string m_upstreamseq;
    std::string m_downstreamseq;
    std::string::size_type m_insertlength;
};


#endif
