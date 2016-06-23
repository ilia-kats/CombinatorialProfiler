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
    static std::pair<std::string::size_type, std::string::size_type> fuzzy_find(NeedleIt nbegin, NeedleIt nend, HaystackIt hbegin, HaystackIt hend)
    {
        auto hsize = std::distance(hbegin, hend);
        auto nsize = std::distance(nbegin, nend);
        auto totest = hsize - nsize + 1;
        if (totest <= 0)
            return std::make_pair(0, UINT16_MAX);
        std::pair<std::string::size_type, std::string::size_type> bestMatch(0, SIZE_MAX);
        for (std::string::size_type i = 0; i < totest; ++i) {
            auto mm = std::inner_product(nbegin,
                                         nend,
                                         hbegin + i,
                                         static_cast<std::string::size_type>(0),
                                         std::plus<std::string::size_type>(),
                                         [](const typename decltype(nbegin)::value_type &n, const typename decltype(hbegin)::value_type &h) -> bool {return static_cast<bool>(n ^ h);});
            if (mm < bestMatch.second) {
                bestMatch.first = i;
                bestMatch.second = mm;
            }
            if (!mm)
                break;
        }
        return bestMatch;
    }
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

class RevBarcodeNode : public BarcodeNode
{
public:
    RevBarcodeNode();
    RevBarcodeNode(std::string, std::vector<std::string>, float mismatches=0);

    virtual BarcodeMatch* match(Read&) const;
};

class FwBarcodeNode : public BarcodeNode
{
public:
    FwBarcodeNode();
    FwBarcodeNode(std::string, std::vector<std::string>, float mismatches=0);

    virtual BarcodeMatch* match(Read&) const;
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
