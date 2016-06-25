#ifndef MATCH_H
#define MATCH_H

#include <string>

class NodeBase;
class InsertNode;

template<class T>
class HammingBarcodeNode;

template<class T>
class SeqlevBarcodeNode;

class MatchBase
{
public:
    MatchBase(const NodeBase*);
    MatchBase(const NodeBase*, bool);
    virtual ~MatchBase();

    operator bool() const;

    const NodeBase* node() const;
    virtual bool perfectMatch() const = 0;

protected:
    const NodeBase *m_node;
    bool m_match;
};

template<typename T>
class Match : public MatchBase
{
public:
    virtual ~Match() {}

    virtual bool perfectMatch() const
    {
        return !m_mismatches;
    }

    template<typename S>
    friend bool operator<(const Match<S> &l, const Match<S> &r);
    template<typename S>
    friend bool operator>(const Match<S> &l, const Match<S> &r);
    template<typename S>
    friend bool operator<=(const Match<S> &l, const Match<S> &r);
    template<typename S>
    friend bool operator>=(const Match<S> &l, const Match<S> &r);
    template<typename S>
    friend bool operator==(const Match<S> &l, const Match<S> &r);
    template<typename S>
    friend bool operator!=(const Match<S> &l, const Match<S> &r);

protected:
    Match(const NodeBase *n, bool m)
    : MatchBase(n, m), m_mismatches(0) {}

    Match(const NodeBase *n, bool m, T mismatches)
    : MatchBase(n, m), m_mismatches(mismatches)
    {}

    T m_mismatches;
};

template<class T>
class HammingBarcodeMatch : public Match<float>
{
public:
    HammingBarcodeMatch(const HammingBarcodeNode<T>*, std::string::size_type, std::string::size_type);
    HammingBarcodeMatch(const HammingBarcodeNode<T>*);
private:
    std::string::size_type m_matchedLength;
    std::string::size_type m_actualMismatches;
};

template<class T>
class SeqlevBarcodeMatch : public Match<std::string::size_type>
{
public:
    SeqlevBarcodeMatch(const SeqlevBarcodeNode<T>*, std::string::size_type);

    SeqlevBarcodeMatch(const SeqlevBarcodeNode<T>*);
};

class InsertMatch : public Match<std::string::size_type>
{
public:
    InsertMatch(const InsertNode*, std::string::size_type, std::string insert = std::string());

    InsertMatch(const InsertNode *n, bool match);

    const std::string& insertSequence() const;

private:
    std::string m_insert;
};

#include "Node.h"
template<class T>
HammingBarcodeMatch<T>::HammingBarcodeMatch(const HammingBarcodeNode<T> *n, std::string::size_type length, std::string::size_type mismatches)
: Match(n, false, (float)mismatches / (float)length), m_matchedLength(length), m_actualMismatches(mismatches)
{
    if (m_mismatches > n->allowedMismatches())
        m_match = false;
    else
        m_match = true;
}

template<class T>
HammingBarcodeMatch<T>::HammingBarcodeMatch(const HammingBarcodeNode<T> *n)
: Match(n, true)
{}

template<class T>
SeqlevBarcodeMatch<T>::SeqlevBarcodeMatch(const SeqlevBarcodeNode<T> *n, std::string::size_type mismatches)
: Match(n, mismatches <= n->allowedMismatches())
{}

template<class T>
SeqlevBarcodeMatch<T>::SeqlevBarcodeMatch(const SeqlevBarcodeNode<T> *n)
: Match(n, true)
{}

#endif
