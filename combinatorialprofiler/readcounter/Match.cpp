#include "Match.h"
#include "Node.h"

MatchBase::MatchBase(const NodeBase *n)
: m_node(n), m_match(true) {}

MatchBase::MatchBase(const NodeBase *n, bool m)
: m_node(n), m_match(m) {}

MatchBase::~MatchBase() {}

MatchBase::operator bool() const
{
    return m_match;
}

const NodeBase* MatchBase::node() const
{
    return m_node;
}

template<typename T>
bool operator<(const Match<T> &l, const Match<T> &r)
{
    if (l and !r)
        return true;
    else if (r and !l)
        return false;
    else
        return l.m_mismatches < r.m_mismatches;
}

template<typename T>
bool operator>(const Match<T> &l, const Match<T> &r)
{
    return r < l;
}

template<typename T>
bool operator<=(const Match<T> &l, const Match<T> &r)
{
    return !(l > r);
}

template<typename T>
bool operator>=(const Match<T> &l, const Match<T> &r)
{
    return !(l < r);
}

template<typename T>
bool operator==(const Match<T> &l, const Match<T> &r)
{
    return l.match && r.match && l.m_mismatches == r.m_mismatches;
}

template<typename T>
bool operator!=(const Match<T> &l, const Match<T> &r)
{
    return !(l == r);
}

template<typename T>
Match<T>::Match(const NodeBase *n, bool m)
: MatchBase(n, m), m_mismatches(0) {}

template<typename T>
Match<T>::Match(const NodeBase *n, bool m, T mismatches)
: MatchBase(n, m), m_mismatches(mismatches)
{}

template<typename T>
Match<T>::~Match() {}

template<typename T>
bool Match<T>::perfectMatch() const
{
    return !m_mismatches;
}

HammingBarcodeMatch::HammingBarcodeMatch(const HammingBarcodeNode *n, std::string::size_type length, std::string::size_type mismatches)
: Match(n, false, (float)mismatches / (float)length), m_matchedLength(length), m_actualMismatches(mismatches)
{
    if (m_mismatches > n->allowedMismatches())
        m_match = false;
    else
        m_match = true;
}

HammingBarcodeMatch::HammingBarcodeMatch(const HammingBarcodeNode *n)
: Match(n, true)
{}

SeqlevBarcodeMatch::SeqlevBarcodeMatch(const SeqlevBarcodeNode *n, std::string::size_type mismatches)
: Match(n, mismatches <= n->allowedMismatches())
{}

SeqlevBarcodeMatch::SeqlevBarcodeMatch(const SeqlevBarcodeNode *n)
: Match(n, true)
{}

InsertMatch::InsertMatch(const InsertNode *n, std::string::size_type mismatches, std::string insert)
: Match(n, mismatches <= n->allowedMismatches(), mismatches), m_insert(std::move(insert))
{}

InsertMatch::InsertMatch(const InsertNode *n, bool match)
: Match(n, match)
{}

const std::string& InsertMatch::insertSequence() const
{
    return m_insert;
}
