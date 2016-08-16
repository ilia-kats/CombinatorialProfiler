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

bool operator<(const MatchBase &l, const MatchBase &r)
{
    if (l and !r)
        return true;
    else if (r and !l)
        return false;
    else
        return l.lessThan(r);
}

bool operator>(const MatchBase &l, const MatchBase &r)
{
    return r < l;
}

bool operator<=(const MatchBase &l, const MatchBase &r)
{
    return !(l > r);
}

bool operator>=(const MatchBase &l, const MatchBase &r)
{
    return !(l < r);
}

bool operator==(const MatchBase &l, const MatchBase &r)
{
    return !(l < r || l > r);
}

bool operator!=(const MatchBase &l, const MatchBase &r)
{
    return !(l == r);
}

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
