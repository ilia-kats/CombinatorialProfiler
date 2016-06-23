#ifndef MATCH_H
#define MATCH_H

#include <string>

class NodeBase;
class BarcodeNode;
class InsertNode;

class MatchBase
{
public:
    MatchBase(const NodeBase*);
    MatchBase(const NodeBase*, bool);
    virtual ~MatchBase();

    operator bool() const;

    const NodeBase* node() const;

protected:
    const NodeBase *m_node;
    bool m_match;
};

template<typename T>
class Match : public MatchBase
{
public:
    virtual ~Match();

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
    Match(const NodeBase *n, bool m);

    Match(const NodeBase *n, bool m, T mismatches);

    T m_mismatches;
};

class BarcodeMatch : public Match<float>
{
public:
    BarcodeMatch(const BarcodeNode*, std::string::size_type, std::string::size_type);

private:
    std::string::size_type m_matchedLength;
    std::string::size_type m_actualMismatches;
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

#endif
