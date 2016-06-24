#ifndef MATCH_H
#define MATCH_H

#include <string>

class NodeBase;
class HammingBarcodeNode;
class SeqlevBarcodeNode;
class InsertNode;

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
    virtual ~Match();

    virtual bool perfectMatch() const;

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

class HammingBarcodeMatch : public Match<float>
{
public:
    HammingBarcodeMatch(const HammingBarcodeNode*, std::string::size_type, std::string::size_type);
    HammingBarcodeMatch(const HammingBarcodeNode*);

private:
    std::string::size_type m_matchedLength;
    std::string::size_type m_actualMismatches;
};

class SeqlevBarcodeMatch : public Match<std::string::size_type>
{
public:
    SeqlevBarcodeMatch(const SeqlevBarcodeNode*, std::string::size_type);
    SeqlevBarcodeMatch(const SeqlevBarcodeNode*);
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
