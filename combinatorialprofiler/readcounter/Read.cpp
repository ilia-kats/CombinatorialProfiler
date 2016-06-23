#include "Read.h"
#include "util.h"

Read::Read() {}

Read::Read(std::string name)
: m_name(std::move(name))
{}

Read::Read(std::string name, std::string sequence)
: m_name(std::move(name)), m_sequence(std::move(sequence))
{}

Read::Read(std::string name, std::string sequence, std::string description)
: m_name(std::move(name)), m_sequence(std::move(sequence)), m_description(std::move(description))
{}

Read::Read(std::string name, std::string sequence, std::string description, std::string quality)
: m_name(std::move(name)), m_sequence(std::move(sequence)), m_description(std::move(description)), m_quality(std::move(quality))
{}

void Read::setName(std::string name)
{
    m_name = std::move(name);
}
const std::string& Read::getName() const
{
    return m_name;
}

void Read::setSequence(std::string sequence)
{
    m_sequence = std::move(sequence);
}
const std::string& Read::getSequence() const
{
    return m_sequence;
}

void Read::setDescription(std::string description)
{
    m_description = std::move(description);
}
const std::string& Read::getDescription() const
{
    return m_description;
}

void Read::setQuality(std::string quality)
{
    m_quality = std::move(quality);
}
const std::string& Read::getQuality() const
{
    return m_quality;
}

Read Read::reverseComplement() const
{
    std::string qual;
    qual.reserve(m_quality.size());
    std::copy(m_quality.crbegin(), m_quality.crend(), std::inserter(qual, qual.begin()));
    return Read(m_name, revCompl(m_sequence), m_description, qual);
}
