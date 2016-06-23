#ifndef READ_H
#define READ_H
#include <string>

class Read
{
public:
    Read();
    Read(std::string);
    Read(std::string, std::string);
    Read(std::string, std::string, std::string);
    Read(std::string, std::string, std::string, std::string);

    void setName(std::string);
    const std::string& getName() const;

    void setSequence(std::string);
    const std::string& getSequence() const;

    void setDescription(std::string);
    const std::string& getDescription() const;

    void setQuality(std::string);
    const std::string& getQuality() const;

    Read reverseComplement() const;

private:
    std::string m_name;
    std::string m_sequence;
    std::string m_description;
    std::string m_quality;
};

#endif
