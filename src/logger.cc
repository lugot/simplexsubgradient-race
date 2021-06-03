#include "../include/logger.h"

#include <iomanip>

#include "../include/globals.h"

void Logger::add(int i) {
    data.push_back(std::to_string(i));
}

void Logger::add(double d) {
    std::stringstream ss;
    ss << d;

    data.push_back(ss.str());
}

void Logger::add(bool b) {
    if (b) {
        data.push_back("1");
    } else {
        data.push_back("0");
    }
}

void Logger::add(std::string s) {
    data.push_back(s);
}

void Logger::add(std::vector<std::string> v) {
    data.insert(data.end(), v.begin(), v.end());
}

std::ostream& operator<<(std::ostream& os, const Logger& l) {
    std::vector<std::string>::const_iterator it;
    for (it = l.data.begin(); next(it) != l.data.end(); ++it) {
        os << *it << DELIM;
    }
    os << *it << std::endl;

    return os;
}
