#ifndef INCLUDE_LOGLINE_H_
#define INCLUDE_LOGLINE_H_

#include <assert.h>

#include <iomanip>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

class LogLine {
   public:
    LogLine() {
        this->delim = ',';
        this->fancy = false;
        this->data = std::vector<std::string>(0);
    }
    LogLine(char delim, bool fancy) {
        this->delim = delim;
        this->fancy = fancy;
        this->data = std::vector<std::string>(0);
    }

    void setDelim(char delim) {
        this->delim = delim;
    }
    void setFancy(bool fancy) {
        this->fancy = fancy;
    }

    friend LogLine& operator>>(LogLine& line, int i) {
        line.data.push_back(std::to_string(i));
        return line;
    }
    friend LogLine& operator>>(LogLine& line, double d) {
        std::stringstream ss;
        ss << d;
        line.data.push_back(ss.str());
        return line;
    }
    friend LogLine& operator>>(LogLine& line, bool b) {
        line.data.push_back(b ? "1" : "0");
        return line;
    }
    friend LogLine& operator>>(LogLine& line, std::string str) {
        line.data.push_back(str);
        return line;
    }
    friend LogLine& operator>>(LogLine& line, std::vector<std::string> v) {
        line.data.insert(line.data.end(), v.begin(), v.end());
        return line;
    }

    friend std::ostream& operator<<(std::ostream& os, const LogLine& line) {
        assert(line.data.size() >= 10);

        if (!line.fancy) {
            std::vector<std::string>::const_iterator it;
            for (it = line.data.begin(); next(it) != line.data.end(); ++it) {
                os << *it << line.delim;
            }
            os << *it << std::endl;
        } else {
            os << "Iteration ";
            os << line.data[1] << "\n";

            os << "\t"
               << "dual: " << line.data[4] << "\n"
               << "\t"
               << "primal: " << line.data[5] << "\n";
            os << "\t"
               << "conditional: " << line.data[2] << " | deflected: " << line.data[3]
               << "\n";
            os << "\t"
               << "lambda: " << line.data[6] << " (" << line.data[9]
               << ") | mu: " << line.data[7] << " | delta: " << line.data[8]
               << "\n";
            if (line.data.size() == 11) {
                os << "\t"
                   << "dist to ustar: " << line.data[10] << std::endl;
            }
        }

        return os;
    }

   private:
    std::vector<std::string> data;
    char delim;
    bool fancy;
};

#endif  // INCLUDE_LOGLINE_H_
