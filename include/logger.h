#ifndef INCLUDE_LOGGER_H_
#define INCLUDE_LOGGER_H_

#include <ostream>
#include <string>
#include <vector>

class Logger {
   public:
    void add(int i);
    void add(double d);
    void add(bool b);
    void add(std::string str);
    void add(std::vector<std::string> v);

    friend std::ostream& operator<<(std::ostream& os, const Logger& l);

   private:
    std::vector<std::string> data;
};

#endif  // INCLUDE_LOGGER_H_
