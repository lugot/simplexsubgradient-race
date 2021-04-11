#ifndef _PARSERS_H_
#define _PARSERS_H_

#include "lp.h"
#include <string>

lp parse_lpmodel(std::string model_name);
lp parse_lpcplex(std::string model_name);

#endif /* _PARSERS_H_*/
