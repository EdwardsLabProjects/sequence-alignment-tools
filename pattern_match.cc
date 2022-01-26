
#include "pattern_match.h"

#if !defined(NO_STD_NAMESPACE)
using namespace std;
#endif

istream & operator>>(istream & is, pattern_list_element & ple) {
  ple.read(is);
  return is;
}

ostream & operator<<(ostream & os, pattern_list_element const & ple) {
  ple.write(os);
  return os;
}

PatternMatch::~PatternMatch() {}
