
#include "types.h"

int compare(FILE_POSITION_TYPE const & a, FILE_POSITION_TYPE const & b) {
  if (a < b) {
    return -1;
  } else if (a > b) {
    return 1;
  } else {
    return 0;
  }
}
