#ifndef TMAC_INCLUDE_RANGE_H
#define TMAC_INCLUDE_RANGE_H

struct Range {
  int start;  // starting index
  int end;    // pass the end index

  Range(int s, int e) : start(s), end(e) {}
};

#endif
