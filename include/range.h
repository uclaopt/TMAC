#ifndef AROCK_INCLUDE_RANGE_H
#define AROCK_INCLUDE_RANGE_H

struct Range {
  int start;  // starting index
  int end;    // pass the end index

  Range(int s, int e) : start(s), end(e) {}
};

#endif
