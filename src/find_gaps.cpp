#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
SEXP find_gaps( std::string seq ) {
  Rcpp::S4 range = Rcpp::S4("IRanges");
  char gap_char = '-';
  int gap_len = 0;
  int gap_pos = 1;
  bool check_gap_start = true;
  bool check_gap_width = false;
  vector<int> gap_start;
  vector<int> gap_width;
  for (auto it = begin(seq); it != end(seq); ++it)
  {  
    if (*it == gap_char)
    {
      if (check_gap_start == true)
      {
        gap_start.push_back(gap_pos);
        check_gap_start = false;
        check_gap_width = true;
      }
      ++gap_len;
    } else {
      if (check_gap_width == true)
      {
        gap_width.push_back(gap_len);
        gap_len = 0;
        check_gap_start = true;
        check_gap_width = false;
      }
    }
    ++gap_pos;
  }
    
  if (check_gap_width == true)
  {
    gap_width.push_back(gap_len);
  }
    
  range.slot("start") = gap_start;
  range.slot("width") = gap_width;
  return range;
}

