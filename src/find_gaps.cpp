#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP find_gaps( std::string seq ) {
  char gap_char = '-';
  int seq_len = seq.length();
  int gap_len = 0;
  
  bool check_gap_start = true;
  bool check_gap_width = false;
  IntegerVector gap_start;
  IntegerVector gap_width;

  for (int i = 0; i < seq_len; ++i) {  
    if (seq[i] == gap_char) {
      if (check_gap_start == true) {
        gap_start.push_back(i+1);
        check_gap_start = false;
        check_gap_width = true;
      }
      gap_len++;
  } else {
    if (check_gap_width == true) {
      gap_width.push_back(gap_len);
      gap_len = 0;
      check_gap_start = true;
      check_gap_width = false;
      }
    }
  }
    
  if (check_gap_width == true) {
    gap_width.push_back(gap_len);
  }
    
  Rcpp::S4 range = Rcpp::S4("IRanges");
  range.slot("start") = gap_start;
  range.slot("width") = gap_width;
  return range;
}
