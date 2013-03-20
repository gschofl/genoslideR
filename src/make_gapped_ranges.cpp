#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

int make_gapped_position(int pos, IntegerVector gapstart, IntegerVector gapwidth);

// [[Rcpp::export]]
Rcpp::S4 make_gapped_ranges(IntegerVector start,
                            IntegerVector end,
                            Rcpp::S4 gaprange)
{
  Rcpp::S4 range = Rcpp::S4("IRanges");
  Rcpp::S4 ranges = gaprange.slot("ranges");
  Rcpp::S4 seqnames = gaprange.slot("seqnames");
  Rcpp::S4 seqinfo = gaprange.slot("seqinfo");
  
  IntegerVector gapstart = ranges.slot("start");
  IntegerVector gapwidth = ranges.slot("width");
  int seqnames_idx = seqnames.slot("values");
  IntegerVector seqlengths = seqinfo.slot("seqlengths");
  int aln_len = seqlengths[ seqnames_idx ];
  
  int n = start.size();
  int len = gapstart.size();
  int gap_end = gapstart[len-1] + gapwidth[len-1] - 1;
  IntegerVector gapped_start(n);
  IntegerVector gapped_end(n);
  IntegerVector gapped_width(n);
  
  for (int i = 0; i < n; ++i)
  {
    gapped_start[i] = make_gapped_position(start[i], gapstart, gapwidth);
    gapped_end[i] = make_gapped_position(end[i], gapstart, gapwidth);
    gapped_width[i] = gapped_end[i] - gapped_start[i] + 1;
    
    if (gapped_end[i] > aln_len)
    {
      if (gapped_end[i] > gap_end)
      {
        gapped_end[i] -= 1;
      } 
    
      if (gapped_start[i] > gap_end)
      {
        gapped_start[i] -= 1;
        gapped_width[i] = 0;
      }
    }
  }

  range.slot("start") = gapped_start;
  range.slot("width") = gapped_width;
  return range;
}


int make_gapped_position(int pos, IntegerVector gapstart, IntegerVector gapwidth)
{
  int i = 0;
  int n = gapstart.size();
  while (pos >= gapstart[i] && i <= n - 1) {
    pos += gapwidth[i];
    i++;
  }
  return pos;
}

