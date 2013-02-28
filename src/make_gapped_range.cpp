#include <Rcpp.h>
using namespace Rcpp;

IntegerVector XOR(LogicalVector x, LogicalVector y);
int make_gapped_position(int pos, IntegerVector gapstart, IntegerVector gapwidth);
int make_ungapped_position(int pos, IntegerVector gapstart, IntegerVector gapwidth, bool right);


// [[Rcpp::export]]
Rcpp::S4 make_ungapped_ranges(IntegerVector start,
                              IntegerVector end,
                              Rcpp::S4 gaprange)
{
  Rcpp::S4 range = Rcpp::S4("IRanges");
  IntegerVector gapstart = gaprange.slot("start");
  IntegerVector gapwidth = gaprange.slot("width");
  IntegerVector ungapped_start(0);
  IntegerVector ungapped_end(0);
  IntegerVector::iterator sit, eit;
  
  for(sit = start.begin(), eit = end.begin(); sit != start.end(); ++sit, ++eit)
  {
    ungapped_start.push_back( make_ungapped_position( *sit, gapstart, gapwidth, true) );
    ungapped_end.push_back( make_ungapped_position( *eit, gapstart, gapwidth, false) );
  }
  
  range.slot("start") = ungapped_start;
  range.slot("width") = ungapped_end - ungapped_start + 1;
  return range;
}


// [[Rcpp::export]]
Rcpp::S4 make_gapped_ranges(IntegerVector start,
                            IntegerVector end,
                            IntegerVector gapstart,
                            IntegerVector gapwidth,
                            int aln_len)
{
  Rcpp::S4 range = Rcpp::S4("IRanges");
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


int make_ungapped_position(int pos, IntegerVector gapstart, IntegerVector gapwidth, bool right)
{
  IntegerVector gapend = gapstart + gapwidth - 1;
  LogicalVector preceding_gapend = gapend < pos;
  LogicalVector preceding_gapstart = gapstart <= pos;
  IntegerVector overlapping = XOR(preceding_gapend, preceding_gapstart);
  
  // sum of gapwidths preceding the overlap
  int sum_preceding_gapend = sum(preceding_gapend);
  IntegerVector gapwidth_preceding( gapwidth.begin(), gapwidth.begin() + sum_preceding_gapend );
  int sum_gapwidth_preceding = sum(gapwidth_preceding);

  int novl = overlapping.size();
  int ovl = overlapping[0];
  int out;
  
  if (novl > 0)
  {
    IntegerVector gapstart_overlapping( gapstart.begin() + ovl, gapstart.begin() + ovl + 1 );
  
    if (right)
    {
      out = pos - sum_gapwidth_preceding - ( pos - gapstart_overlapping[0]);
    }
    else
    {
      out = pos - sum_gapwidth_preceding - ( pos - gapstart_overlapping[0] + 1);
    }
  }
  else
  {
    out = pos - sum_gapwidth_preceding;
  }
  return out;
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

IntegerVector XOR(LogicalVector x, LogicalVector y)
{
  int n = x.size();
  IntegerVector out(0);
  for (int i = 0; i < n; i++)
  {
    if (x[i] != y[i])
    {
      out.push_back(i);
    }
  }
  return out;
}
