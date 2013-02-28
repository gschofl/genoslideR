// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// find_gaps
SEXP find_gaps(std::string seq);
RcppExport SEXP genoslideR_find_gaps(SEXP seqSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    std::string seq = Rcpp::as<std::string >(seqSEXP);
    SEXP __result = find_gaps(seq);
    return Rcpp::wrap(__result);
END_RCPP
}
// make_ungapped_ranges
Rcpp::S4 make_ungapped_ranges(IntegerVector start, IntegerVector end, Rcpp::S4 gaprange);
RcppExport SEXP genoslideR_make_ungapped_ranges(SEXP startSEXP, SEXP endSEXP, SEXP gaprangeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    IntegerVector start = Rcpp::as<IntegerVector >(startSEXP);
    IntegerVector end = Rcpp::as<IntegerVector >(endSEXP);
    Rcpp::S4 gaprange = Rcpp::as<Rcpp::S4 >(gaprangeSEXP);
    Rcpp::S4 __result = make_ungapped_ranges(start, end, gaprange);
    return Rcpp::wrap(__result);
END_RCPP
}
// make_gapped_ranges
Rcpp::S4 make_gapped_ranges(IntegerVector start, IntegerVector end, IntegerVector gapstart, IntegerVector gapwidth, int aln_len);
RcppExport SEXP genoslideR_make_gapped_ranges(SEXP startSEXP, SEXP endSEXP, SEXP gapstartSEXP, SEXP gapwidthSEXP, SEXP aln_lenSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    IntegerVector start = Rcpp::as<IntegerVector >(startSEXP);
    IntegerVector end = Rcpp::as<IntegerVector >(endSEXP);
    IntegerVector gapstart = Rcpp::as<IntegerVector >(gapstartSEXP);
    IntegerVector gapwidth = Rcpp::as<IntegerVector >(gapwidthSEXP);
    int aln_len = Rcpp::as<int >(aln_lenSEXP);
    Rcpp::S4 __result = make_gapped_ranges(start, end, gapstart, gapwidth, aln_len);
    return Rcpp::wrap(__result);
END_RCPP
}
