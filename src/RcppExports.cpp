// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// find_gaps
SEXP find_gaps(std::string seq);
RcppExport SEXP genoslideR_find_gaps(SEXP seqSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< std::string >::type seq(seqSEXP );
        SEXP __result = find_gaps(seq);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// make_gapped_ranges
Rcpp::S4 make_gapped_ranges(IntegerVector start, IntegerVector end, Rcpp::S4 gaprange);
RcppExport SEXP genoslideR_make_gapped_ranges(SEXP startSEXP, SEXP endSEXP, SEXP gaprangeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type start(startSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type end(endSEXP );
        Rcpp::traits::input_parameter< Rcpp::S4 >::type gaprange(gaprangeSEXP );
        Rcpp::S4 __result = make_gapped_ranges(start, end, gaprange);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// make_ungapped_ranges
Rcpp::S4 make_ungapped_ranges(IntegerVector start, IntegerVector end, CharacterVector nm, Rcpp::S4 gaprange);
RcppExport SEXP genoslideR_make_ungapped_ranges(SEXP startSEXP, SEXP endSEXP, SEXP nmSEXP, SEXP gaprangeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type start(startSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type end(endSEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type nm(nmSEXP );
        Rcpp::traits::input_parameter< Rcpp::S4 >::type gaprange(gaprangeSEXP );
        Rcpp::S4 __result = make_ungapped_ranges(start, end, nm, gaprange);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// update_genomic_position_cpp
void update_genomic_position_cpp(Rcpp::S4 cut_ranges, Rcpp::S4 genomic_hit_ranges, Rcpp::S4 alignment_hit_ranges, Rcpp::IntegerVector strand);
RcppExport SEXP genoslideR_update_genomic_position_cpp(SEXP cut_rangesSEXP, SEXP genomic_hit_rangesSEXP, SEXP alignment_hit_rangesSEXP, SEXP strandSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::S4 >::type cut_ranges(cut_rangesSEXP );
        Rcpp::traits::input_parameter< Rcpp::S4 >::type genomic_hit_ranges(genomic_hit_rangesSEXP );
        Rcpp::traits::input_parameter< Rcpp::S4 >::type alignment_hit_ranges(alignment_hit_rangesSEXP );
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type strand(strandSEXP );
        update_genomic_position_cpp(cut_ranges, genomic_hit_ranges, alignment_hit_ranges, strand);
    }
    return R_NilValue;
END_RCPP
}
// update_alignment_position_cpp
void update_alignment_position_cpp(Rcpp::S4 cut_ranges, Rcpp::S4 genomic_hit_ranges, Rcpp::S4 alignment_hit_ranges, Rcpp::IntegerVector strand);
RcppExport SEXP genoslideR_update_alignment_position_cpp(SEXP cut_rangesSEXP, SEXP genomic_hit_rangesSEXP, SEXP alignment_hit_rangesSEXP, SEXP strandSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< Rcpp::S4 >::type cut_ranges(cut_rangesSEXP );
        Rcpp::traits::input_parameter< Rcpp::S4 >::type genomic_hit_ranges(genomic_hit_rangesSEXP );
        Rcpp::traits::input_parameter< Rcpp::S4 >::type alignment_hit_ranges(alignment_hit_rangesSEXP );
        Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type strand(strandSEXP );
        update_alignment_position_cpp(cut_ranges, genomic_hit_ranges, alignment_hit_ranges, strand);
    }
    return R_NilValue;
END_RCPP
}
