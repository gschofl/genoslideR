#include <Rcpp.h>

void reverse_position_cpp( Rcpp::IntegerVector &cut_start,
                           Rcpp::IntegerVector &cut_width,
                           Rcpp::IntegerVector &hit_start,
                           Rcpp::IntegerVector &hit_width,
                           Rcpp::IntegerVector &strand );

// [[Rcpp::export]]
void update_genomic_position_cpp( Rcpp::S4 cut_ranges,
                                  Rcpp::S4 genomic_hit_ranges,
                                  Rcpp::S4 alignment_hit_ranges,
                                  Rcpp::IntegerVector strand)
{
    Rcpp::IntegerVector cut_start = cut_ranges.slot("start");
    Rcpp::IntegerVector cut_width = cut_ranges.slot("width");
    Rcpp::IntegerVector gen_start = genomic_hit_ranges.slot("start");
    Rcpp::IntegerVector gen_width = genomic_hit_ranges.slot("width");
    //Rcpp::Rcout << "Cut_start before rev: " << cut_start[0] << std::endl;
    if ( is_true(any( strand == 2 )) ) {
      reverse_position_cpp( cut_start, cut_width, gen_start, gen_width, strand ); 
    };
    //Rcpp::Rcout << "Cut_start after rev: " << cut_start[0] << std::endl;
    Rcpp::IntegerVector aln_start = alignment_hit_ranges.slot("start");
    Rcpp::IntegerVector shift = aln_start - gen_start;
    //Rcpp::Rcout << "Shift: " << shift[0] << std::endl;
    cut_ranges.slot("start") =  cut_start + shift;
}


// [[Rcpp::export]]
void update_alignment_position_cpp( Rcpp::S4 cut_ranges,
                                    Rcpp::S4 genomic_hit_ranges,
                                    Rcpp::S4 alignment_hit_ranges,
                                    Rcpp::IntegerVector strand)
{
    Rcpp::IntegerVector cut_start = cut_ranges.slot("start");
    Rcpp::IntegerVector cut_width = cut_ranges.slot("width");
    Rcpp::IntegerVector aln_start = alignment_hit_ranges.slot("start");
    Rcpp::IntegerVector aln_width = alignment_hit_ranges.slot("width");
    
    if ( is_true(any( strand == 2 )) ) {
      reverse_position_cpp( cut_start, cut_width, aln_start, aln_width, strand ); 
    };
    
    Rcpp::IntegerVector gen_start = genomic_hit_ranges.slot("start");
    Rcpp::IntegerVector gen_width = genomic_hit_ranges.slot("width");
    Rcpp::IntegerVector shift = gen_start + gen_width - aln_start - aln_width;
    cut_ranges.slot("start") =  cut_start + shift;
}


void reverse_position_cpp( Rcpp::IntegerVector &cut_start,
                           Rcpp::IntegerVector &cut_width,
                           Rcpp::IntegerVector &hit_start,
                           Rcpp::IntegerVector &hit_width,
                           Rcpp::IntegerVector &strand )
{
    for (auto cs = cut_start.begin(), cw = cut_width.begin(),
         hs = hit_start.begin(), hw = hit_width.begin(),
         s = strand.begin();
         cs != cut_start.end(); ++cs, ++cw, ++hs, ++hw, ++s)
    {
        if (*s == 2)
        {
            *cs = *cs + ( 2 * *hs + *hw - 2 * *cs - *cw );
        } 
    }
}  

