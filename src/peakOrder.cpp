#include <algorithm>
#include <Rcpp.h>

struct IntervalOrder {

  Rcpp::NumericVector chr;
  Rcpp::NumericVector left;
  Rcpp::NumericVector right;


  bool operator()(int a,int b) {
    bool less;

    if (chr[a] < chr[b]) {
      less = true;
    } else if (chr[a] > chr[b]) {
      less = false;
    } else if (left[a] < left[b]) {
      less = true;
    } else if (left[a] > left[b]) {
      less = false;
    } else if (right[a] < right[b]) {
      less = true;
    } else {
      less = false;
    }
    return less;
  }
};

// [[Rcpp::export]]
Rcpp::RObject peakOrder(SEXP schrom,SEXP sleft,SEXP sright) {
  Rcpp::NumericVector chrom(schrom);
  Rcpp::NumericVector left(sleft);
  Rcpp::NumericVector right(sright);
  int vlen = chrom.size();
  Rcpp::NumericVector index(vlen);
  struct IntervalOrder cmp;
  
  cmp.chr = chrom;
  cmp.left = left;
  cmp.right = right;
  for (int i=0;i<vlen;i++) {
    index[i] = i;
  }
  std::sort(index.begin(),index.end(),cmp);
  for (int i=0;i<vlen;i++) {
    index[i] += 1;
  }
  return index;
}
  
