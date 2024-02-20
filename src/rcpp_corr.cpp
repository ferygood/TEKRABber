#include <Rcpp.h>
using namespace Rcpp;

//' Estimate the correlation between genes and transposable elements
//' 
//' @param df1 First dataframe
//' @param df2 Second dataframe
//' @param Method correlation method
//' @return a dataframe containing correlation results
// [[Rcpp::export]]
DataFrame rcpp_corr(DataFrame df1, DataFrame df2, StringVector Method) {
    
    // get columns name
    StringVector name1 = df1.names();
    StringVector name2 = df2.names();
    
    // create empty vectors for building dataframe
    StringVector gene_name;
    StringVector te_name;
    NumericVector p_value;
    NumericVector coef;
    
    #pragma omp parallel for
    for (int i = 0; i < df1.size(); ++i) {
        
        for (int j = 0; j < df2.size(); ++j) {
            
            NumericVector value1 = df1[i];
            NumericVector value2 = df2[j];
            
            // append gene and te name
            gene_name.push_back(name1[i]);
            te_name.push_back(name2[j]);
            
            // calculate selected correlation method and append results
            List result;
            Function cor("cor.test");
            result = cor(value1, value2, Named("method")=Method);
            p_value.push_back(result["p.value"]);
            coef.push_back(result["estimate"]);
            
        }
    }
    
    DataFrame df = DataFrame::create(Named("GeneName") = gene_name,
                                     Named("TEName") = te_name,
                                     Named("pvalue") = p_value,
                                     Named("coef") = coef);
    
    return df;
    
}
