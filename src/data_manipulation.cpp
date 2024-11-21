
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <unordered_map>
#include <fstream>
#include <string>
#include <Rinternals.h>
#include <iomanip>
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <ctime>
#include <Rcpp.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace RcppArmadillo;
//[[Rcpp::depends(RcppArmadillo)]]




// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix fastOrder(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(), 
    n = Br.nrow(),
    k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  arma::colvec An =  sum(square(A),1);
  arma::colvec Bn =  sum(square(B),1);
  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  C = sqrt(C);
  arma::umat ordMat = arma::umat(C.n_rows, C.n_cols);
  for(int a = 0; a < C.n_cols; a = a + 1 ) {
    ordMat.col(a) = arma::sort_index(C.col(a), "ascend") + 1;
  }
  return wrap(ordMat); 
}// end func

// [[Rcpp::export]]
NumericMatrix fastPDist(NumericMatrix Ar, NumericMatrix Br) {
  int m = Ar.nrow(), n = Br.nrow(), k = Ar.ncol();
  arma::mat A = arma::mat(Ar.begin(), m, k, false);
  arma::mat B = arma::mat(Br.begin(), n, k, false); 
  arma::colvec An =  sum(square(A),1);
  arma::colvec Bn =  sum(square(B),1);
  arma::mat C = -2 * (A * B.t());
  C.each_col() += An;
  C.each_row() += Bn.t();
  return wrap(sqrt(C)); 
}// end func

// [[Rcpp::export]]
List MCAStep1(NumericMatrix X) {
  arma::mat AM = arma::mat(X.begin(), X.rows(), X.cols(), true);
  arma::colvec rmin = arma::min(AM,1);
  arma::colvec rmax = arma::max(AM,1);
  arma::colvec range = (rmax -rmin);
  AM.each_col() -= rmin;
  AM.each_col() /= range;
  arma::mat FM = join_cols(AM, 1 - AM);
  AM.clear();
  long total = arma::accu(FM);
  arma::rowvec colsum = arma::sum(FM,0);
  arma::colvec rowsum = arma::sum(FM,1);
  FM.each_row() /= sqrt(colsum);
  FM.each_col() /= sqrt(rowsum);
  arma::colvec Dc = 1/(sqrt(rowsum/total));
  return List::create(Named("Z") = wrap(FM), Named("Dc") = wrap(Dc));
}// end func

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List MCAStep2(NumericMatrix Z, NumericMatrix V, NumericVector Dc) {
  arma::mat AV = arma::mat(V.begin(), V.rows(), V.cols(),  true);
  arma::mat AZ = arma::mat(Z.begin(), Z.rows(), Z.cols(),  true);
  arma::colvec ADc = arma::colvec(Dc);
  arma::mat FeaturesCoordinates =  AZ * AV;
  int AZcol = AZ.n_cols;
  AZ.clear();
  FeaturesCoordinates.each_col() %= ADc;
  ADc.clear();
  return List::create(Named("cellsCoordinates") = wrap(std::sqrt(AZcol) * AV), Named("featuresCoordinates") = wrap(FeaturesCoordinates.head_rows(FeaturesCoordinates.n_rows/2)));
}//end func


template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S> &v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}// end func



// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISampling(Eigen::SparseMatrix<double> data, int sample_val, bool upsample = false, bool display_progress=true){
    Progress p(data.outerSize(), display_progress);
    Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
    for (int k=0; k < data.outerSize(); ++k){
      p.increment();
      for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
        double entry = it.value();
        if( (upsample) || (colSums[k] > sample_val)){
          entry = entry * double(sample_val) / colSums[k];
          if (fmod(entry, 1) != 0){
            double rn = R::runif(0,1);
            if(fmod(entry, 1) <= rn){
              it.valueRef() = floor(entry);
            }
            else{
              it.valueRef() = ceil(entry);
            }
          }
          else{
            it.valueRef() = entry;
          }
        }
      }
    }
  return(data);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISamplingPerCell(Eigen::SparseMatrix<double> data, NumericVector sample_val, bool upsample = false, bool display_progress=true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      double entry = it.value();
      if( (upsample) || (colSums[k] > sample_val[k])){
        entry = entry * double(sample_val[k]) / colSums[k];
        if (fmod(entry, 1) != 0){
          double rn = R::runif(0,1);
          if(fmod(entry, 1) <= rn){
            it.valueRef() = floor(entry);
          }
          else{
            it.valueRef() = ceil(entry);
          }
        }
        else{
          it.valueRef() = entry;
        }
      }
    }
  }
  return(data);
}


typedef Eigen::Triplet<double> T;
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> RowMergeMatrices(Eigen::SparseMatrix<double, Eigen::RowMajor> mat1, Eigen::SparseMatrix<double, Eigen::RowMajor> mat2, std::vector< std::string > mat1_rownames,
                                             std::vector< std::string > mat2_rownames, std::vector< std::string > all_rownames){


  // Set up hash maps for rowname based lookup
  std::unordered_map<std::string, int> mat1_map;
  for(unsigned int i = 0; i < mat1_rownames.size(); i++){
    mat1_map[mat1_rownames[i]] = i;
  }
  std::unordered_map<std::string, int> mat2_map;
  for(unsigned int i = 0; i < mat2_rownames.size(); i++){
    mat2_map[mat2_rownames[i]] = i;
  }

  // set up tripletList for new matrix creation
  std::vector<T> tripletList;
  int num_rows = all_rownames.size();
  int num_col1 = mat1.cols();
  int num_col2 = mat2.cols();


  tripletList.reserve(mat1.nonZeros() + mat2.nonZeros());
  for(int i = 0; i < num_rows; i++){
    std::string key = all_rownames[i];
    if (mat1_map.count(key)){
      for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it1(mat1, mat1_map[key]); it1; ++it1){
        tripletList.emplace_back(i, it1.col(), it1.value());
      }
    }
    if (mat2_map.count(key)){
      for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it2(mat2, mat2_map[key]); it2; ++it2){
        tripletList.emplace_back(i, num_col1 + it2.col(), it2.value());
      }
    }
  }
  Eigen::SparseMatrix<double> combined_mat(num_rows, num_col1 + num_col2);
  combined_mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return combined_mat;
}

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> LogLibSizeFactorNorm(Eigen::SparseMatrix<double> data, int scale_factor, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      it.valueRef() = log1p(double(it.value()) / colSums[k] * scale_factor);
      //it.valueRef() = log(( (double(it.value()) / colSums[k]) +pseudo_count)* scale_factor);
    }
  }
  return data;
}// end func


// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> LogSizeFactorNorm(Eigen::SparseMatrix<double> data, Eigen::VectorXd size_factor, int pseudo_count, bool display_progress = true){
  Progress p(data.outerSize(), display_progress);
  
  //Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
  
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      it.valueRef() = log(double(it.value()) / size_factor[k] + pseudo_count);
    }
  }// end for
  return data;
}// end func


/* Performs column scaling and/or centering. Equivalent to using scale(mat, TRUE, apply(x,2,sd)) in R.
 Note: Doesn't handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export(rng = false)]]
NumericMatrix Standardize(Eigen::Map<Eigen::MatrixXd> mat, bool display_progress = true){
  Progress p(mat.cols(), display_progress);
  NumericMatrix std_mat(mat.rows(), mat.cols());
  for(int i=0; i < mat.cols(); ++i){
    p.increment();
    Eigen::ArrayXd r = mat.col(i).array();
    double colMean = r.mean();
    double colSdev = sqrt((r - colMean).square().sum() / (mat.rows() - 1));
    NumericMatrix::Column new_col = std_mat(_, i);
    for(int j=0; j < new_col.size(); j++) {
      new_col[j] = (r[j] - colMean) / colSdev;
    }
  }
  return std_mat;
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastSparseRowScale(Eigen::SparseMatrix<double> mat, bool scale = true, bool center = true,
                                   double scale_max = 10, bool display_progress = true){
  mat = mat.transpose();
  Progress p(mat.outerSize(), display_progress);
  Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double colMean = 0;
    double colSdev = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      colMean += it.value();
    }
    colMean = colMean / mat.rows();
    if (scale == true){
      int nnZero = 0;
      if(center == true){
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          nnZero += 1;
          colSdev += pow((it.value() - colMean), 2);
        }
        colSdev += pow(colMean, 2) * (mat.rows() - nnZero);
      }
      else{
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
        {
          colSdev += pow(it.value(), 2);
        }
      }
      colSdev = sqrt(colSdev / (mat.rows() - 1));
    }
    else{
      colSdev = 1;
    }
    if(center == false){
      colMean = 0;
    }
    Eigen::VectorXd col = Eigen::VectorXd(mat.col(k));
    scaled_mat.col(k) = (col.array() - colMean) / colSdev;
    for(int s=0; s<scaled_mat.col(k).size(); ++s){
      if(scaled_mat(s,k) > scale_max){
        scaled_mat(s,k) = scale_max;
      }
    }
  }
  return scaled_mat.transpose();
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastSparseRowScaleWithKnownStats(Eigen::SparseMatrix<double> mat, NumericVector mu, NumericVector sigma, bool scale = true, bool center = true,
                                   double scale_max = 10, bool display_progress = true){
    mat = mat.transpose();
    Progress p(mat.outerSize(), display_progress);
    Eigen::MatrixXd scaled_mat(mat.rows(), mat.cols());
    for (int k=0; k<mat.outerSize(); ++k){
        p.increment();
        double colMean = 0;
        double colSdev = 1;
        if (scale == true){
            colSdev = sigma[k];
        }
        if(center == true){
            colMean = mu[k];
        }
        Eigen::VectorXd col = Eigen::VectorXd(mat.col(k));
        scaled_mat.col(k) = (col.array() - colMean) / colSdev;
        for(int s=0; s<scaled_mat.col(k).size(); ++s){
            if(scaled_mat(s,k) > scale_max){
                scaled_mat(s,k) = scale_max;
            }
        }
    }
    return scaled_mat.transpose();
}

/* Note: May not handle NA/NaNs in the same way the R implementation does, */

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastCov(Eigen::MatrixXd mat, bool center = true){
  if (center) {
    mat = mat.rowwise() - mat.colwise().mean();
  }
  Eigen::MatrixXd cov = (mat.adjoint() * mat) / double(mat.rows() - 1);
  return(cov);
}

// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastCovMats(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center = true){
  if(center){
    mat1 = mat1.rowwise() - mat1.colwise().mean();
    mat2 = mat2.rowwise() - mat2.colwise().mean();
  }
  Eigen::MatrixXd cov = (mat1.adjoint() * mat2) / double(mat1.rows() - 1);
  return(cov);
}

/* Note: Faster than the R implementation but is not in-place */
// [[Rcpp::export(rng = false)]]
Eigen::MatrixXd FastRBind(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2){
  Eigen::MatrixXd mat3(mat1.rows() + mat2.rows(), mat1.cols());
  mat3 << mat1, mat2;
  return(mat3);
}

/* Calculates the row means of the logged values in non-log space */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastExpMean(Eigen::SparseMatrix<double> mat, bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowmeans(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene means" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += expm1(it.value());
    }
    rm = rm / ncols;
    rowmeans[k] = log1p(rm);
  }
  return(rowmeans);
}


/* use this if you know the row means */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat,
                            NumericVector mu,
                            bool display_progress){
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variances" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  NumericVector allVars = no_init(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      nZero -= 1;
      colSum += pow(it.value() - mu[k], 2);
    }
    colSum += pow(mu[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

/* standardize matrix rows using given mean and standard deviation,
   clip values larger than vmax to vmax,
   then return variance for each row */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat,
                              NumericVector mu,
                              NumericVector sd,
                              double vmax,
                              bool display_progress){
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating feature variances of standardized and clipped values" << std::endl;
  }
  mat = mat.transpose();
  NumericVector allVars(mat.cols());
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    if (sd[k] == 0) continue;
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      nZero -= 1;
      colSum += pow(std::min(vmax, (it.value() - mu[k]) / sd[k]), 2);
    }
    colSum += pow((0 - mu[k]) / sd[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

/* Calculate the variance to mean ratio (VMR) in non-logspace (return answer in
log-space) */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastLogVMR(Eigen::SparseMatrix<double> mat,  bool display_progress){
  int ncols = mat.cols();
  Eigen::VectorXd rowdisp(mat.rows());
  mat = mat.transpose();
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variance to mean ratios" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    double v = 0;
    int nnZero = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += expm1(it.value());
    }
    rm = rm / ncols;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      v += pow(expm1(it.value()) - rm, 2);
      nnZero += 1;
    }
    v = (v + (ncols - nnZero) * pow(rm, 2)) / (ncols - 1);
    rowdisp[k] = log(v/rm);

  }
  return(rowdisp);
}

/* Calculates the variance of rows of a matrix */
// [[Rcpp::export(rng = false)]]
NumericVector FastRowVar(Eigen::Map<Eigen::MatrixXd> x){
  NumericVector out(x.rows());
  for(int i=0; i < x.rows(); ++i){
    Eigen::ArrayXd r = x.row(i).array();
    double rowMean = r.mean();
    out[i] = (r - rowMean).square().sum() / (x.cols() - 1);
  }
  return out;
}// end func

// [[Rcpp::export(rng = false)]]
NumericVector FastRowMean(Eigen::Map<Eigen::MatrixXd> x){
  NumericVector out(x.rows());
  for(int i=0; i < x.rows(); ++i){
    Eigen::ArrayXd r = x.row(i).array();
    out[i] =  r.mean();
  }// end for
  return out;
}// end func

/* Calculate the variance in non-logspace (return answer in non-logspace) */
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastSparseRowVar(Eigen::SparseMatrix<double> mat){
  int ncols = mat.cols();
  Eigen::VectorXd rowdisp(mat.rows());
  mat = mat.transpose();
  
  bool display_progress = true;
  if(display_progress == true){
    Rcpp::Rcerr << "Calculating gene variances" << std::endl;
  }
  Progress p(mat.outerSize(), display_progress);
  for (int k=0; k<mat.outerSize(); ++k){
    p.increment();
    double rm = 0;
    double v = 0;
    int nnZero = 0;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      rm += (it.value());
    }
    rm = rm / ncols;
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
      v += pow((it.value()) - rm, 2);
      nnZero += 1;
    }
    v = (v + (ncols - nnZero) * pow(rm, 2)) / (ncols - 1);
    rowdisp[k] = v;
  }
  return(rowdisp);
}// end func


//Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.rows());
// [[Rcpp::export(rng = false)]]
Eigen::VectorXd FastSparseRowMean(Eigen::SparseMatrix<double> mat){
  Eigen::VectorXd rowSums = mat * Eigen::VectorXd::Ones(mat.cols());
  return(rowSums/mat.cols());
}// end func

//cols_idx should be 0-indexed
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> ReplaceColsC(Eigen::SparseMatrix<double> mat, NumericVector col_idx, Eigen::SparseMatrix<double> replacement){
  int rep_idx = 0;
  for(auto const &ci : col_idx){
    mat.col(ci) = replacement.col(rep_idx);
    rep_idx += 1;
  }
  return(mat);
}


// [[Rcpp::export(rng = false)]]
List GraphToNeighborHelper(Eigen::SparseMatrix<double> mat) {
  mat = mat.transpose();
  //determine the number of neighbors
  int n = 0;
  for(Eigen::SparseMatrix<double>::InnerIterator it(mat, 0); it; ++it) {
    n += 1;
  }
  Eigen::MatrixXd nn_idx(mat.rows(), n);
  Eigen::MatrixXd nn_dist(mat.rows(), n);

  for (int k=0; k<mat.outerSize(); ++k){
    int n_k = 0;
    std::vector<double> row_idx;
    std::vector<double> row_dist;
    row_idx.reserve(n);
    row_dist.reserve(n);
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      if (n_k > (n-1)) {
        Rcpp::stop("Not all cells have an equal number of neighbors.");
      }
      row_idx.push_back(it.row() + 1);
      row_dist.push_back(it.value());
      n_k += 1;
    }
    if (n_k != n) {
      Rcpp::Rcout << n << ":::" << n_k << std::endl;
      Rcpp::stop("Not all cells have an equal number of neighbors.");
    }
    //order the idx based on dist
    std::vector<size_t> idx_order = sort_indexes(row_dist);
    for(int i = 0; i < n; ++i) {
      nn_idx(k, i) = row_idx[idx_order[i]];
      nn_dist(k, i) = row_dist[idx_order[i]];
    }
  }
  List neighbors = List::create(nn_idx, nn_dist);
  return(neighbors);
}

/* Calculates the variance of rows of a matrix */
// [[Rcpp::export(rng = false)]]
NumericVector RowVarRcpp(Eigen::Map<Eigen::MatrixXd> x){
  NumericVector out(x.rows());
  for(int i=0; i < x.rows(); ++i){
    Eigen::ArrayXd r = x.row(i).array();
    double rowMean = r.mean();
    out[i] = (r - rowMean).square().sum() / (x.cols() - 1);
  }
  return out;
}// end func


///#######################################################################
typedef Eigen::Triplet<double> T;

//' ComputeSNN of a given matrix
//'
//' @param nn_ranked A symmetric matrix
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> ComputeSNN(Eigen::MatrixXd nn_ranked, double prune) {
  std::vector<T> tripletList;
  int k = nn_ranked.cols();
  tripletList.reserve(nn_ranked.rows() * nn_ranked.cols());
  for(int j=0; j<nn_ranked.cols(); ++j){
    for(int i=0; i<nn_ranked.rows(); ++i) {
      tripletList.push_back(T(i, nn_ranked(i, j) - 1, 1));
    }
  }
  Eigen::SparseMatrix<double> SNN(nn_ranked.rows(), nn_ranked.rows());
  SNN.setFromTriplets(tripletList.begin(), tripletList.end());
  SNN = SNN * (SNN.transpose());
  for (int i=0; i < SNN.outerSize(); ++i){
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it){
      it.valueRef() = it.value()/(k + (k - it.value()));
      if(it.value() < prune){
        it.valueRef() = 0;
      }
    }
  }
  SNN.prune(0.0); // actually remove pruned values
  return SNN;
}



// [[Rcpp::export(rng = false)]]
void WriteEdgeFile(Eigen::SparseMatrix<double> snn, String filename, bool display_progress){
  if (display_progress == true) {
    Rcpp::Rcerr << "Writing SNN as edge file" << std::endl;
  }
  // Write out lower triangle
  std::ofstream output;
  output.open(filename);
  Progress p(snn.outerSize(), display_progress);
  for (int k=0; k < snn.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, k); it; ++it){
      if(it.col() >= it.row()){
        continue;
      }
      output << std::setprecision(15) << it.col() << "\t" << it.row() << "\t" << it.value() << "\n";
    }
  }
  output.close();
}

// Wrapper function so that we don't have to go back into R before writing to file
// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> DirectSNNToFile(Eigen::MatrixXd nn_ranked,
                                            double prune, bool display_progress,
                                            String filename) {
  Eigen::SparseMatrix<double> SNN = ComputeSNN(nn_ranked, prune);
  WriteEdgeFile(SNN, filename, display_progress);
  return SNN;
}


// [[Rcpp::export]]
std::vector<double> SNN_SmallestNonzero_Dist(
    Eigen::SparseMatrix<double> snn,
    Eigen::MatrixXd mat,
    int n,
    std::vector<double> nearest_dist
) {
  std::vector<double> results;
  for (int i=0; i < snn.outerSize(); ++i){
    // create vectors to store the nonzero snn elements and their indices
    std::vector<double> nonzero;
    std::vector<size_t> nonzero_idx;
    for (Eigen::SparseMatrix<double>::InnerIterator it(snn, i); it; ++it) {
      nonzero.push_back(it.value());
      nonzero_idx.push_back(it.row());
    }
    std::vector<size_t> nonzero_order = sort_indexes(nonzero);
    int n_i = n;
    if (n_i > nonzero_order.size()) n_i = nonzero_order.size();
    std::vector<double> dists;
    for (int j = 0; j < nonzero_order.size(); ++j) {
      // compute euclidean distances to cells with small edge weights
      // if multiple entries have same value as nth element, calc dist to all
      size_t cell = nonzero_idx[nonzero_order[j]];
      if(dists.size() < n_i  || nonzero[nonzero_order[j]] == nonzero[nonzero_order[n_i-1]]) {
        double res = (mat.row(cell) - mat.row(i)).norm();
        if (nearest_dist[i] > 0) {
          res = res - nearest_dist[i];
          if (res < 0) res = 0;
        }
        dists.push_back(res);
      } else {
        break;
      }
    }
    double avg_dist;
    if (dists.size() > n_i) {
      std::sort(dists.rbegin(), dists.rend());
      avg_dist = std::accumulate(dists.begin(), dists.begin() + n_i, 0.0) / n_i;
    } else {
      avg_dist = std::accumulate(dists.begin(), dists.end(), 0.0) / dists.size();
    }
    results.push_back(avg_dist);
  }
  return results;
}// end funcs


//' Graph Laplacian calculation
//'
//' Calculate graph Laplacian of a symmetrix matrix
//'
//' @param A symmetric matrix
//' @export
// [[Rcpp::export]]
arma::mat LaplacianMat(arma::mat A, double sigma, double alpha) {
  // gaussian kernel 
  A = exp(-A/(sigma*A.max()) );
  //A = exp(-A/(2*sigma) );
  //A.elem( find(A > 0.5) ).ones();
  // quantity D^{1/2}
  arma::rowvec D_row = pow(sum(A), -0.5);
  A.each_row() %= D_row;
  arma::colvec D_col = conv_to< colvec >::from(D_row);
  A.each_col() %= D_col;
  arma::mat B = eye(A.n_cols, A.n_cols) - alpha*A;
  return(B);
  //arma::vec eigval;
  //arma::mat eigvec;
  
  //arma::eig_sym(eigval, eigvec, B);
  
  //return(reverse(eigvec, 1));
}// end func




// [[Rcpp::export]]
SEXP MultiplyMat(SEXP Xin, SEXP Yin){
  try{	
    mat X = as<mat>(Xin);
    mat Y = as<mat>(Yin);
    
    return List::create(Named("XY") = X*Y );
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
  }
  return R_NilValue;
}



//' Compute Euclidean distance matrix by columns
//'
//' Used in sc3-funcs.R distance matrix calculation
//' and within the consensus clustering.
//'
//' @param x A numeric matrix.
// [[Rcpp::export]]
Rcpp::NumericMatrix ED2(const Rcpp::NumericMatrix & x) {
  unsigned int outcols = x.ncol(), i = 0, j = 0;
  double d;
  Rcpp::NumericMatrix out(outcols, outcols);
  
  for (j = 0; j < outcols - 1; j++) {
    Rcpp::NumericVector v1 = x.column(j);
    for (i = j + 1; i < outcols; i++) {
      d = sqrt(sum(pow(v1 - x.column(i), 2.0)));
      out(i, j) = d;
      out(j, i) = d;
    }// end for
  }//end for
  
  return out;
}// end func


//' Euclidean distance calculation
//'
//' Calculate graph Laplacian of a symmetrix matrix
//'
//' @param A symmetric matrix
//' @export
// [[Rcpp::export]]
arma::mat parsesED2(arma::mat X, int num_core) {
  // gaussian kernel
  unsigned int num_col = X.n_cols;
  
  // main loop
  arma::mat out = arma::zeros<arma::mat>(0.5*num_col*(num_col - 1), 3);
  size_t irow = 0;
  
  // main loop
  #pragma omp parallel for num_threads(num_core)
  for(size_t j = 0; j < num_col - 1; j++) {
    arma::vec v1 = X.col(j);
    for (size_t i = j + 1; i < num_col; i++) {
      double d = sqrt(sum(pow(v1 - X.col(i), 2.0)));
      out(irow, 0) = j;
      out(irow, 1) = i;
      out(irow, 2) = d;
      irow++;
    }// end for
  }//end for
    return(out);
}// end func


/* Calculate the euclidean distance for data X and data Y */
// [[Rcpp::export]]
SEXP DistEuc(SEXP Xin, SEXP Yin){// *
  try{
    // X -- N x M
    // Y -- K x M
    // D -- N x K
    mat X = as<mat>(Xin);
    mat Y = as<mat>(Yin);
    const int N = X.n_rows;
    const int K = Y.n_rows;
    vec XX = zeros<vec>(N);
    vec YY = zeros<vec>(K);
    
    mat XtY = zeros<mat>(N,K);
    mat D = zeros<mat>(N,K);
    
    XX = sum(X%X,1);
    YY = sum(Y%Y,1);
    
    XtY = X*Y.t();
    D = kron(XX, ones<mat>(1,K));
    D = D + kron(ones<vec>(N),YY.t());
    D = D - 2*XtY;
    
    return List::create(Named("D") = sqrt(D) );
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
  }
  return R_NilValue;
}// end function

//' Consensus matrix computation
//'
//' Computes consensus matrix given cluster labels
//'
//' @param dat a matrix containing clustering solutions in columns
// [[Rcpp::export]]
arma::mat consmx(const arma::mat dat) {
  
  mat res = dat.n_cols * eye<mat>( dat.n_rows, dat.n_rows );
  
  int i, j, k;
  for (j = 0; j < dat.n_cols; j++) {
    for (i = 0; i < dat.n_rows; i++) {
      for (k = i + 1; k < dat.n_rows; k++) {
        if (dat(i, j) == dat(k, j)) {
          res(i, k)++;
          res(k, i)++;
        }
      }
    }
  }
  res /= dat.n_cols;
  return res;
}// end func


//' FastMap for dimension reduction, SIGMOD '95: Proceedings of the 1995 ACM SIGMOD international conference on Management of data
//' June 1995 Pages 163–174; DOI: https://doi.org/10.1145/223784.223812
//'
//' Computes reduced dimension with linear computation time
//'
//' @param D Distance matrix
//' 
//' @param dat a matrix containing clustering solutions in columns
// [[Rcpp::export]]
arma::mat FastMapEculidean(arma::mat D, int num_dim) {
  
  int num_col = D.n_cols; // number of samples
  //int num_row = D.n_rows; // number of genes
  // compute the distance and store as sparse format
  //X = parsesED2(X, num_core);
  // end fi
  
  arma::mat Dold = D;
  arma::mat Dnew = arma::zeros<arma::mat>(num_col, num_col);
  arma::mat DR = arma::zeros<arma::mat>(num_col, num_dim);
  
  for(size_t k = 0; k < num_dim; k++){
    int ida = max(Dold, 0).index_max(); 
    int idb = Dold.col(ida).index_max();

    double dab2 = D(ida, idb)*D(ida, idb);
    if(dab2 > (sqrt(123*2.220446e-16))){
      double dab = sqrt(dab2);
      arma::vec ida_vec = D.col(ida);
      arma::vec idb_vec = D.col(idb);

      DR.col(k) = (ida_vec%ida_vec + dab2 - idb_vec%idb_vec)/(2.0*dab);
    }//end fi
    
    // update Dnew
    for(size_t i = 0; i < num_col; i++){
      for(size_t j = 0; j < num_col; j++){
        double theval = sqrt( abs(Dold(i, j)*Dold(i, j) - pow(DR(i, k) - DR(j, k), 2.0)) );
        Dnew(i, j) = theval;
      }// end for
    }// end for
    Dold = Dnew;
  }// end for
  
  // return dimension reduction results
  return(DR);
}// end func




/* Calculate the euclidean distance for data X and data Y */
// [[Rcpp::export]]
SEXP DistEigenTestRcpp(SEXP Xin){// *
  try{
    // X -- N x M
    // Y -- K x M
    // D -- N x K
    mat D = as<mat>(Xin);
    //arma::mat D = ED2(X);
    
    arma::vec eigval;
    arma::mat eigvec;
    
    arma::eig_sym(eigval, eigvec, D);
    
    return List::create(Named("eigval") = reverse(eigval), Named("eigvec") = eigvec, Named("eigvec2") = reverse(eigvec, 1));
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
  }
  return R_NilValue;
}// end function
