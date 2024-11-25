#ifndef DATA_MANIPULATION
#define DATA_MANIPULATION

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>
#include <iomanip>

using namespace Rcpp;

//----------------------------------------------------
Eigen::SparseMatrix<double> RunUMISampling(Eigen::SparseMatrix<double> data, int sample_val,
                                           bool upsample, bool display_progress);
Eigen::SparseMatrix<double> RunUMISamplingPerCell(Eigen::SparseMatrix<double> data,
                                                  NumericVector sample_val, bool upsample,
                                                  bool display_progress);
Eigen::SparseMatrix<double> RowMergeMatrices(Eigen::SparseMatrix<double, Eigen::RowMajor> mat1,
                                             Eigen::SparseMatrix<double, Eigen::RowMajor> mat2,
                                             std::vector< std::string > mat1_rownames,
                                             std::vector< std::string > mat2_rownames,
                                             std::vector< std::string > all_rownames);
Eigen::SparseMatrix<double> LogLibSizeFactorNorm(Eigen::SparseMatrix<double> data, int scale_factor,
                                    bool display_progress );
Eigen::SparseMatrix<double> LogSizeFactorNorm(Eigen::SparseMatrix<double> data, 
                                              Eigen::VectorXd size_factor, 
                                              int pseudo_count, 
                                              bool display_progress );
NumericMatrix Standardize(const Eigen::Map<Eigen::MatrixXd> mat, bool display_progress);
Eigen::MatrixXd FastSparseRowScale(Eigen::SparseMatrix<double> mat, bool scale, bool center,
                                   double scale_max, bool display_progress);
Eigen::MatrixXd FastCov(Eigen::MatrixXd mat, bool center);
Eigen::MatrixXd FastCovMats(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2, bool center);
Eigen::MatrixXd FastRBind(Eigen::MatrixXd mat1, Eigen::MatrixXd mat2);
Eigen::VectorXd FastExpMean(Eigen::MatrixXd mat, bool display_progress);
Eigen::VectorXd FastRowMean(Eigen::MatrixXd mat, bool display_progress);
Eigen::VectorXd FastLogVMR(Eigen::SparseMatrix<double> mat, bool display_progress);
Eigen::VectorXd FastExpVar(Eigen::SparseMatrix<double> mat, bool display_progress);
Eigen::VectorXd SparseRowVar(Eigen::SparseMatrix<double> mat, bool display_progress);
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat,
                            NumericVector mu,
                            bool display_progress);
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat,
                              NumericVector mu,
                              NumericVector sd,
                              double vmax,
                              bool display_progress);
NumericVector RowVarRcpp(Eigen::Map<Eigen::MatrixXd> x);
template <typename S>
std::vector<size_t> sort_indexes(const std::vector<S> &v);
List GraphToNeighborHelper(Eigen::SparseMatrix<double> mat);
//----------------------------------------------------
NumericMatrix fastOrder(NumericMatrix Ar, NumericMatrix Br);
NumericMatrix fastPDist(NumericMatrix Ar, NumericMatrix Br);
List MCAStep1(NumericMatrix X);
List MCAStep2(NumericMatrix Z, NumericMatrix V, NumericVector Dc);
//----------------------------------------------------
Eigen::SparseMatrix<double> ComputeSNN(Eigen::MatrixXd nn_ranked);
void WriteEdgeFile(Eigen::SparseMatrix<double> snn, String filename, bool display_progress);
Eigen::SparseMatrix<double> DirectSNNToFile(Eigen::MatrixXd nn_ranked, double prune, bool display_progress, String filename);
//----------------------------------------------------

#endif//DATA_MANIPULATION

