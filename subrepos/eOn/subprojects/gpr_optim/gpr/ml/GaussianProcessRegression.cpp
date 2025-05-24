/*
 * GaussianProcessRegression.cpp
 *
 *  Created on: 3 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "GaussianProcessRegression.h"

#include <stdio.h>

#define EIGEN_NO_DEBUG
#include <Eigen/Core>
#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "../../backend/DistibutionFunctions.h"
#include "../auxiliary/Distance.h"

namespace gpr {

GaussianProcessRegression::GaussianProcessRegression()
{
    sigma2 = 1e-7;
    jitter_sigma2 = 1e-6;
    optimization_alg = SCG_opt;

    lik_gaussian = new LikGaussian;
    const_cov_fun = new ConstantCF;
    sexpat_cov_func = new SexpatCF;

    is_training_cov_matrix_evaluated = false;
    is_decomposed_succesfully = false;

    num_of_potential_calls = 0;
    failedOptimizer = false;
}

GaussianProcessRegression::~GaussianProcessRegression()
{
    if (lik_gaussian != nullptr) {
        delete lik_gaussian;
        lik_gaussian = nullptr;
    }
    if (const_cov_fun != nullptr) {
        delete const_cov_fun;
        const_cov_fun = nullptr;
    }
    if (sexpat_cov_func != nullptr) {
        delete sexpat_cov_func;
        sexpat_cov_func = nullptr;
    }
}

void GaussianProcessRegression::initialize(const InputParameters& parameters,
                                           const AtomsConfiguration& conf_info)
{
    PriorBase prior_parameters;

    sigma2 = parameters.gp_sigma2.value;
    jitter_sigma2 = parameters.jitter_sigma2.value;
    report_level = parameters.report_level.value;

    if (parameters.optimization_alg.value == "SCG_opt") {
        optimization_alg = SCG_opt;
    }
    // Add "else if" here for yet another algorithm
    else {
        optimization_alg = SCG_opt;
    }

    // Initialize Gaussian likelihood
    prior_parameters.setMu(parameters.prior_mu.value);
    prior_parameters.setNu(parameters.prior_nu.value);
    prior_parameters.setS2(parameters.prior_s2.value);
    lik_gaussian->setSigma2(parameters.sigma2.value);

    // Initialize constant covariance function
    //    Field<double> lengthScale; <----- ? is it calculated or is in initial
    //    data?

    sexpat_cov_func->setMagnSigma2(parameters.magnSigma2.value);
    sexpat_cov_func->setConfInfo(conf_info);
    sexpat_cov_func->setPriorParametersSqrtt(prior_parameters);
    sexpat_cov_func->setPriorParametersGaussian(prior_parameters);

    // Initialize constant covariance function
    const_cov_fun->setConstSigma2(parameters.constSigma2.value);

    if (parameters.check_derivative.value == "true")
        opt_alg_settings.check_derivative = true;
    else
        opt_alg_settings.check_derivative = false;
    opt_alg_settings.report_level = parameters.report_level.value;
    opt_alg_settings.max_iter = parameters.max_iter.value;
    opt_alg_settings.tolerance_func = parameters.tolerance_func.value;
    opt_alg_settings.tolerance_sol = parameters.tolerance_sol.value;
    opt_alg_settings.lambda_limit = parameters.lambda_limit.value;
    opt_alg_settings.lambda = parameters.lambda.value;
}

// used to be update()
// calculation for R_all2 should be moved out
void GaussianProcessRegression::setHyperparameters(
    const Observation& all_obs, const AtomsConfiguration& conf_info,
    const bool update_sexpat_cf_param, const bool update_const_cf_param,
    const bool update_sqrt_prior_param)
{
    aux::Distance distance;
    Field<double> dist;
    Field<double> dummy(1, 1);
    double mean_y;
    double range_x;
    double range_y;
    double norm_inv;
    PriorBase prior_parameters;

    math::DistributionFunctions distr_func;
    aux::AuxiliaryFunctionality aux_func;
    Coord R_all2;

    mean_y = all_obs.E.getMean();
    range_y = all_obs.E.getMaxElt() - all_obs.E.getMinElt();

    // range_x = max(max(dist_at(R_all,R_all,conf_info,1)));
    dummy.set(1);
    distance.dist_at(all_obs.R, all_obs.R, conf_info, dummy, dist);
    range_x = dist.getMaxElt();

#ifndef NDEBUG
    assertMsg(const_cov_fun != nullptr,
              "Object of constant covariance function is not allocated!");
    assertMsg(sexpat_cov_func != nullptr,
              "Object of sexpAt covariance function is not allocated!");
#endif

    if (update_const_cf_param)
        const_cov_fun->setConstSigma2(std::max(1., mean_y * mean_y));

    if (update_sexpat_cf_param) {
        if (update_const_cf_param) {
            norm_inv = distr_func.normalCDFInverse(0.75, 0., range_y / 3.);
            sexpat_cov_func->setMagnSigma2(norm_inv * norm_inv);
        }

        norm_inv = distr_func.normalCDFInverse(0.75, 0., range_x / 3.);
        aux_func.repmatConst(1, conf_info.n_pt, norm_inv,
                             sexpat_cov_func->getLengthScaleRef());
    }

    if (update_sqrt_prior_param) {
        prior_parameters = sexpat_cov_func->getPriorParametersSqrtt();
        prior_parameters.setS2(std::max(1., (range_y / 3.) * (range_y / 3.)));
        sexpat_cov_func->setPriorParametersSqrtt(prior_parameters);
    }

    prior_parameters = sexpat_cov_func->getPriorParametersGaussian();
    prior_parameters.setS2(std::max(1., (range_x / 3.) * (range_x / 3.)));
    sexpat_cov_func->setPriorParametersGaussian(prior_parameters);
}

void GaussianProcessRegression::calculateVariance(Observation& image1)
{
    // NOTE: C is a truncated covariance!
    Index_t n = image1.R.getNumRows();
    EigenMatrix R_mod;
    Eigen::VectorXd R_mod_ind;
    EigenMatrix cov_matrix, KK;
    Eigen::VectorXd V;
    EigenMatrix v;
    Eigen::VectorXd VarEG_R;
    Index_t counter = 0;
    aux::AuxiliaryFunctionality aux_func;

    // R2 = [repmat(R,D+1,1),reshape(repmat(0:D,n,1),[],1)];
    aux_func.assembleMatrixOfRepetitiveCoordinates(image1.R, R_mod, R_mod_ind);

    // KK =
    // gp_cov(gp,R_all2,[repmat(R,D+1,1),reshape(repmat(0:D,n,1),[],1)]);
    evaluateCovarianceMatrix(R_matrix, R_mod, R_indices, R_mod_ind, KK);

    // Cov = gp_trcov(gp,R2);
    evaluateTrainingCovarianceMatrix(R_mod, R_mod_ind, cov_matrix);

    // V = diag(Cov);
    V = cov_matrix.diagonal();

    // v = L\KK;
    // v = L.triangularView<Eigen::Upper>().solve(KK);
    v = L.inverse() * KK;

    // VarEG_R = V - sum(v'.*v',2);
    // VarEG_R.resize(v.rows());
    // std::cout << "\nv is " << v.rows() << " " << v.cols() << std::endl;
    // std::cout << "\nL is " << L.rows() << " " << L.cols() << std::endl;
    // std::cout << "\nKK is " << KK.rows() << " " << KK.cols() << std::endl;
    // std::cout << "\ncov_matrix is " << cov_matrix.rows() << " "
    //           << cov_matrix.cols() << std::endl;
    // VarEG_R = V - sum(v'.*v',2);
    VarEG_R.resize(V.size());
    for (Index_t i = 0; i < V.size(); ++i) {
        double tmp = 0.;
        for (Index_t j = 0; j < v.rows(); ++j) {
            tmp += v(j, i) * v(j, i);
        }
        VarEG_R(i) = V(i) - tmp;
    }
    // std::cout << VarEG_R;
    // VarE_R = VarEG_R(1:n,:);
    image1.E.resize(1, n);
    for (Index_t idx = 0; idx < n; ++idx) {
        image1.E[idx] = VarEG_R(idx);
    }

    // G_R = reshape(EG_R((n+1):end,1), n, D);
    // G_R = reshape(EG_R((n+1):end,:), n, size(EG_R,1)/n-1);
    counter = n;
    image1.G.resize(n, image1.R.getNumCols());
    for (Index_t j = 0; j < image1.G.getNumCols(); ++j) {
        for (Index_t i = 0; i < image1.G.getNumRows(); ++i) {
            image1.G(i, j) = VarEG_R(counter++);
        }
    }
}

void GaussianProcessRegression::evaluateTrainingCovarianceMatrix(
    const EigenMatrix& x, const Eigen::VectorXd& x_ind, EigenMatrix& cov_matrix)
{
    Index_t n = (Index_t)x.rows();
    Eigen::VectorXd uDdim;  // unique derivatives' dimensions
    EigenMatrix Ktemp;

    cov_matrix.resize(n, n);
    cov_matrix.setZero();

    extractUniqueIndices(x_ind, uDdim);

    // Apply all covariance functions
    applyCovarianceFunction(x, x_ind, uDdim, *const_cov_fun, Ktemp);
    cov_matrix = cov_matrix + Ktemp;

    applyCovarianceFunction(x, x_ind, uDdim, *sexpat_cov_func, Ktemp);
    cov_matrix = cov_matrix + Ktemp;

    // TODO: check
    // The following lines from the MATLAB code make no sense!
    //    n = size(K,1);
    //    n1 = n+1;
    //    K(1:n1:end)=K(1:n1:end) + gp.jitterSigma2;

    //    C = C + gp.lik.fh.trcov(gp.lik, x1);
    Field<double> C;
    Coord x_coord;
    x_coord.resize((Index_t)x.rows(), (Index_t)x.cols());
    for (Index_t i = 0; i < x_coord.getNumRows(); ++i) {
        for (Index_t j = 0; j < x_coord.getNumCols(); ++j) {
            x_coord(i, j) = x(i, j);
        }
    }

#ifndef NDEBUG
    assertMsg(lik_gaussian != nullptr,
              "Object of Gausian likelihood is not allocated!");
#endif
    lik_gaussian->evaluateTrainingCovarianceMatrix(x_coord, C);

    for (Index_t i = 0; i < C.getNumRows(); ++i) {
        for (Index_t j = 0; j < C.getNumCols(); ++j) {
            cov_matrix(i, j) += C(i, j);
        }
    }

    // log_man << "Training covariance matrix is " << cov_matrix.size() << "
    // over ( " << cov_matrix.rows() << ", " << cov_matrix.cols() << " )" <<
    // "\n";
    is_training_cov_matrix_evaluated = true;
}

void GaussianProcessRegression::evaluateCovarianceMatrix(
    const EigenMatrix& x1, const EigenMatrix& x2, const Eigen::VectorXd& x1_ind,
    const Eigen::VectorXd& x2_ind, EigenMatrix& C)
{
    Index_t n = (Index_t)x1.rows();
    Index_t n1 = (Index_t)x2.rows();
    EigenMatrix Ktemp;
    Eigen::VectorXd uDdim, uDdim2;

    C.resize(n, n1);

    C.setZero();
    Ktemp.setZero();

    extractUniqueIndices(x1_ind, uDdim);
    extractUniqueIndices(x2_ind, uDdim2);

    /* Evaluate all covariance matrices */
    evaluateCovarianceFunction(x1, x2, x1_ind, uDdim, x2_ind, uDdim2,
                               *sexpat_cov_func, Ktemp);
    C += Ktemp;

    evaluateCovarianceFunction(x1, x2, x1_ind, uDdim, x2_ind, uDdim2,
                               *const_cov_fun, Ktemp);
    C += Ktemp;
}

// void GaussianProcessRegression::extractCoordinatesByIndex(
//    const EigenMatrix& x, const Eigen::VectorXd& ind_Ddim,
//    const Index_t ind, Coord& x_loc)
//{
//    // We know that ind_Ddim is a vector of repetitive indices, starting from
//    0,
//    // so we can use this information to determine the start/end rows of x_loc
//    Index_t Ddim_rows = (Index_t)ind_Ddim.rows();
//    Index_t num_rep_Ddim = Ddim_rows / (ind_Ddim(Ddim_rows - 1) + 1);
//
//    Index_t i_start = ind * num_rep_Ddim;
//    Index_t i_end = (ind + 1) * num_rep_Ddim;
//
//    Index_t counter_x_loc = 0;
////    Index_t counter_x = i_start * (Index_t)x.cols();
//
//    for (Index_t i = i_start; i < i_end; ++i) {
//        // Here we assume that the Eigen::Matrix is stored in a row-major
//        // order
//        std::copy(x.data() + i * x.cols(), x.data() + (i + 1) * x.cols(),
//                  x_loc.getInternalVector().data() + counter_x_loc *
//                  x_loc.getNj());
//        x_loc(counter_x_loc, x_loc.getNj() - 1) = ind;
//        ++counter_x_loc;
////        for (Index_t j = 0; j < x.cols(); ++j) {
////            // Here we assume that the Eigen::Matrix is stored in a
/// row-major /            // order /            x_loc[counter_x_loc++] =
/// x.data()[counter_x++]; // instead of x(i, j); /        } /
/// x_loc[counter_x_loc++] = ind;
//    }
//}
//
void GaussianProcessRegression::assignBlockToMatrix(
    const Eigen::VectorXd& ind_Ddim1, const Eigen::VectorXd& ind_Ddim2,
    const Index_t row_val, const Index_t col_val, const Field<double>& field,
    EigenMatrix& matrix, bool transpose_field)
{
    // We know that ind_Ddim1 and ind_Ddim2 are vectors of repetitive indices,
    // starting from 0, so we can use this information to determine the
    // start/end rows and columns of the matrix
    Index_t Ddim1_rows = (Index_t)ind_Ddim1.rows();
    Index_t Ddim2_rows = (Index_t)ind_Ddim2.rows();
    Index_t num_rep_Ddim1 = Ddim1_rows / (ind_Ddim1(Ddim1_rows - 1) + 1);
    Index_t num_rep_Ddim2 = Ddim2_rows / (ind_Ddim2(Ddim2_rows - 1) + 1);

    Index_t i_start = row_val * num_rep_Ddim1;
    Index_t i_end = (row_val + 1) * num_rep_Ddim1;
    Index_t j_start = col_val * num_rep_Ddim2;
    Index_t j_end = (col_val + 1) * num_rep_Ddim2;
    Index_t i_conter = 0;

    if (!transpose_field) {
        for (Index_t i = i_start; i < i_end; ++i) {
            Index_t j_conter = 0;
            for (Index_t j = j_start; j < j_end; ++j) {
                matrix(i, j) = field(i_conter, j_conter);
                ++j_conter;
            }
            ++i_conter;
        }
    } else {
        for (Index_t i = i_start; i < i_end; ++i) {
            Index_t j_conter = 0;
            for (Index_t j = j_start; j < j_end; ++j) {
                matrix(i, j) = field(j_conter, i_conter);
                ++j_conter;
            }
            ++i_conter;
        }
    }
}

void GaussianProcessRegression::extractCoordinatesByIndex(
    const EigenMatrix& x, const Eigen::VectorXd& ind_Ddim, const Index_t ind,
    Coord& x_loc)
{
    Index_t row_counter = 0;

    for (Index_t i = 0; i < x.rows(); ++i) {
        if (ind_Ddim[i] == ind) {
            for (Index_t j = 0; j < x.cols(); ++j) {
                x_loc(row_counter, j) = x(i, j);
            }
            ++row_counter;
        }
    }
}

// void GaussianProcessRegression::assignBlockToMatrix(
//    const Eigen::VectorXd& ind_Ddim1, const Eigen::VectorXd& ind_Ddim2,
//    const Index_t row_val, const Index_t col_val, const Field<double>& field,
//    EigenMatrix& matrix, bool transpose_field)
//{
//    // FIXME: this is a forehead approach to find and replace a block in
//    // a matrix. Should be changed!
//    Index_t i_conter = 0;
//    for (Index_t i = 0; i < matrix.rows(); ++i) {
//        if (ind_Ddim1(i) == row_val) {
//            Index_t j_conter = 0;
//            for (Index_t j = 0; j < matrix.cols(); ++j) {
//                if (ind_Ddim2(j) == col_val) {
//                    if (!transpose_field)
//                        matrix(i, j) = field(i_conter, j_conter);
//                    else
//                        matrix(i, j) = field(j_conter, i_conter);
//                    ++j_conter;
//                }
//            }
//            ++i_conter;
//        }
//    }
//}

void GaussianProcessRegression::removeColumn(const Index_t column_ind,
                                             EigenMatrix& matrix)
{
    Index_t num_rows = (Index_t)matrix.rows();
    Index_t num_cols = (Index_t)matrix.cols() - 1;

    if (column_ind < num_cols)
        matrix.block(0, column_ind, num_rows, num_cols - column_ind) =
            matrix.block(0, column_ind + 1, num_rows, num_cols - column_ind);

    matrix.conservativeResize(num_rows, num_cols);
}

void GaussianProcessRegression::extractUniqueIndices(
    const Eigen::VectorXd& input_ind, Eigen::VectorXd& output_ind)
{
    // TODO: check if this algorithm can be optimized
    std::vector<double> unique_ind_loc;

    unique_ind_loc.insert(unique_ind_loc.begin(), input_ind.data(),
                          input_ind.data() + input_ind.size());

    // erase all non positive values
    // FIXME: optimize
    auto it = std::remove_if(unique_ind_loc.begin(), unique_ind_loc.end(),
                             [](const double value) { return value <= 0; });
    unique_ind_loc.erase(it, unique_ind_loc.end());

    // erase all duplicates and sort
    sort(unique_ind_loc.begin(), unique_ind_loc.end());
    unique_ind_loc.erase(unique(unique_ind_loc.begin(), unique_ind_loc.end()),
                         unique_ind_loc.end());

    Eigen::Map<Eigen::VectorXd> map(unique_ind_loc.data(),
                                    unique_ind_loc.size());
    output_ind = map;
}

void GaussianProcessRegression::evaluateEnergyAndGradient(
    const Eigen::VectorXd& w, const EigenMatrix& x,
    const Eigen::VectorXd& x_ind, const Eigen::VectorXd& y,
    EnergyAndGradient& energy_and_gradient)
{
    setParameters(w);

    *energy_and_gradient.energy = evaluateEnergy(x, x_ind, y);

    if (fabs(*energy_and_gradient.energy) <= DBL_EPSILON) {
        io::ErrorManager err;
        err << "Energy field is empty!"
            << "\n";
        energy_and_gradient.gradient->setZero();
    } else {
        evaluateGradient(x, x_ind, y, *energy_and_gradient.gradient);
    }
}

double GaussianProcessRegression::evaluateEnergy(const EigenMatrix& x,
                                                 const Eigen::VectorXd& x_ind,
                                                 const Eigen::VectorXd& y)
{
    Index_t n = (Index_t)x.rows();
    double zc = 0., edata = 0., eprior = 0.;
    double energy = 0.;

    is_decomposed_succesfully = decomposeCovarianceMatrix(x, x_ind);

    // Test if the matrix is positive definite
    if (!is_decomposed_succesfully) {
        log_man << "Warning! Matrix L is not positive definite"
                << "\n";
        energy = std::numeric_limits<double>::quiet_NaN();
        ;
    } else {
        calculateMeanPrediction(y);

        zc = 0.;
        for (Index_t n = 0; n < L.rows(); ++n) {
            zc += log(L(n, n));
        }
        edata = 0.5 * n * log(2 * M_PI) + zc +
                0.5 * b.dot(b);  // Line 7, zc should be negative

        // Evaluate the prior contribution to the error from covariance
        // functions
#ifndef NDEBUG
        assertMsg(const_cov_fun != nullptr,
                  "Object of Gausian likelihood is not allocated!");
        assertMsg(sexpat_cov_func != nullptr,
                  "Object of Gausian likelihood is not allocated!");
#endif

        eprior -= const_cov_fun->calculateLogPrior();

        eprior -= sexpat_cov_func->calculateLogPrior();

        // Evaluate the prior contribution to the error from Gaussian likelihood
        eprior -= lik_gaussian->evaluateLogPrior();

        energy = edata + eprior;
    }
    //    std::cout << "\n in evaluateEnergy(): \n" << y << "\n\n";
    return energy;
}

void GaussianProcessRegression::evaluateGradient(const EigenMatrix& x,
                                                 const Eigen::VectorXd& x_ind,
                                                 const Eigen::VectorXd& y,
                                                 Eigen::VectorXd& gradient)
{
    EigenMatrix invC;
    Eigen::VectorXd b;
    Eigen::VectorXd uDdim;
    Field<double> gdata;
    Field<double> gprior;

    // Note: the calculateTrainingCovarianceMatrix() method is always called in
    // the evaluateEnergy() method, and we always call for the
    // evaluateGradient() method after evaluateEnergy() being called (see
    // evaluateEnergyAndGradient(). So, there is no point to re-evaluate
    // covariance metrix here. However, we keep it hear guarded by an
    // if-statement for sake of unit tests. Otherwise, the matrix C will not be
    // allocated.
    if (!is_training_cov_matrix_evaluated) {
        evaluateTrainingCovarianceMatrix(x, x_ind, C);
    }

    invC = C.inverse();

    b = C.lu().solve(y);

    // Gradient with respect to covariance function parameters
    extractUniqueIndices(x_ind, uDdim);

    calculateGradientWithCovFunc(x, x_ind, uDdim, b, invC, *const_cov_fun,
                                 gdata, gprior);

    calculateGradientWithCovFunc(x, x_ind, uDdim, b, invC, *sexpat_cov_func,
                                 gdata, gprior);

    // Gradient with respect to Gaussian likelihood function parameters
    // Evaluate the gradient from Gaussian likelihood
    Field<double> gprior_lik;

#ifndef NDEBUG
    assertMsg(lik_gaussian != nullptr,
              "Object of Gausian likelihood is not allocated!");
#endif
    lik_gaussian->evaluateLogPriorGradient(gprior_lik);

    if (!gprior_lik.isEmpty()) {
        gprior_lik *= -1.;
        gprior.append(gprior_lik);
    }

    if (!gprior.isEmpty()) {
        if (gradient.rows() != gprior.getSize())
            gradient.resize(gprior.getSize());
        for (Index_t n = 0; n < gradient.rows(); ++n)
            gradient[n] = gdata[n] + gprior[n];
    }
}

void GaussianProcessRegression::setParameters(const Eigen::VectorXd& w)
{
    if (w.rows() == 0) return;

#ifndef NDEBUG
    assertMsg(sexpat_cov_func != nullptr,
              "Object of Gausian likelihood is not allocated!");
    assertMsg(const_cov_fun != nullptr,
              "Object of Gausian likelihood is not allocated!");
    assertMsg(lik_gaussian != nullptr,
              "Object of Gausian likelihood is not allocated!");
#endif

    sexpat_cov_func->setParameters(w);
    const_cov_fun->setParameters(w);
    lik_gaussian->setParameters(w);
}

bool GaussianProcessRegression::decomposeCovarianceMatrix(
    const EigenMatrix& x, const Eigen::VectorXd& x_ind)
{
    Eigen::LLT<EigenMatrix> llt;

    evaluateTrainingCovarianceMatrix(x, x_ind, C);
    llt.compute(C);
    L = llt.matrixL();

    return llt.info() == Eigen::NumericalIssue ? false : true;
}

void GaussianProcessRegression::calculatePosteriorMeanPrediction()
{
    // This function mimics the following MATLAB code
    // a = L'\(L\[E_all;G_all(:)]);
    a = L.transpose().triangularView<Eigen::Upper>().solve(b);
}

void GaussianProcessRegression::calculateMeanPrediction(
    const Eigen::VectorXd& y)
{
    // This function mimics the following MATLAB code
    // b = L\[E_all;G_all(:)];
    b = L.triangularView<Eigen::Lower>().solve(y);
}

void GaussianProcessRegression::calculatePotential(Observation& image1)
{
    // NOTE: C is a truncated covariance!

    Index_t N_im = image1.R.getNumRows();
    EigenMatrix R_mod;
    Eigen::VectorXd R_mod_ind;
    EigenMatrix KK;
    Eigen::VectorXd EG_R;
    Index_t counter = 0;
    aux::AuxiliaryFunctionality aux_func;

    // [repmat(R,D+1,1),reshape(repmat(0:D,N_im,1),[],1)]
    aux_func.assembleMatrixOfRepetitiveCoordinates(image1.R, R_mod, R_mod_ind);

    // KK =
    // gp_cov(gp,R_all2,[repmat(R,D+1,1),reshape(repmat(0:D,N_im,1),[],1)]);
    evaluateCovarianceMatrix(R_matrix, R_mod, R_indices, R_mod_ind, KK);
    // log_man << "KK is " << KK.size() << " over ( " << KK.rows() << ", " <<
    // KK.cols() << " )" << "\n";

    // EG_R = KK'*a;
    EG_R = KK.transpose() * a;

    // E_R = EG_R(1:N_im,1);
    image1.E.resize(1, N_im);
    for (Index_t n = 0; n < N_im; ++n)
        image1.E[n] = EG_R(n);

    // G_R = reshape(EG_R((N_im+1):end,1),N_im,D);
    counter = N_im;
    image1.G.resize(N_im, image1.R.getNumCols());
    for (Index_t j = 0; j < image1.G.getNumCols(); ++j) {
        for (Index_t i = 0; i < image1.G.getNumRows(); ++i) {
            image1.G(i, j) = EG_R(counter++);
        }
    }

    ++num_of_potential_calls;
}

void GaussianProcessRegression::optimize(const Observation& observation)
{
    // FIXME: add check for the type of minimization algorithm
    aux::AuxiliaryFunctionality aux_func;
    funcmin::SCG scg;
    Eigen::VectorXd parameters_cf;  // assembled parameters of covariance
                                    // functions and likelihood

    aux_func.assembleMatrixOfRepetitiveCoordinates(observation.R, R_matrix,
                                                   R_indices);
    aux_func.assembleVectorFromEnergyAndGradient(observation,
                                                 energy_and_gradient);

    // The order is taken from gp_pak.m from MATLAB.
    // Note that only sexpat_cov_func actually return some meaningful data
    // parameters_cf = const_cov_fun->combineParameters();
    parameters_cf = sexpat_cov_func->combineParameters();
    Eigen::VectorXd old_parameters_cf = sexpat_cov_func->combineParameters();
    // parameters_cf = lik_gaussian->combineParameters();

    // log_man << "Parameters before SCG: " << parameters_cf << "\n";
    scg.setAlgorithmSettings(opt_alg_settings);

    std::chrono::time_point<std::chrono::steady_clock> start;

    if (this->report_level >= 2) {
        start = std::chrono::steady_clock::now();
    }

    scg.optimize(R_matrix, R_indices, energy_and_gradient, parameters_cf,
                 &gpr::GaussianProcessRegression::evaluateEnergyAndGradient,
                 *this);

    if (this->report_level >= 2) {
        log_man << "R_matrix size for optimization " << R_matrix.size()
                << " over ( " << R_matrix.rows() << ", " << R_matrix.cols()
                << " )"
                << "\n";
    }

    if (this->report_level >= 2) {
        std::chrono::duration<double> elp_time =
            std::chrono::steady_clock::now() - start;
        log_man << "optimize time: " << elp_time.count() << "s\n";
    }

    if (scg.failedOptim) {
        log_man << "Set optim failure in GP\n";
        failedOptimizer = true;
    } else {
        log_man << "Optimize success \n";
        failedOptimizer = false;
    }

    setParameters(parameters_cf);
    if (this->report_level >= 2) {
        log_man << "Parameter size " << parameters_cf.size() << "\n";
        log_man << "magnSigma2: " << sexpat_cov_func->getMagnSigma2() << "\n";
        log_man << "lengthScales:\n"
                << sexpat_cov_func->getLengthScaleRef().extractEigenVector()
                << "\n";
        // log_man << "Difference: " << old_parameters_cf - parameters_cf <<
        // "\n";
        log_man << "Parameter Norm: "
                << (old_parameters_cf - parameters_cf).norm() << "\n";
    }
    calculatePosteriorMeanPrediction();
}

} /* namespace gpr */
