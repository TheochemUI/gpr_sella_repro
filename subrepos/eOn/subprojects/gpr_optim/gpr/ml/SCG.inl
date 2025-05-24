//
//  SCG.inl
//  gpr_dimer
//
//  Created by Maxim Masterov on 24/11/2020.
//

#ifndef SCG_h
#define SCG_h

#include <algorithm>
#include <cfloat>
#include <cmath>

#include "../../managers/io//LogManager.h"
#include "../../structures/Structures.h"

namespace funcmin {

template <typename ClassName, typename FuncName>
void SCG::optimize(const gpr::EigenMatrix& x, const Eigen::VectorXd& x_ind,
                   const Eigen::VectorXd& y, Eigen::VectorXd& w,
                   FuncName func_to_min, ClassName& holder)
{
    double sigma0 = 1e-4;
    gpr::Index_t func_count = 1;
    gpr::Index_t grad_count = 1;
    gpr::Index_t fail_count = 0;
    double lambdamin = 1.0e-15;
    double lambdamax = 1.0e100;
    gpr::Index_t iter_counter = 0;
    bool success = true;
    this->failedOptim = false;
    gpr::Index_t nsuccess = 0;
    gpr::Index_t nparams = (gpr::Index_t)w.rows();
    gpr::EnergyAndGradient energy_and_gradient;
    gpr::io::LogManager log_man;

    Eigen::VectorXd r, r_new, r_old;  //, r_plus;
    Eigen::VectorXd p;
    Eigen::VectorXd x_new, x_plus;

    double delta = 0., gamma = 0., kappa = 0., Delta = 0.;
    double alpha = 0., mu = 0., beta = 0., sigma = 0.;

    double f_new = 0., f_old = 0.;
    double tmp = 0.;

    uint8_t exit_code = 0;

    // Initial function value and gradient
    energy_and_gradient.energy = &f_old;
    energy_and_gradient.gradient = &r_old;
    (holder.*func_to_min)(w, x, x_ind, y, energy_and_gradient);

    r = r_old;
    p = -r;

    if (settings.report_level >= 1) log_man << "\n";

    // Main optimization loop
    while (iter_counter < settings.max_iter) {
        if (settings.report_level >= 1)
            log_man << " SCG iteration: " << iter_counter << "\n";

        // Step 2
        if (success) {
            mu = p.dot(r);
            if (mu >= 0.) {
                p = -r;
                mu = p.dot(r);
            }

            kappa = p.dot(p);

            if (kappa < DBL_EPSILON) {
                if (settings.report_level >= 2) {
                    log_man << " Gradient smaller than machine precission"
                            << "\n";
                }
                this->failedOptim = true;
                exit_code = 1;
                break;
            }

            sigma = sigma0 / sqrt(kappa);
            x_plus = w + sigma * p;
            energy_and_gradient.energy = &tmp;
            energy_and_gradient.gradient = &r_new;
            (holder.*func_to_min)(x_plus, x, x_ind, y, energy_and_gradient);
            while ((isInfCoeff(r_new) || isNanCoeff(r_new)) &&
                   !std::isnan(f_old)) {
                sigma = 2. * sigma;
                kappa = 0.25 * kappa;
                x_plus = w + sigma * p;
                // [tmp,gplus] = fun(xplus);
                energy_and_gradient.energy = &tmp;
                energy_and_gradient.gradient = &r_new;
                (holder.*func_to_min)(x_plus, x, x_ind, y, energy_and_gradient);
            }
            ++func_count;
            ++grad_count;
            gamma = (p.dot(r_new - r)) / sigma;
        }

        // Increase effective curvature and evaluate step size alpha
        delta = gamma + settings.lambda * kappa;

        if (delta <= 0.) {
            delta = settings.lambda * kappa;
            settings.lambda -= gamma / kappa;
        }

        // Step 5
        alpha = -mu / delta;

        // Calculate the comparison ratio
        x_new = w + alpha * p;
        energy_and_gradient.energy = &f_new;
        energy_and_gradient.gradient = &r_new;
        (holder.*func_to_min)(x_new, x, x_ind, y, energy_and_gradient);
        ++func_count;
        ++grad_count;

        while (std::isinf(f_new) || std::isnan(f_new)) {
            log_man << " Warning! Function value at xnew not finite "
                       "or a number"
                    << "\n";

            settings.lambda = std::min(4. * settings.lambda, lambdamax);
            delta = gamma + settings.lambda * kappa;
            if (delta <= 0) {
                delta = settings.lambda * kappa;
                settings.lambda = settings.lambda - gamma / kappa;
            }
            alpha = -mu / delta;
            x_new = w + alpha * p;
            energy_and_gradient.energy = &f_new;
            energy_and_gradient.gradient = &r_new;
            (holder.*func_to_min)(x_new, x, x_ind, y, energy_and_gradient);
            ++func_count;
            ++grad_count;
            fail_count++;
            if (fail_count > 100) {
            log_man << " Critical! Too many failures"
                    << "\n";
                this->failedOptim = true;
                break;
            }
        }

        // Step 6
        // check data types of f_new and f_old
        Delta = 2 * (f_new - f_old) / (alpha * mu);

        // Step 7
        if (Delta >= 0) {
            success = true;
            ++nsuccess;
            w = x_new;
        } else {
            success = false;
        }

        if (success) {
            if (findMaxAbsValue(p * alpha) < settings.tolerance_sol) {
                exit_code = 2;
                this->failedOptim = false;
                break;
            } else if (std::fabs(f_new - f_old) < settings.tolerance_func) {
                exit_code = 3;
                this->failedOptim = false;
                break;
            } else {
                // Update variables for new position
                f_old = f_new;
                r_old = r;
                r = r_new;

                // If the gradient is zero then we are done.
                if (r.dot(r) < DBL_EPSILON) {  //  && all(isreal(grad_new)
                    this->failedOptim = false;
                    exit_code = 1;
                    break;
                }
            }
        }

        // Step 8
        // Adjust lambda according to comparison ratio.
        if (Delta < 0.25) {
            settings.lambda = std::min(4.0 * settings.lambda, lambdamax);
        }
        if (Delta > 0.75) {
            settings.lambda = std::max(0.5 * settings.lambda, lambdamin);
        }

        // If scale parameter is at its limit, stop optimization
        if (settings.lambda >= settings.lambda_limit) {
            this->failedOptim = true;
            exit_code = 0;
            if (settings.report_level >= 1) {
                log_man << " Warning: Optimization stopped because lambda "
                           "parameter reached limit. Check that the analytic "
                           "gradients are correct!"
                        << "\n";
            }
            break;
        }

        // Update search direction using Polak-Ribiere formula, or re-start
        // in direction of negative gradient after nparams steps.
        if (nsuccess == nparams) {
            p = -r;
            nsuccess = 0;
        } else {
            if (success) {
                beta = (r_old - r).dot(r) / mu;
                p = beta * p - r;
            }
        }

        ++iter_counter;
    }

    // If we get here, then we haven't terminated in the given number of
    // iterations.
    exit_code = 0;
    if (settings.report_level >= 1 && iter_counter == settings.max_iter) {
        this->failedOptim = true;
        log_man << " Maximum number of iterations has been exceeded."
                << "\n";
    }

    if (settings.report_level >= 2) {
        log_man << " Func-count " << func_count << ". Final f(x)=" << f_new
                << "."
                << "\n";
    }

    if (settings.report_level >= 1) log_man << "\n";
}

inline double SCG::findMaxAbsValue(const Eigen::VectorXd& vec)
{
    double res = 0.;

    for (gpr::Index_t n = 0; n < vec.rows(); ++n) {
        double abs_val = fabs(vec(n));
        if (abs_val > res) {
            res = abs_val;
        }
    }

    return res;
}

inline bool SCG::isNanCoeff(const Eigen::VectorXd& vec)
{
    for (gpr::Index_t n = 0; n < vec.rows(); ++n) {
        if (std::isnan(vec(n))) return true;
    }

    return false;
}

inline bool SCG::isInfCoeff(const Eigen::VectorXd& vec)
{
    for (gpr::Index_t n = 0; n < vec.rows(); ++n) {
        if (std::isinf(vec(n))) return true;
    }

    return false;
}

inline void SCG::setAlgorithmSettings(
    const gpr::OptimizationAlgorithmSettings& _settings)
{
    settings = _settings;
}

} /* namespace funcmin */

#endif /* SCG_h */
