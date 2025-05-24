/*
 * Distance.cpp
 *
 *  Created on: 30 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "Distance.h"

#include <cfloat>
#include <cmath>

namespace aux {

Distance::Distance() { }

Distance::~Distance() { }

void Distance::dist_max1Dlog(const gpr::Coord& x1, const gpr::Coord& x2,
                             const gpr::AtomsConfiguration& conf_info,
                             gpr::Field<double>& dist)
{
    gpr::Index_t n1 = x1.getNumRows();
    gpr::Index_t n2 = x2.getNumRows();
    gpr::Index_t N_mov = x1.getNumCols() / 3;
    gpr::Index_t N_fro = conf_info.atoms_froz_active.positions.getNumCols() / 3;

    dist.resize(n1, n2);

    // Distances between moving atoms
    if (N_mov > 1) {
        for (gpr::Index_t i = 0; i < N_mov - 1; ++i) {
            for (gpr::Index_t j = i + 1; j < N_mov; ++j) {
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double r_ij_1 = (x1.at(n, j) - x1.at(n, i)).length();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double r_ij_2 = (x2.at(m, j) - x2.at(m, i)).length();
                        double tmp = fabs(log(r_ij_2 / r_ij_1));
                        if (tmp > dist(n, m)) dist(n, m) = tmp;
                    }
                }
            }
        }
    }

    // Distances from moving atoms to active frozen atoms
    if (N_fro > 0) {
        for (gpr::Index_t j = 0; j < N_mov; ++j) {
            for (gpr::Index_t i = 0; i < N_fro; ++i) {
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double r_ij_1 =
                        (x1.at(n, j) -
                         conf_info.atoms_froz_active.positions.at(0, i))
                            .length();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double r_ij_2 =
                            (x2.at(m, j) -
                             conf_info.atoms_froz_active.positions.at(0, i))
                                .length();
                        double tmp = fabs(log(r_ij_2 / r_ij_1));
                        if (tmp > dist(n, m)) dist(n, m) = tmp;
                    }
                }
            }
        }
    }
}

void Distance::dist_at(const gpr::Coord& x1, const gpr::Coord& x2,
                       const gpr::AtomsConfiguration& conf_info,
                       const gpr::Field<double>& lengthscale,
                       gpr::Field<double>& dist)
{
    gpr::Index_t n1 = x1.getNumRows();
    gpr::Index_t n2 = x2.getNumRows();
    gpr::Index_t N_mov = conf_info.atoms_mov.type.getSize();
    gpr::Index_t N_fro = conf_info.atoms_froz_active.type.getSize();
    gpr::Field<double> s2;
    gpr::Index_t size_s2 = lengthscale.getSize();

    dist.resize(n1, n2);

    s2.resize(1, size_s2);
    for (gpr::Index_t n = 0; n < s2.getSize(); ++n)
        s2[n] = 1. / (lengthscale[n] * lengthscale[n]);

    // If ARD is not used make s a vector of equal elements
    if (size_s2 == 1) {
        double ref_s2 = s2(0, 0);
        size_s2 = conf_info.n_pt;
        s2.resize(1, size_s2);
        for (gpr::Index_t n = 0; n < size_s2; ++n) {
            s2(0, n) = ref_s2;
        }
    }

    // Distances between moving atoms
    if (N_mov > 1) {
        for (gpr::Index_t j = 0; j < N_mov - 1; ++j) {
            for (gpr::Index_t i = j + 1; i < N_mov; ++i) {
                double s2_val =
                    s2(0, conf_info.pairtype(conf_info.atoms_mov.type(0, i),
                                             conf_info.atoms_mov.type(0, j)));
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double invr_ij_1 = (x1.at(n, j) - x1.at(n, i)).rlength();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double invr_ij_2 =
                            (x2.at(m, j) - x2.at(m, i)).rlength();
                        double invr_diff = invr_ij_1 - invr_ij_2;
                        dist(n, m) += 2. * s2_val * invr_diff * invr_diff;
                    }
                }
            }
        }
    }

    // Distances from moving atoms to active frozen atoms
    if (N_fro > 0) {
        for (gpr::Index_t j = 0; j < N_mov; ++j) {
            for (gpr::Index_t i = 0; i < N_fro; ++i) {
                double s2_val = s2(
                    0,
                    conf_info.pairtype(conf_info.atoms_froz_active.type(0, i),
                                       conf_info.atoms_mov.type(0, j)));
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double invr_ij_1 =
                        (x1.at(n, j) -
                         conf_info.atoms_froz_active.positions.at(0, i))
                            .rlength();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double invr_ij_2 =
                            (x2.at(m, j) -
                             conf_info.atoms_froz_active.positions.at(0, i))
                                .rlength();
                        double invr_diff = invr_ij_1 - invr_ij_2;
                        dist(n, m) += 2. * s2_val * invr_diff * invr_diff;
                    }
                }
            }
        }
    }

    for (gpr::Index_t n = 0; n < dist.getSize(); ++n) {
        dist[n] = sqrt(dist[n]);
    }
}

void Distance::dist_at_vec(const gpr::Coord& x1, const gpr::Coord& x2,
                           const gpr::AtomsConfiguration& conf_info,
                           const gpr::Field<double>& lengthscale,
                           std::vector<gpr::Field<double> >& dist)
{
    gpr::Index_t n1 = x1.getNumRows();
    gpr::Index_t n2 = x2.getNumRows();
    gpr::Index_t N_mov = conf_info.atoms_mov.type.getSize();
    gpr::Index_t N_fro = conf_info.atoms_froz_active.type.getSize();
    gpr::Field<double> s2;
    gpr::Index_t size_s2 = lengthscale.getSize();

    dist.resize(conf_info.n_pt);
    for (gpr::Index_t n = 0; n < dist.size(); ++n)
        dist[n].resize(n1, n2);

    s2.resize(1, size_s2);
    for (gpr::Index_t n = 0; n < size_s2; ++n)
        s2(0, n) = 1. / (lengthscale(0, n) * lengthscale(0, n));

    // If ARD is not used make s a vector of equal elements
    if (size_s2 == 1) {
        double ref_s2 = s2(0, 0);
        size_s2 = conf_info.n_pt;
        s2.resize(1, size_s2);
        for (gpr::Index_t n = 0; n < size_s2; ++n) {
            s2(0, n) = ref_s2;
        }
    }

    // Distances between moving atoms
    if (N_mov > 1) {
        for (gpr::Index_t j = 0; j < N_mov - 1; ++j) {
            for (gpr::Index_t i = j + 1; i < N_mov; ++i) {
                gpr::Index_t pt =
                    conf_info.pairtype(conf_info.atoms_mov.type(0, i),
                                       conf_info.atoms_mov.type(0, j));
                double s2_val = s2(0, pt);
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double invr_ij_1 = (x1.at(n, j) - x1.at(n, i)).rlength();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double invr_ij_2 =
                            (x2.at(m, j) - x2.at(m, i)).rlength();
                        double invr_diff = invr_ij_1 - invr_ij_2;
                        dist[pt](n, m) -= 2. * s2_val * invr_diff * invr_diff;
                    }
                }
            }
        }
    }

    // Distances from moving atoms to active frozen atoms
    if (N_fro > 0) {
        for (gpr::Index_t j = 0; j < N_mov; ++j) {
            for (gpr::Index_t i = 0; i < N_fro; ++i) {
                gpr::Index_t pt =
                    conf_info.pairtype(conf_info.atoms_froz_active.type(0, i),
                                       conf_info.atoms_mov.type(0, j));
                double s2_val = s2(0, pt);
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double invr_ij_1 =
                        (x1.at(n, j) -
                         conf_info.atoms_froz_active.positions.at(0, i))
                            .rlength();
                    for (gpr::Index_t m = 0; m < n2; ++m) {
                        double invr_ij_2 =
                            (x2.at(m, j) -
                             conf_info.atoms_froz_active.positions.at(0, i))
                                .rlength();
                        double invr_diff = invr_ij_1 - invr_ij_2;
                        dist[pt](n, m) -= 2. * s2_val * invr_diff * invr_diff;
                    }
                }
            }
        }
    }

    //    for(gpr::Index_t n = 0; n < n1; ++n) {
    //        for(gpr::Index_t m = 0; m < n2; ++m) {
    //            dist(n, m) = sqrt(dist(n, m));
    //        }
    //    }
}

void Distance::mindist_interatomic(const gpr::Coord& x,
                                   const gpr::AtomsConfiguration& conf_info,
                                   gpr::Field<double>& dist)
{
    gpr::Index_t n1 = x.getNumRows();
    gpr::Index_t N_mov = x.getSize() / 3;
    gpr::Index_t N_fro = conf_info.atoms_froz_active.positions.getSize() / 3;

    dist.resize(n1, N_mov);
    dist.set(DBL_MAX);

    // Distances between moving atoms
    if (N_mov > 1) {
        for (gpr::Index_t i = 0; i < N_mov - 1; ++i) {
            for (gpr::Index_t j = i + 1; j < N_mov; ++j) {
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double r_ij = (x.at(n, j) - x.at(n, i)).length();
                    for (gpr::Index_t m = 0; m < N_mov; ++m) {
                        if (r_ij < dist(n, i)) dist(n, i) = r_ij;
                        if (r_ij < dist(n, j)) dist(n, j) = r_ij;
                    }
                }
            }
        }
    }

    // Distances from moving atoms to active frozen atoms
    if (N_fro > 0) {
        for (gpr::Index_t j = 0; j < N_mov; ++j) {
            for (gpr::Index_t i = 0; i < N_fro; ++i) {
                for (gpr::Index_t n = 0; n < n1; ++n) {
                    double r_ij =
                        (x.at(n, j) -
                         conf_info.atoms_froz_active.positions.at(0, i))
                            .length();
                    for (gpr::Index_t m = 0; m < N_mov; ++m) {
                        if (r_ij < dist(n, m)) dist(n, m) = r_ij;
                    }
                }
            }
        }
    }
}

} /* namespace aux */
