/*
 * Distance.h
 *
 *  Created on: 30 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_DISTANCE_H_
#define GPR_DISTANCE_H_

#include <vector>

#include "../../data_types/Coord.h"
#include "../../data_types/Field.h"
#include "../../structures/Structures.h"

namespace aux {

/**
 * @brief Methods to calculate distances between atoms.
 */
class Distance {
public:
    Distance();
    virtual ~Distance();

    /**
     * @brief Give the maximum difference in the logarithmic inter-atomic
     *        distances between atomic configurations C and C'.
     *
     * dist(C,C') = MAX_ij{|log(r_ij')-log(r_ij)|} = MAX_ij{|log(r_ij'/r_ij)|},
     * where r_ij and r_ij' are the distances between atoms i and j in
     * configurations C and C', respectively.
     *
     * The input vectors x1 and x2 are assumed to be row vectors including the
     * gpr::Coordinates of the moving atoms: [x_1,y_1,z_1,x_2,y_2,z_2,...].
     *
     * @param x1 gpr::Coordinates of the moving atoms in configurations C (n1 x
     * 3*N_mov)
     * @param x2 gpr::Coordinates of the moving atoms in configurations C' (n2 x
     * 3*N_mov)
     * @param conf_info Structure array including information about the
     * configurations necessary for the GP model
     *                  - conf_info.conf_fro: gpr::Coordinates of active frozen
     * atoms (N_fro x 3)
     *                  - conf_info.atomtype_mov: atomtype indices for moving
     * atoms (1 x N_mov)
     *                  - conf_info.atomtype_fro: pairtype indices for active
     * frozen atoms (1 x N_fro)
     *                  - conf_info.pairtype: pairtype indices for pairs of
     * atomtypes (n_at x n_at)
     *                  - conf_info.n_pt: number of active pairtypes
     *
     * @param dist Matrix including the "distances" from configurations C to
     * configurations C' (n1 x n2)
     */
    void dist_max1Dlog(const gpr::Coord& x1, const gpr::Coord& x2,
                       const gpr::AtomsConfiguration& conf_info,
                       gpr::Field<double>& dist);

    /**
     * @brief Give the distance between two atomic configurations as defined in
     *        the special GPstuff covariance function 'gpcf_matern32at'.
     *
     * The distance between configurations C and C' is based on the changes of
     * the inter-atomic distances:
     *
     * dist(C,C') = sqrt(SUM_ij{[(1/r_ij-1/r_ij')/l_ij]^2}), where r_ij and
     * r_ij' are the distances between atoms i and j in configurations C and
     * C', respectively, and l_ij is the lengthscale of the corresponding
     * atom pair type.
     *
     * The input vectors x1 and x2 are assumed to be row vectors including the
     * gpr::Coordinates of the moving atoms: [x_1,y_1,z_1,x_2,y_2,z_2,...].
     *
     * The parameter 'conf_info' is a structure array including necessary
     * information about the configurations: conf_info.conf_fro:
     * gpr::Coordinates of active frozen atoms (N_fro x 3)
     * conf_info.atomtype_mov: atomtype indices for moving atoms (1 x N_mov)
     * conf_info.atomtype_fro: atomtype indices for active frozen atoms (1 x
     * N_fro) Atomtypes must be indexed as 1,2,...,n_at (may include also
     * inactive atomtypes). conf_info.pairtype: pairtype indices for pairs of
     * atomtypes (n_at x n_at) conf_info.n_pt: number of active pairtypes Active
     * pairtypes are indexed as 0,1,2,...,n_pt-1. Inactive pairtypes are given
     * index -1.
     *
     * @param x1           gpr::Coordinates of the moving atoms in
     * configurations C (n1 x 3*N_mov)
     * @param x2           gpr::Coordinates of the moving atoms in
     * configurations C' (n2 x 3*N_mov)
     * @param conf_info    structure array including information about the
     * configurations necessary for the GP model
     * @param lengthscale  lengthscales for different atom pair types (1 x n_pt)
     * @param dist         matrix including the distances between configurations
     * C and C' (n1 x n2)
     */
    void dist_at(const gpr::Coord& x1, const gpr::Coord& x2,
                 const gpr::AtomsConfiguration& conf_info,
                 const gpr::Field<double>& lengthscale,
                 gpr::Field<double>& dist);

    /**
     * @brief Similar to \e dist_at but calculates vector of distances for every
     *        pair of moving and frozen atoms.
     *
     * @param x1
     * @param x2
     * @param conf_info
     * @param lengthscale
     * @param dist
     */
    void dist_at_vec(const gpr::Coord& x1, const gpr::Coord& x2,
                     const gpr::AtomsConfiguration& conf_info,
                     const gpr::Field<double>& lengthscale,
                     std::vector<gpr::Field<double> >& dist);

    /**
     * @brief Give the distance from each moving atom to its nearest neighbour
     *        atom in configuration C.
     *
     * The input vectors in x are assumed to be row vectors including the
     * gpr::Coordinates of the moving atoms: [x_1,y_1,z_1,x_2,y_2,z_2,...].
     *
     * @param x            gpr::Coordinates of the moving atoms in
     * configurations C (n x 3*N_mov)
     * @param conf_info    structure array including information about the
     * configurations necessary for the GP model
     *                      - conf_info.conf_fro: gpr::Coordinates of active
     * frozen atoms (N_fro x 3)
     *                      - conf_info.atomtype_mov: atomtype indices for
     * moving atoms (1 x N_mov)
     *                      - conf_info.atomtype_fro: pairtype indices for
     * active frozen atoms (1 x N_fro)
     *                      - conf_info.pairtype: pairtype indices for pairs of
     * atomtypes (n_at x n_at)
     *                      - conf_info.n_pt: number of active pairtypes
     *
     * @param dist         matrix including the minimum interatomic distances
     * for each moving atom in configurations C (n x N_mov)
     */
    void mindist_interatomic(const gpr::Coord& x,
                             const gpr::AtomsConfiguration& conf_info,
                             gpr::Field<double>& dist);
};

} /* namespace aux */

#endif /* GPR_DISTANCE_H_ */
