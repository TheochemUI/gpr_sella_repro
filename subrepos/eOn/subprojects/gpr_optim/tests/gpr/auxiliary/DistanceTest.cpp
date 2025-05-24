/*
 * DistanceTest.cpp
 *
 *  Created on: 30 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "DistanceTest.h"

namespace gpr {
namespace tests {

DistanceTest::DistanceTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.

    conf_info.atoms_froz_active.positions.resize(1, 26 * 3);
    conf_info.atoms_froz_active.type.resize(1, 26);
    conf_info.atoms_mov.type.resize(1, 2);
    conf_info.pairtype.resize(2, 2);

    io_manager.readFromPlainFile(
        "tests/reference/gpr/auxiliary/ConfigurationFrozen.dat",
        conf_info.atoms_froz_active.positions.getInternalVector());

    conf_info.atoms_mov.type.set(0);
    conf_info.atoms_froz_active.type.set(1);

    conf_info.pairtype(0, 0) = 0;
    conf_info.pairtype(0, 1) = 1;
    conf_info.pairtype(1, 0) = 1;
    conf_info.pairtype(1, 1) = EMPTY;

    conf_info.n_pt = 2;
}

DistanceTest::~DistanceTest() { }

TEST_F(DistanceTest, DistMax1Dlog)
{
    gpr::Coord x1, x2;
    gpr::Field<double> dist;
    std::vector<double> dist_ref;

    x1.resize(1, 2 * 3);
    x2.resize(6, 2 * 3);
    dist_ref.resize(6);

    x1(0, 0) = 8.984799577839603;  // .x
    x1(0, 1) = 9.946773537381196;  // .y
    x1(0, 2) = 7.883158760864689;  // .z
    x1(0, 3) = 7.648053020353811;  // .x
    x1(0, 4) = 9.947018098632240;  // .y
    x1(0, 5) = 7.884341447664946;  // .z

    x2(0, 0) = 8.982373164830571;  // .x
    x2(0, 1) = 9.937230835772036;  // .y
    x2(0, 2) = 7.894416323850492;  // .z
    x2(0, 3) = 7.652483227274959;  // .x
    x2(0, 4) = 9.955905494573976;  // .y
    x2(0, 5) = 7.877879589983660;  // .z
                                   //
    x2(1, 0) = 8.978562773030582;  // .x
    x2(1, 1) = 9.932116280672286;  // .y
    x2(1, 2) = 7.898827614144257;  // .z
    x2(1, 3) = 7.648887496635559;  // .x
    x2(1, 4) = 9.955170515128863;  // .y
    x2(1, 5) = 7.872742150466696;  // .z
                                   //
    x2(2, 0) = 8.976993347159784;  // .x
    x2(2, 1) = 9.935874361177897;  // .y
    x2(2, 2) = 7.894533919329630;  // .z
    x2(2, 3) = 7.644176926015588;  // .x
    x2(2, 4) = 9.955573130196294;  // .y
    x2(2, 5) = 7.878193601315998;  // .z
                                   //
    x2(3, 0) = 8.982700234844764;  // .x
    x2(3, 1) = 9.938213927286730;  // .y
    x2(3, 2) = 7.892109355815186;  // .z
    x2(3, 3) = 7.643018372071518;  // .x
    x2(3, 4) = 9.956281713287362;  // .y
    x2(3, 5) = 7.875909963059283;  // .z
                                   //
    x2(4, 0) = 8.989491283507103;  // .x
    x2(4, 1) = 9.936555712126957;  // .y
    x2(4, 2) = 7.892532513641982;  // .z
    x2(4, 3) = 7.645972793452610;  // .x
    x2(4, 4) = 9.955825824010615;  // .y
    x2(4, 5) = 7.876166184685688;  // .z
                                   //
    x2(5, 0) = 8.983330780846810;  // .x
    x2(5, 1) = 9.944154132051644;  // .y
    x2(5, 2) = 7.875729042237787;  // .z
    x2(5, 3) = 7.650585914617185;  // .x
    x2(5, 4) = 9.929441987442432;  // .y
    x2(5, 5) = 7.882144032513800;  // .z

    dist_ref = {0.005675600658215, 0.008024540882219, 0.005889796208082,
                0.004873773574595, 0.006608214559330, 0.007841185655827};

    distance.dist_max1Dlog(x1, x2, conf_info, dist);

    for (gpr::Index_t n = 0; n < dist.getSize(); ++n) {
        EXPECT_LE(std::fabs(dist(0, n) - dist_ref[n]), threshold)
            << "Elements of the matrix are not equal to the expected ones.";
    }
}

TEST_F(DistanceTest, DistAt)
{
    gpr::Coord x1;
    gpr::Coord x2;
    gpr::Field<double> lengthscale;
    gpr::Field<double> dist;
    gpr::Field<double> dist_ref(6, 6);

    x1.resize(6, 2 * 3);
    x2.resize(6, 2 * 3);
    dist.resize(5, 5);
    lengthscale.resize(1, 1);

    x1.set(0, 0, {8.982373164830571, 9.937230835772036, 7.894416323850492});
    x1.set(1, 0, {8.978562773030582, 9.932116280672286, 7.898827614144257});
    x1.set(2, 0, {8.976993347159784, 9.935874361177897, 7.894533919329630});
    x1.set(3, 0, {8.982700234844764, 9.938213927286730, 7.892109355815186});
    x1.set(4, 0, {8.989491283507103, 9.936555712126957, 7.892532513641982});
    x1.set(5, 0, {8.983330780846810, 9.944154132051644, 7.875729042237787});

    x1.set(0, 1, {7.652483227274959, 9.955905494573976, 7.877879589983660});
    x1.set(1, 1, {7.648887496635559, 9.955170515128863, 7.872742150466696});
    x1.set(2, 1, {7.644176926015588, 9.955573130196294, 7.878193601315998});
    x1.set(3, 1, {7.643018372071518, 9.956281713287362, 7.875909963059283});
    x1.set(4, 1, {7.645972793452610, 9.955825824010615, 7.876166184685688});
    x1.set(5, 1, {7.650585914617185, 9.929441987442432, 7.882144032513800});

    x2.set(0, 0, {8.982373164830571, 9.937230835772036, 7.894416323850492});
    x2.set(1, 0, {8.978562773030582, 9.932116280672286, 7.898827614144257});
    x2.set(2, 0, {8.976993347159784, 9.935874361177897, 7.894533919329630});
    x2.set(3, 0, {8.982700234844764, 9.938213927286730, 7.892109355815186});
    x2.set(4, 0, {8.989491283507103, 9.936555712126957, 7.892532513641982});
    x2.set(5, 0, {8.983330780846810, 9.944154132051644, 7.875729042237787});

    x2.set(0, 1, {7.652483227274959, 9.955905494573976, 7.877879589983660});
    x2.set(1, 1, {7.648887496635559, 9.955170515128863, 7.872742150466696});
    x2.set(2, 1, {7.644176926015588, 9.955573130196294, 7.878193601315998});
    x2.set(3, 1, {7.643018372071518, 9.956281713287362, 7.875909963059283});
    x2.set(4, 1, {7.645972793452610, 9.955825824010615, 7.876166184685688});
    x2.set(5, 1, {7.650585914617185, 9.929441987442432, 7.882144032513800});

    lengthscale(0, 0) = 1.;

    dist_ref(0, 0) = 0.;
    dist_ref(0, 1) = 0.003953892010280;
    dist_ref(0, 2) = 0.004070202425876;
    dist_ref(0, 3) = 0.008460855559739;
    dist_ref(0, 4) = 0.011285766759530;
    dist_ref(0, 5) = 0.014672259642108;
    dist_ref(1, 0) = 0.003953892010280;
    dist_ref(1, 1) = 0.;
    dist_ref(1, 2) = 0.004423202945557;
    dist_ref(1, 3) = 0.009080320379056;
    dist_ref(1, 4) = 0.011918448534057;
    dist_ref(1, 5) = 0.016476334213804;
    dist_ref(2, 0) = 0.004070202425876;
    dist_ref(2, 1) = 0.004423202945557;
    dist_ref(2, 2) = 0.;
    dist_ref(2, 3) = 0.005989481708574;
    dist_ref(2, 4) = 0.009472766608845;
    dist_ref(2, 5) = 0.014771352424969;
    dist_ref(3, 0) = 0.008460855559739;
    dist_ref(3, 1) = 0.009080320379056;
    dist_ref(3, 2) = 0.005989481708574;
    dist_ref(3, 3) = 0.;
    dist_ref(3, 4) = 0.003964897696479;
    dist_ref(3, 5) = 0.015433995102774;
    dist_ref(4, 0) = 0.011285766759530;
    dist_ref(4, 1) = 0.011918448534057;
    dist_ref(4, 2) = 0.009472766608845;
    dist_ref(4, 3) = 0.003964897696479;
    dist_ref(4, 4) = 0.;
    dist_ref(4, 5) = 0.016861382029173;
    dist_ref(5, 0) = 0.014672259642108;
    dist_ref(5, 1) = 0.016476334213804;
    dist_ref(5, 2) = 0.014771352424969;
    dist_ref(5, 3) = 0.015433995102774;
    dist_ref(5, 4) = 0.016861382029173;
    dist_ref(5, 5) = 0.;

    distance.dist_at(x1, x2, conf_info, lengthscale, dist);

    for (gpr::Index_t n = 0; n < dist.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < dist.getNumCols(); ++m) {
            EXPECT_LE(std::fabs(dist(n, m) - dist_ref(n, m)), threshold)
                << "Elements of the matrix are not equal to the expected ones.";
        }
    }
}

TEST_F(DistanceTest, MinDistInteratomic)
{
    gpr::Coord x;
    gpr::Field<double> dist;
    gpr::Field<double> dist_ref(1, 2);

    x.resize(1, 3 * 2);

    x.set(0, 0, {8.982538896311098, 9.937443825779654, 7.894092264564180});
    x.set(0, 1, {7.652281581062163, 9.955691736800297, 7.877995340784432});

    dist_ref(0, 0) = 1.330479846515945;
    dist_ref(0, 1) = 1.330479846515945;

    distance.mindist_interatomic(x, conf_info, dist);

    for (gpr::Index_t n = 0; n < dist.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < dist.getNumCols(); ++m) {
            EXPECT_LE(std::fabs(dist(n, m) - dist_ref(n, m)), threshold)
                << "Elements of the matrix are not equal to the expected ones.";
        }
    }
}

} /* namespace tests */
} /* namespace gpr */
