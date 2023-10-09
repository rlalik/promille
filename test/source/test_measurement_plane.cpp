#include <gtest/gtest.h>
#include <promille/promille.hpp>

#include "dummy_alignment_model.hpp"

static double accuracy = 0.00001;

using ROOT::Math::Rotation3D;
using ROOT::Math::RotationZYX;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

TEST(MeasurementPlanes, StsDerivatives)
{
    std::vector<std::pair<std::pair<std::array<double, 12>, std::array<double, 7>>, std::array<double, 10>>> data = {
        {{{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0}}, {1, 0, 0, 0, 0, 0, -1, 0, 0, 0}},
        {{{0, 1, 2, 3, 4, 5, 1, 2, 3, 0, 0, 0}, {6, 7, 8, 9, 10, 11, 12}},
         {-0.436231, 0.330121, 0.837092, -4.71041, 1.5327, 1.85027, 0.436231, -0.330121, -0.316471, 0.239491}}};

    for (const auto& d : data) {
        auto derivs = mb_tests::dummy_residual_model<double, promille::euler::zyz>(
            d.first.first[0], d.first.first[1], d.first.first[2], d.first.first[3], d.first.first[4], d.first.first[5]);

        derivs.set_params(d.first.first[6], d.first.first[7], d.first.first[8], d.first.first[9], d.first.first[10], d.first.first[11]);

        derivs.residual(XYZPoint(d.first.second[3], d.first.second[4], 0),
                        XYZVector(d.first.second[5], d.first.second[6], 1),
                        XYZPoint(d.first.second[0], d.first.second[1], d.first.second[2]),
                        XYZVector(0, 0, 1));

        EXPECT_NEAR(derivs.global_derivative(1), d.second[0], accuracy);
        EXPECT_NEAR(derivs.global_derivative(2), d.second[1], accuracy);
        EXPECT_NEAR(derivs.global_derivative(3), d.second[2], accuracy);

        EXPECT_NEAR(derivs.global_derivative(4), d.second[3], accuracy);
        EXPECT_NEAR(derivs.global_derivative(5), d.second[4], accuracy);
        EXPECT_NEAR(derivs.global_derivative(6), d.second[5], accuracy);

        EXPECT_NEAR(derivs.local_derivative(1), d.second[6], accuracy);
        EXPECT_NEAR(derivs.local_derivative(2), d.second[7], accuracy);

        EXPECT_NEAR(derivs.local_derivative(3), d.second[8], accuracy);
        EXPECT_NEAR(derivs.local_derivative(4), d.second[9], accuracy);
    }
}
