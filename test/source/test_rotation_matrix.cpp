#include <StrawAlignment/StrawAlignment.hpp>
#include <gtest/gtest.h>

TEST(RotationMatrix, Components)
{
    float psi = 0.0;
    float theta = 0.0;
    float phi = 0.0;

    EXPECT_DOUBLE_EQ(SA::helper::R11(psi, theta, phi), 1.0);
    EXPECT_DOUBLE_EQ(SA::helper::R12(psi, theta, phi), 0.0);
    EXPECT_DOUBLE_EQ(SA::helper::R13(psi, theta, phi), 0.0);

    EXPECT_DOUBLE_EQ(SA::helper::R21(psi, theta, phi), 0.0);
    EXPECT_DOUBLE_EQ(SA::helper::R22(psi, theta, phi), 1.0);
    EXPECT_DOUBLE_EQ(SA::helper::R23(psi, theta, phi), 0.0);

    EXPECT_DOUBLE_EQ(SA::helper::R31(psi, theta, phi), 0.0);
    EXPECT_DOUBLE_EQ(SA::helper::R32(psi, theta, phi), 0.0);
    EXPECT_DOUBLE_EQ(SA::helper::R33(psi, theta, phi), 1.0);
}

TEST(RotationMatrix, Matrix)
{
    float psi = 0.0;
    float theta = 0.0;
    float phi = 0.0;

    auto R = SA::make_rotation_matrix(psi, theta, phi);

    Rotation3D Rt(1, 0, 0, 0, 1, 0, 0, 0, 1);
    EXPECT_EQ(R, Rt);
}

TEST(RotationMatrix, Rotating)
{
    double psi = 0.0;
    double theta = 0.0;
    double phi = 0.0;

    auto R = SA::make_rotation_matrix(psi, theta, phi);

    XYZVector v(1.0, 0.0, 0.0);

    auto res = R * v;

    EXPECT_EQ(res, v);
}
