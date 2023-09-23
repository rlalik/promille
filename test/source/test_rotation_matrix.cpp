#include <StrawAlignment/StrawAlignment.hpp>
#include <gtest/gtest.h>

TEST(RotationMatrix, Components)
{
    float phi = 0.0;
    float theta = 0.0;
    float ypsilon = 0.0;

    EXPECT_DOUBLE_EQ(SA::helper::R11(phi, theta, ypsilon), 1.0);
    EXPECT_DOUBLE_EQ(SA::helper::R12(phi, theta, ypsilon), 0.0);
    EXPECT_DOUBLE_EQ(SA::helper::R13(phi, theta, ypsilon), 0.0);

    EXPECT_DOUBLE_EQ(SA::helper::R21(phi, theta, ypsilon), 0.0);
    EXPECT_DOUBLE_EQ(SA::helper::R22(phi, theta, ypsilon), 1.0);
    EXPECT_DOUBLE_EQ(SA::helper::R23(phi, theta, ypsilon), 0.0);

    EXPECT_DOUBLE_EQ(SA::helper::R31(phi, theta, ypsilon), 0.0);
    EXPECT_DOUBLE_EQ(SA::helper::R32(phi, theta, ypsilon), 0.0);
    EXPECT_DOUBLE_EQ(SA::helper::R33(phi, theta, ypsilon), 1.0);
}

TEST(RotationMatrix, Matrix)
{
    float phi = 0.0;
    float theta = 0.0;
    float ypsilon = 0.0;

    auto R = SA::make_rotation_matrix(phi, theta, ypsilon);

    Rotation3D Rt(1, 0, 0, 0, 1, 0, 0, 0, 1);
    EXPECT_EQ(R, Rt);
}

TEST(RotationMatrix, Rotating)
{
    double phi = 0.0;
    double theta = 0.0;
    double ypsilon = 0.0;

    auto R = SA::make_rotation_matrix(phi, theta, ypsilon);

    XYZVector v(1.0, 0.0, 0.0);

    auto res = R * v;

    EXPECT_EQ(res, v);
}
