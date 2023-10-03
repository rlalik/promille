#include <gtest/gtest.h>
#include <mille_builder/euler_angles.hpp>
#include <mille_builder/mille_builder.hpp>

TEST(EulerAngles, ZYZ)
{
    float psi = 1.0;
    float theta = 1.0;
    float phi = 1.0;

    mb::euler::zyz<double> s_zyz(1, 1, 1);

    auto R = mb::euler::make_rotation_matrix(s_zyz);

    mb::XYZVector v(1, 2, 3);
    auto res = R * v;
    auto exp = mb::XYZVector(-0.586992, 2.64322, 2.58241);

    EXPECT_NEAR(res.x(), exp.x(), 0.0001);
    EXPECT_NEAR(res.y(), exp.y(), 0.0001);
    EXPECT_NEAR(res.z(), exp.z(), 0.0001);
}
