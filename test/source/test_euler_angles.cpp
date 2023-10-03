#include <gtest/gtest.h>
#include <mille_builder/euler_angles.hpp>
#include <mille_builder/mille_builder.hpp>

static const double comp_error = 0.00001;
static const double comp_error_det = 0.00001;
TEST(EulerAngles, ZYZ_angles)
{
    std::vector<std::pair<std::array<double, 3>, std::pair<std::array<double, 3>, std::array<double, 3>>>> data = {
        // clang-format off
        {{0, 0, 0}, {{0, 0, 0}, {1, 2, 3}}},
        {{1, 1, 1}, {{1, 1, 1}, {-0.586992, 2.64322, 2.58241}}}
        // clang-format on
    };

    const auto rot_I = ROOT::Math::Rotation3D(1, 0, 0, 0, 1, 0, 0, 0, 1);
    std::array<ROOT::Math::Rotation3D::Scalar, 9> components_I;
    rot_I.GetComponents(components_I.begin());

    for (const auto& d : data) {
        const mb::euler::zyz<double> s_zyz(d.first[0], d.first[1], d.first[2]);

        EXPECT_NEAR(s_zyz.a1, d.second.first[0], comp_error);
        EXPECT_NEAR(s_zyz.a2, d.second.first[1], comp_error);
        EXPECT_NEAR(s_zyz.a3, d.second.first[2], comp_error);

        EXPECT_NEAR(s_zyz.determinat(), 1.0, comp_error);

        const auto rot = mb::euler::make_rotation_matrix(s_zyz);
        const auto rot_t = rot.Inverse();

        auto rrT = rot * rot_t;
        decltype(components_I) components_rrT;
        rrT.GetComponents(components_rrT.begin());

        for (size_t i = 0; i < 9; ++i)
            EXPECT_NEAR(components_rrT[i], components_I[i], comp_error_det);

        const mb::XYZVector track(1, 2, 3);
        const auto res = rot * track;

        EXPECT_NEAR(res.x(), d.second.second[0], comp_error);
        EXPECT_NEAR(res.y(), d.second.second[1], comp_error);
        EXPECT_NEAR(res.z(), d.second.second[2], comp_error);
    }
}

TEST(EulerAngles, ZYZ_matrix)
{
    std::vector<std::pair<std::array<double, 9>, std::pair<std::array<double, 3>, std::array<double, 3>>>> data = {
        // clang-format off
        {{1, 0, 0, 0, 1, 0, 0, 0, 1}, {{0, 0, 0}, {1, 2, 3}}},
        {   // 10, 20, 30 deg
            {0.71461, -0.613092, 0.336824, 0.633718, 0.771281, 0.0593912, -0.296198, 0.17101, 0.939693},
            {{0.17453293, 0.34906585, 0.52359878}, {0.4988980, 2.35445, 2.86490}}
        }
        // clang-format on
    };

    const auto rot_I = ROOT::Math::Rotation3D(1, 0, 0, 0, 1, 0, 0, 0, 1);
    std::array<ROOT::Math::Rotation3D::Scalar, 9> components_I;
    rot_I.GetComponents(components_I.begin());

    for (const auto& d : data) {
        const mb::euler::zyz<double> s_zyz(
            d.first[0], d.first[1], d.first[2], d.first[3], d.first[4], d.first[5], d.first[6], d.first[7], d.first[8]);

        EXPECT_NEAR(s_zyz.a1, d.second.first[0], comp_error);
        EXPECT_NEAR(s_zyz.a2, d.second.first[1], comp_error);
        EXPECT_NEAR(s_zyz.a3, d.second.first[2], comp_error);

        EXPECT_NEAR(s_zyz.determinat(), 1.0, comp_error);

        const auto rot = mb::euler::make_rotation_matrix(s_zyz);
        const auto rot_t = rot.Inverse();

        auto rrT = rot * rot_t;
        decltype(components_I) components_rrT;
        rrT.GetComponents(components_rrT.begin());

        for (size_t i = 0; i < 9; ++i)
            EXPECT_NEAR(components_rrT[i], components_I[i], comp_error_det);

        const mb::XYZVector track(1, 2, 3);
        const auto res = rot * track;

        EXPECT_NEAR(res.x(), d.second.second[0], comp_error);
        EXPECT_NEAR(res.y(), d.second.second[1], comp_error);
        EXPECT_NEAR(res.z(), d.second.second[2], comp_error);
    }
}
