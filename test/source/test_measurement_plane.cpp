#include <gtest/gtest.h>
#include <promille/promille.hpp>

#include "dummy_residual_model.hpp"

TEST(PlaneAndModels, MeasurementPlane)
{
    using the_type = double;
    auto gp1 = promille::global_parameter<the_type>(21, 0, "GP1");
    auto gp2 = promille::global_parameter<the_type>(22, 0, "GP2");
    auto gp3 = promille::global_parameter<the_type>(23, 0, "GP3");

    auto plane1 = promille::measurement_plane<the_type, promille_tests::dummy_residual_model<the_type, the_type>>(gp1, gp2, gp3);

    // plane1.print();

    EXPECT_EQ(&gp1, &plane1.globals_kind[0].get());
    EXPECT_EQ(&gp2, &plane1.globals_kind[1].get());
    EXPECT_EQ(&gp3, &plane1.globals_kind[2].get());

    plane1.set_globals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE);
    // plane1.print();

    EXPECT_TRUE(gp1.is_anywere_free);
    EXPECT_FALSE(gp2.is_anywere_free);

    plane1.set_globals_configuration({promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE});
    // plane1.print();

    plane1.set_globals_configuration(promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED);
    plane1.set_globals_configuration({promille::Kind::FIXED, promille::Kind::FIXED});

    plane1.set_locals_configuration(promille::Kind::FIXED, promille::Kind::FIXED);
    plane1.set_locals_configuration({promille::Kind::FIXED, promille::Kind::FIXED});
}

TEST(PlaneAndModels, ModelPlanes)
{
    using the_type = double;
    auto gp1 = promille::global_parameter<the_type>(21, 0, "GP1");
    auto gp2 = promille::global_parameter<the_type>(22, 0, "GP2");
    auto gp3 = promille::global_parameter<the_type>(23, 0, "GP3");

    auto model_plane1 = promille::model_planes<the_type, promille_tests::dummy_residual_model<the_type, the_type>>(nullptr, false);

    EXPECT_THROW(model_plane1.plane(10), std::runtime_error);
}
