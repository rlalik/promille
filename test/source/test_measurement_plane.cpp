#include <gtest/gtest.h>
#include <promille/promille.hpp>

#include "dummy_residual_model.hpp"

TEST(MeasurementPlanes, MeasurementCreation)
{
    using the_type = double;
    auto gp1 = promille::global_parameter<the_type>(21, 0, "GP1");
    auto gp2 = promille::global_parameter<the_type>(22, 0, "GP2");
    auto gp3 = promille::global_parameter<the_type>(23, 0, "GP3");

    auto plane1 = promille::alignment_plane<the_type, 2>(
        gp1.make_state_from_parameter(), gp2.make_state_from_parameter(), gp3.make_state_from_parameter());

    EXPECT_EQ(&gp1, &plane1.globals[0].get());
    EXPECT_EQ(&gp2, &plane1.globals[1].get());
    EXPECT_EQ(&gp3, &plane1.globals[2].get());

    plane1.set_globals_configuration({promille::Kind::FREE, promille::Kind::FREE, promille::Kind::FREE});
    plane1.set_globals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE);

    EXPECT_TRUE(gp1.is_anywere_free);
    EXPECT_FALSE(gp2.is_anywere_free);

    auto model1 = promille_tests::dummy_residual_model<the_type, the_type>(0, 1, 2);

    auto measurement1 = promille::measurement_plane<the_type, decltype(plane1), 2>(plane1, model1);

    measurement1.get_alignment().set_globals_configuration(promille::Kind::FIXED, promille::Kind::FIXED, promille::Kind::FIXED);
    measurement1.get_alignment().set_globals_configuration({promille::Kind::FIXED, promille::Kind::FIXED});

    measurement1.get_alignment().set_locals_configuration(promille::Kind::FIXED, promille::Kind::FIXED);
    measurement1.get_alignment().set_locals_configuration({promille::Kind::FIXED, promille::Kind::FIXED});
}
