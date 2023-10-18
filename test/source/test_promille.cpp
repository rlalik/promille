#include <gtest/gtest.h>
#include <promille/promille.hpp>

#include "dummy_residual_model.hpp"

using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

TEST(Mille, SingleModel)
{
    promille::promille mille("test_", "test.bin");

    mille.add_global_parameter(1, -10, "P1");
    mille.add_global_parameter(2, -20, "P2");
    mille.add_global_parameter(3, -30, "P3");

    auto straw_planes = mille.make_model_planes<promille_tests::dummy_residual_model<float, float>>();

    straw_planes.add_plane(1, 1, 2, 3)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED)
        .set_globals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE);

    straw_planes.plane(1).add_measurement(0);

    mille.end();
}

TEST(Mille, DoubleModel)
{
    promille::promille mille("test_", "test.bin");

    mille.set_verbose(2);

    mille.add_global_parameter(1, -10, "P1");
    mille.add_global_parameter(2, -20, "P2");
    mille.add_global_parameter(3, -30, "P3");

    auto straw_planes1 = mille.make_model_planes<promille_tests::dummy_residual_model<float, float>>();

    auto straw_planes2 = mille.make_model_planes<promille_tests::dummy_residual_model_2<float, float>>();

    straw_planes1.add_plane(1, 1, 2, 3)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED)
        .set_globals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE);

    straw_planes1.plane(1).add_measurement(0);

    straw_planes2.add_plane(2, 1, 2, 3, mille.add_global_parameter(4, -40, "P4"))
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED)
        .set_globals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FREE);

    // straw_planes2.plane(2).add_measurement(0);

    mille.end();
}
