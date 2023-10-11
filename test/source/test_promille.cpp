#include <gtest/gtest.h>
#include <promille/promille.hpp>

#include "dummy_residual_model.hpp"

using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

TEST(Mille, CreateAndConfigure)
{
    promille::promille<promille_tests::dummy_residual_model<float, float>> mille("test_", "test");

    mille.add_global_parameter(1, -10, "P1");
    mille.add_global_parameter(2, -20, "P2");
    mille.add_global_parameter(3, -30, "P3");

    mille.add_plane(1, 1, 2, 3)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED)
        .set_globals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE);

    mille.add_measurement(1, 0, 0, 0);

    mille.end();
}
