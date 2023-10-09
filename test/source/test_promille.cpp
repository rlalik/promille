#include <gtest/gtest.h>
#include <promille/promille.hpp>

#include "dummy_residual_model.hpp"

using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

TEST(Mille, AddLocals)
{
    promille::promille<mb_tests::dummy_residual_model<float, float>> mille("test_", "test");
    mille.set_verbose(2);

    mille.add_global_parameter(1, -10, "P1");
    mille.add_global_parameter(2, -20, "P2");
    mille.add_global_parameter(3, -30, "P3");
    mille.add_global_parameter(4, -40, "P4");
    mille.add_global_parameter(5, -50, "P5");
    mille.add_global_parameter(6, -60, "P6");

    mille.add_plane(0, 1, 2, 3, 4, 5, 6)
        .set_locals_configuration(promille::Kind::FREE, promille::Kind::FIXED, promille::Kind::FREE, promille::Kind::FREE)
        .set_globals_configuration();

    mille.print(true);

    mille.add_measurement(0, 0, 0, 0, 0);

    mille.end();
}
