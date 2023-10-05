#include <gtest/gtest.h>
#include <mille_builder/alignment_model.hpp>
#include <mille_builder/mille_builder.hpp>

using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

TEST(Mille, AddLocals)
{
    mb::mille_builder<mb::sts_residua<float, mb::euler::zyz>> mille("test_", "test");
    mille.set_verbose(2);

    mille.add_global_parameter(1, -10, "P1");
    mille.add_global_parameter(2, -20, "P2");
    mille.add_global_parameter(3, -30, "P3");
    mille.add_global_parameter(4, -40, "P4");
    mille.add_global_parameter(5, -50, "P5");
    mille.add_global_parameter(6, -60, "P6");

    // mille.add_globals_set(0,
    //                       {1, 0, mb::Kind::FIXED},
    //                       {2, 0, mb::Kind::FIXED},
    //                       {3, 0, mb::Kind::FIXED},
    //                       {4, 0, mb::Kind::FIXED},
    //                       {5, 0, mb::Kind::FIXED},
    //                       {6, 0, mb::Kind::FIXED},
    //                       0,
    //                       0,
    //                       0,
    //                       0,
    //                       0,
    //                       0,
    //                       mb::Kind::FIXED,
    //                       mb::Kind::FIXED,
    //                       mb::Kind::FIXED,
    //                       mb::Kind::FIXED);

    mille.add_plane(0, 1, 2, 3, 4, 5, 6);

    mille.add_measurement(0, 0.1f, 0.15f, XYZPoint(0, 0, 0), XYZVector(0, 0, 1), XYZPoint(0, 0, 0), XYZVector(0, 0, 1));

    mille.end();
}
