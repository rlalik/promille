#include <StrawAlignment/StrawAlignment.hpp>
#include <gtest/gtest.h>

TEST(Mille, AddLocals)
{
    auto mille = SA::MilleBuilder<SA::euler::zyz>("test_", "test");

    mille.add_planes_globals({0, SA::Kind::FIXED},
                             {0, SA::Kind::FIXED},
                             {0, SA::Kind::FIXED},
                             {0, SA::Kind::FIXED},
                             {0, SA::Kind::FIXED},
                             {0, SA::Kind::FIXED},
                             0,
                             0,
                             0,
                             0,
                             0,
                             0);

    mille.add_local(0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.15);

    mille.end();
}
