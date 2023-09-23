#include <StrawAlignment/StrawAlignment.hpp>
#include <gtest/gtest.h>

TEST(Mille, AddLocals)
{
    auto mille = SA::MilleBuilder("test");

    mille.add_planes_global(0, 0, 0, 0, 0, 0);

    mille.add_local(0, 0, 0, 0, 0, 0);
}
