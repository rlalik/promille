#include <gtest/gtest.h>
#include <mille_builder/mille_builder.hpp>

TEST(Operations, Components)
{
    {
        auto b1 = mb::XYZPoint(0, 0, 0);
        auto d1 = mb::XYZVector(0, 0, 1);

        auto b2 = mb::XYZPoint(1, 0, 0);
        auto d2 = mb::XYZVector(0, 1, 0);

        auto dist = mb::geom::distance(b1, d1, b2, d2);
        EXPECT_FLOAT_EQ(dist, 1);
    }

    {
        auto b1 = mb::XYZPoint(0, 0, 0);
        auto d1 = mb::XYZVector(0, 0, 1);

        auto b2 = mb::XYZPoint(2, 0, 0);
        auto d2 = mb::XYZVector(0, 1, 0);

        auto dist = mb::geom::distance(b1, d1, b2, d2);
        EXPECT_FLOAT_EQ(dist, 2);
    }
}
