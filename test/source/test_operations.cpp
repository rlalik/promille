#include <gtest/gtest.h>
#include <mille_builder/mille_builder.hpp>

TEST(Operations, Components)
{
    {
        const auto b1 = mb::XYZPoint(0, 0, 0);
        const auto d1 = mb::XYZVector(0, 0, 1);

        const auto b2 = mb::XYZPoint(1, 0, 0);
        const auto d2 = mb::XYZVector(0, 1, 0);

        const auto dist = mb::geom::distance(b1, d1, b2, d2);
        EXPECT_FLOAT_EQ(dist, 1);
    }

    {
        const auto b1 = mb::XYZPoint(0, 0, 0);
        const auto d1 = mb::XYZVector(0, 0, 1);

        const auto b2 = mb::XYZPoint(2, 0, 0);
        const auto d2 = mb::XYZVector(0, 1, 0);

        const auto dist = mb::geom::distance(b1, d1, b2, d2);
        EXPECT_FLOAT_EQ(dist, 2);
    }
}
