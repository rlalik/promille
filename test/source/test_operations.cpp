#include <StrawAlignment/StrawAlignment.hpp>
#include <gtest/gtest.h>

TEST(Operations, Components)
{
    {
        auto b1 = SA::make_point(0, 0, 0);
        auto d1 = SA::make_vector(0, 0, 1);

        auto b2 = SA::make_point(1, 0, 0);
        auto d2 = SA::make_vector(0, 1, 0);

        auto dist = SA::distance(b1, d1, b2, d2);
        EXPECT_FLOAT_EQ(dist, 1);
    }

    {
        auto b1 = SA::make_point(0, 0, 0);
        auto d1 = SA::make_vector(0, 0, 1);

        auto b2 = SA::make_point(2, 0, 0);
        auto d2 = SA::make_vector(0, 1, 0);

        auto dist = SA::distance(b1, d1, b2, d2);
        EXPECT_FLOAT_EQ(dist, 2);
    }
}
