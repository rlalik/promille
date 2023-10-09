#include <gtest/gtest.h>
#include <promille/promille.hpp>

TEST(ParameterStates, GlobalParameters)
{
    auto param = promille::global_parameter<double>(1, 10.0);

    auto state1 = param.make_state_from_parameter();
    EXPECT_EQ(state1.get_kind(), promille::Kind::FIXED);

    auto state2 = param.make_state_from_parameter(promille::Kind::FIXED);
    EXPECT_EQ(state2.get_kind(), promille::Kind::FIXED);

    auto state3 = param.make_state_from_parameter(promille::Kind::FREE);
    EXPECT_EQ(state3.get_kind(), promille::Kind::FREE);

    state3.set_kind(promille::Kind::FIXED);
    EXPECT_EQ(state3.get_kind(), promille::Kind::FIXED);
}
