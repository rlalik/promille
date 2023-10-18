#include <gtest/gtest.h>
#include <promille/promille.hpp>

TEST(Parameters, Parameters)
{
    auto param = promille::global_parameter<double>(1, 10.0);

    EXPECT_EQ(param.id, 1);
    EXPECT_EQ(param.value, 10.0);
    EXPECT_EQ(param.description, "");
    EXPECT_EQ(param.is_anywere_free, false);

    EXPECT_EQ(param, 10.0);
}

TEST(Parameters, ParameterKind)
{
    auto param = promille::global_parameter<double>(1, 10.0);

    auto state1 = param.make_parameter_kind();
    EXPECT_EQ(state1.get_kind(), promille::Kind::FIXED);

    auto state2 = param.make_parameter_kind(promille::Kind::FIXED);
    EXPECT_EQ(state2.get_kind(), promille::Kind::FIXED);

    auto state3 = param.make_parameter_kind(promille::Kind::FREE);
    EXPECT_EQ(state3.get_kind(), promille::Kind::FREE);

    state3.set_kind(promille::Kind::FIXED);
    EXPECT_EQ(state3.get_kind(), promille::Kind::FIXED);
}
