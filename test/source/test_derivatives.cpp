#include <StrawAlignment/StrawAlignment.hpp>
#include <gtest/gtest.h>

TEST(Derivatives, Global)
{
    SA::derivatives<double> derivs(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13);

    EXPECT_FLOAT_EQ(derivs.dr_dpsi(), 8.22009);
    EXPECT_FLOAT_EQ(derivs.dr_dtheta(), 7.17463);
    EXPECT_FLOAT_EQ(derivs.dr_dphi(), 0.794929);

    EXPECT_FLOAT_EQ(derivs.dr_dXa(), -0.746818);
    EXPECT_FLOAT_EQ(derivs.dr_dYa(), 0.0);
    EXPECT_FLOAT_EQ(derivs.dr_dZa(), 0.665028);
}

TEST(Derivatives, Local)
{
    SA::derivatives<float> derivs(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13);

    EXPECT_FLOAT_EQ(derivs.dr_dbx(), -0.72752);
    EXPECT_FLOAT_EQ(derivs.dr_dby(), 0.67908);
    EXPECT_FLOAT_EQ(derivs.dr_dbz(), -0.0977986);

    EXPECT_FLOAT_EQ(derivs.dr_dtx(), 0.0849619);
    EXPECT_FLOAT_EQ(derivs.dr_dty(), -0.0793049);
}
