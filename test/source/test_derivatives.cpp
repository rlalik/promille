#include <StrawAlignment/StrawAlignment.hpp>
#include <gtest/gtest.h>

static double accuracy = 0.00001;

TEST(Derivatives, Global)
{
    std::vector<std::pair<std::pair<std::array<double, 12>, std::array<double, 8>>, std::array<double, 11>>> data = {
        {{{0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}},
         {0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0}},
         {{{0, 1, 2, 3, 4, 5, 1, 2, 3, 0, 0, 0}, {6, 7, 8, 9, 10, 11, 12, 13}},
         {0, 0, 0, -3.77752, 1.54111, 1.3419, 0.434116, -0.33644, -0.835674, -0.199162, 0.15435}}};

    for (const auto& d : data) {
        SA::derivatives<double, SA::euler::zyz> derivs(d.first.first[0],
                                                       d.first.first[1],
                                                       d.first.first[2],
                                                       d.first.first[3],
                                                       d.first.first[4],
                                                       d.first.first[5],
                                                       d.first.first[6],
                                                       d.first.first[7],
                                                       d.first.first[8],
                                                       d.first.first[9],
                                                       d.first.first[10],
                                                       d.first.first[11]);

        derivs.update(d.first.second[0],
                      d.first.second[1],
                      d.first.second[2],
                      d.first.second[3],
                      d.first.second[4],
                      d.first.second[5],
                      d.first.second[6],
                      d.first.second[7]);

        EXPECT_NEAR(derivs.dr_dgx(), d.second[0], accuracy);
        EXPECT_NEAR(derivs.dr_dgy(), d.second[1], accuracy);
        EXPECT_NEAR(derivs.dr_dgz(), d.second[2], accuracy);

        EXPECT_NEAR(derivs.dr_dga(), d.second[3], accuracy);
        EXPECT_NEAR(derivs.dr_dgb(), d.second[4], accuracy);
        EXPECT_NEAR(derivs.dr_dgc(), d.second[5], accuracy);

        EXPECT_NEAR(derivs.dr_dbx(), d.second[6], accuracy);
        EXPECT_NEAR(derivs.dr_dby(), d.second[7], accuracy);
        EXPECT_NEAR(derivs.dr_dbz(), d.second[8], accuracy);

        EXPECT_NEAR(derivs.dr_dtx(), d.second[9], accuracy);
        EXPECT_NEAR(derivs.dr_dty(), d.second[10], accuracy);
    }
}

TEST(Derivatives, Local)
{
    // SA::derivatives<float, SA::euler::zyz> derivs(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13);
    //
    // EXPECT_NEAR(derivs.dr_dbx(), -0.72752, accuracy);
    // EXPECT_NEAR(derivs.dr_dby(), 0.67908, accuracy);
    // EXPECT_NEAR(derivs.dr_dbz(), -0.0977986, accuracy);
    //
    // EXPECT_NEAR(derivs.dr_dtx(), 0.0849619, accuracy);
    // EXPECT_NEAR(derivs.dr_dty(), -0.0793049, accuracy);
}
