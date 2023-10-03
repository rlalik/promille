#include <gtest/gtest.h>
#include <mille_builder/mille_builder.hpp>

TEST(Mille, AddLocals)
{
    mb::mille_builder<mb::euler::zyz> mille("test_", "test");

    mille.add_planes_globals({0, mb::Kind::FIXED},
                             {0, mb::Kind::FIXED},
                             {0, mb::Kind::FIXED},
                             {0, mb::Kind::FIXED},
                             {0, mb::Kind::FIXED},
                             {0, mb::Kind::FIXED},
                             0,
                             0,
                             0,
                             0,
                             0,
                             0,
                             mb::Kind::FIXED,
                             mb::Kind::FIXED,
                             mb::Kind::FIXED,
                             mb::Kind::FIXED);

    mille.add_local(0, 0, 0, 0, 0, 0, 0, 0, 0.1f, 0.15f);

    mille.end();
}
