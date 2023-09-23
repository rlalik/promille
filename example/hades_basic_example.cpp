#include <StrawAlignment/StrawAlignment.hpp>
#include <TMath.h>

auto main() -> int
{
    auto mille = SA::MilleBuilder("hades_basic_example_alignment_data.bin");

    mille.add_planes_global({0, SA::Kind::FREE},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {TMath::DegToRad() * 0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED});

    mille.add_planes_global({0, SA::Kind::FIXED},
                            {0, SA::Kind::FREE},
                            {0, SA::Kind::FIXED},
                            {TMath::DegToRad() * 90, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED});

    mille.add_planes_global({0, SA::Kind::FIXED},
                            {0, SA::Kind::FREE},
                            {0, SA::Kind::FIXED},
                            {TMath::DegToRad() * 90, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED});

    mille.add_planes_global({0, SA::Kind::FREE},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {TMath::DegToRad() * 0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED});

    mille.add_planes_global({0, SA::Kind::FREE},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {TMath::DegToRad() * 0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED});

    mille.add_planes_global({0, SA::Kind::FIXED},
                            {0, SA::Kind::FREE},
                            {0, SA::Kind::FIXED},
                            {TMath::DegToRad() * 90, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED});

    mille.add_planes_global({0, SA::Kind::FREE},
                            {0, SA::Kind::FREE},
                            {0, SA::Kind::FIXED},
                            {TMath::DegToRad() * 45, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED});

    mille.add_planes_global({0, SA::Kind::FREE},
                            {0, SA::Kind::FREE},
                            {0, SA::Kind::FIXED},
                            {TMath::DegToRad() * -45, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED},
                            {0, SA::Kind::FIXED});

    mille.add_local(0, 0, 1, 2, 3, 4, 0, 0, 0, 0.1, 0.15);
    mille.add_local(1, 1, 2, 3, 4, 5, 0, 0, 0, 0.1, 0.15);
    mille.add_local(2, 2, 3, 4, 5, 6, 0, 0, 0, 0.1, 0.15);
    mille.add_local(3, 0, 1, 2, 3, 4, 0, 0, 0, 0.1, 0.15);
    mille.add_local(4, 1, 2, 3, 4, 5, 0, 0, 0, 0.1, 0.15);
    mille.add_local(5, 2, 3, 4, 5, 6, 0, 0, 0, 0.1, 0.15);
    mille.end();

    mille.add_local(0, 1, 2, 3, 4, 5, 0, 0, 0, 0.1, 0.15);
    mille.add_local(1, 2, 3, 4, 5, 6, 0, 0, 0, 0.1, 0.15);
    mille.add_local(2, 3, 4, 5, 6, 7, 0, 0, 0, 0.1, 0.15);
    mille.end();

    return 0;
}
