#pragma once

// internal
#include <promille/euler_angles.hpp>
#include <promille/promille.hpp>

// ROOT
#include <Math/Point3D.h>
#include <Math/Rotation3D.h>
#include <Math/Vector3D.h>

/*
 * Rotation matrix, see https://mathworld.wolfram.com/EulerAngles.html (48..50)
 * for definitions.
 */

namespace promille_tests
{

using ROOT::Math::Rotation3D;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

template<typename T, typename U, size_t Ng = 3, size_t Nl = 2>
struct dummy_residual_model final : promille::residual_model_base<T, Ng, Nl, U, U>
{
    // gloal translational corrections
    T gp1 {0};
    T gp2 {0};
    T gp3 {0};

    dummy_residual_model(T gp1, T gp2, T gp3)
        : gp1(gp1)
        , gp2(gp2)
        , gp3(gp3)
    {
    }

    auto calc_model(U track_base, U track_dir) -> T override
    {
        auto b_x = track_base;
        auto t_x = track_dir;

        this->template set_global_derivative<0>() = 1;

        this->template set_global_derivative<1>() = 2;

        this->template set_global_derivative<2>() = 3;

        this->template set_local_derivative<0>() = 0.1;

        this->template set_local_derivative<1>() = 0.2;

        return 100;
    }
};

}  // namespace promille_tests
