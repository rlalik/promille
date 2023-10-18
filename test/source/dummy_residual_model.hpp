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

template<typename T, typename U, size_t Nl = 2, size_t Ng = 3>
struct dummy_residual_model final : promille::residual_model_base<T, Nl, Ng>
{
    // gloal translational corrections
    T gp1 {0};
    T gp2 {0};
    T gp3 {0};

    dummy_residual_model(T gp1, T gp2, T gp3)
        : promille::residual_model_base<T, Nl, Ng>()
        , gp1(gp1)
        , gp2(gp2)
        , gp3(gp3)
    {
    }

    auto residual() const -> T override { return 100.; }

    auto compute(U track_base, U track_dir) -> void
    {
        auto b_x = track_base;
        auto t_x = track_dir;

        this->template set_global_derivative<0>() = gp1;

        this->template set_global_derivative<1>() = gp2;

        this->template set_global_derivative<2>() = gp3;

        this->template set_local_derivative<0>() = 0.1;

        this->template set_local_derivative<1>() = 0.2;

        this->residual = 100;
    }
};

template<typename T, typename U, size_t Nl = 2, size_t Ng = 4>
struct dummy_residual_model_2 final : promille::residual_model_base<T, Nl, Ng>
{
    // gloal translational corrections
    T gp1 {0};
    T gp2 {0};
    T gp3 {0};
    T gp4 {0};

    dummy_residual_model_2(T gp1, T gp2, T gp3, T gp4)
        : promille::residual_model_base<T, Nl, Ng>()
        , gp1(gp1)
        , gp2(gp2)
        , gp3(gp3)
        , gp4(gp4)
    {
    }

    auto residual() const -> T override { return 100.; }

    auto compute(U track_base, U track_dir) -> void
    {
        auto b_x = track_base;
        auto t_x = track_dir;

        this->template set_global_derivative<0>() = gp1;

        this->template set_global_derivative<1>() = gp2;

        this->template set_global_derivative<2>() = gp3;

        this->template set_local_derivative<0>() = 0.1;

        this->template set_local_derivative<1>() = 0.2;

        this->residual = 100;
    }
};

}  // namespace promille_tests
