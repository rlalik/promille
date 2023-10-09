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

namespace mb_tests
{

using ROOT::Math::Rotation3D;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

template<typename T, typename U>
struct dummy_residual_model final : promille::residual_model_base<T, 6, 4, U, U, U>
{
    // gloal translational corrections
    T g_x {0};
    T g_y {0};
    T g_z {0};

    // gloal rotational corrections
    T g_a {0};
    T g_b {0};
    T g_c {0};

    // translational alignment of straws
    T a_x {0};
    T a_y {0};
    T a_z {0};

    U current_straw_base;
    U current_straw_dir;

    dummy_residual_model(T gx, T gy, T gz, T ga, T gb, T gc)
        : g_x(gx)
        , g_y(gy)
        , g_z(gz)
        , g_a(ga)
        , g_b(gb)
        , g_c(gc)
    {
    }

    auto set_params(T ax, T ay, T az, T alpha, T beta, T gamma) -> void
    {
        a_x = ax;
        a_y = ay;
        a_z = az;
    }

    auto update_extras(U /* unused */) -> void override {}

    auto calc_residual(U track_base, U track_dir) -> T override { return 100; }

    auto calc_derivatives(U track_base, U track_dir) -> void override
    {
        auto b_x = track_base;
        auto t_x = track_dir;

        // dgx
        this->template set_global_derivative<1>() = 1;

        // dgy
        this->template set_global_derivative<2>() = 2;

        // dgz
        this->template set_global_derivative<3>() = 3;

        // dga
        this->template set_global_derivative<4>() = 4;

        // dgb
        this->template set_global_derivative<5>() = 5;

        // dgc
        this->template set_global_derivative<6>() = 6;

        // dbx
        this->template set_local_derivative<1>() = 0.1;

        // dby
        this->template set_local_derivative<2>() = 0.2;

        // dtx
        this->template set_local_derivative<3>() = 0.3;

        // dty
        this->template set_local_derivative<4>() = 0.4;
    }
};

}  // namespace mb_tests
