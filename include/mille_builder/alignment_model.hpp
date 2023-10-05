#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include <Math/EulerAngles.h>
#include <Math/Point3D.h>
#include <Math/Rotation3D.h>
#include <Math/RotationZYX.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <mille_builder/euler_angles.hpp>

/*
 * Rotation matrix, see https://mathworld.wolfram.com/EulerAngles.html (48..50)
 * for definitions.
 */

namespace mb
{

using ROOT::Math::Rotation3D;
using ROOT::Math::RotationZYX;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

template<typename T, std::size_t Nglobal, std::size_t Nlocal, typename PointType, typename VectorType, typename... ExtraArgs>
struct alignment_model_base
{
    static const int n_globals = Nglobal;
    static const int n_locals = Nlocal;

    auto set_extras(ExtraArgs... args) -> void { update_extras(args...); }

    auto get_distance(PointType base, VectorType track) -> T
    {
        calc_derivatives(base, track);
        return calc_distance(base, track);
    }

    auto get_distance(PointType base, VectorType track, ExtraArgs... args) -> T
    {
        update_extras(args...);
        calc_derivatives(base, track);
        return calc_distance(base, track);
    }

    auto transform(PointType point) const -> PointType { return calc_transform(point); }

    auto global_derivative(int variable_number) const -> T
    {
        assert(variable_number > 0 and variable_number <= Nglobal);
        return global_derivatives[variable_number];
    }

    auto local_derivative(int variable_number) const -> T
    {
        assert(variable_number > 0 and variable_number <= Nlocal);
        return local_derivatives[variable_number];
    }

    auto print() const -> void
    {
        for (int i = 1; i <= Nglobal; ++i) {
            std::cout << " Global derivative " << i << " = " << global_derivative(i) << '\n';
        }
        // std::cout << "#1 Gt" << ROOT::Math::XYZVector(g_x, g_y, g_z) << "  Gr" << ROOT::Math::XYZVector(g_a, g_b, g_c) << "  At"
        //           << ROOT::Math::XYZVector(a_x, a_y, a_z) << "  WM: " << euler::make_rotation_matrix(wm) << "#2 Sl"
        //           << ROOT::Math::XYZVector(s_x, s_y, s_z) << "  Bt" << ROOT::Math::XYZVector(b_x, b_y, 0) << "  Tt"
        //           << ROOT::Math::XYZVector(t_x, t_y, 1) << '\n';
    }

  protected:
    template<int variable_number>
    auto set_global_derivative() -> T&
    {
        static_assert(variable_number > 0 and variable_number <= Nglobal, "Index out of bounds");
        return global_derivatives[variable_number];
    }

    template<int variable_number>
    auto set_local_derivative() -> T&
    {
        static_assert(variable_number > 0 and variable_number <= Nlocal, "Index out of bounds");
        return local_derivatives[variable_number];
    }

  protected:
    constexpr alignment_model_base() = default;

  private:
    std::array<T, Nglobal + 1> global_derivatives {0};
    std::array<T, Nlocal + 1> local_derivatives {0};

    virtual auto update_extras(ExtraArgs... args) -> void = 0;
    virtual auto calc_distance(PointType base, VectorType track) -> T = 0;

    virtual auto calc_derivatives(PointType base, VectorType track) -> void = 0;

    virtual auto calc_transform(PointType point) const -> PointType = 0;
};

template<typename T, template<class> class R>
struct sts_residua final : alignment_model_base<T, 6, 4, XYZPoint, XYZVector, XYZPoint, XYZVector>
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

    XYZPoint current_straw_base;
    XYZVector current_straw_dir;

    mb::euler::euler_base<T> wm;

    sts_residua(T gx, T gy, T gz, T ga, T gb, T gc)
        : g_x(gx)
        , g_y(gy)
        , g_z(gz)
        , g_a(ga)
        , g_b(gb)
        , g_c(gc)
        , wm(R<T>(0, 0, 0))
    {
    }

    auto set_params(T ax, T ay, T az, T alpha, T beta, T gamma) -> void
    {
        a_x = ax;
        a_y = ay;
        a_z = az;
        wm = R<T>(alpha, beta, gamma);
    }

    auto update_extras(XYZPoint straw_base, XYZVector straw_dir) -> void override
    {
        current_straw_base = straw_base;
        current_straw_dir = straw_dir;
    }

    auto calc_distance(XYZPoint track_base, XYZVector track_dir) -> T override
    {
        return std::abs((track_base - current_straw_base).Dot(track_dir.Cross(current_straw_dir)) / track_dir.Cross(current_straw_dir).R());
    }

    auto calc_transform(XYZPoint point) const -> XYZPoint { return XYZPoint(); }

    auto calc_derivatives(XYZPoint track_base, XYZVector track_dir) -> void override
    {
        // T sx_, T sy_, T sz_, T bx_, T by_, T tx_, T ty_
        auto b_x = track_base.X();
        auto b_y = track_base.Y();
        auto t_x = track_dir.X();
        auto t_y = track_dir.Y();
        auto s_x = current_straw_base.X();
        auto s_y = current_straw_base.Y();
        auto s_z = current_straw_base.Z();

        // Manually optimized shortcuts for frequently appearing expressions.
        // Touch it on your own risk.
        auto common_0 = wm.R12 - g_b * t_x * wm.R12 - g_c * wm.R22 + g_a * t_x * wm.R22 + g_b * wm.R32 - t_x * wm.R32;
        auto common_1 = g_c * wm.R12 + g_b * t_y * wm.R12 - wm.R22 - g_a * t_y * wm.R22 - g_a * wm.R32 + t_y * wm.R32;
        auto common_2 = -(g_c * t_x * wm.R12) - t_y * wm.R12 + t_x * wm.R22 + g_c * t_y * wm.R22 + g_a * t_x * wm.R32 - g_b * t_y * wm.R32;
        auto common_3 = -a_x + b_x - g_x - s_x * (wm.R11 - g_c * wm.R21 + g_b * wm.R31) - s_y * (wm.R12 - g_c * wm.R22 + g_b * wm.R32)
            - s_z * (wm.R13 - g_c * wm.R23 + g_b * wm.R33);
        auto common_4 = -a_y + b_y - g_y - s_x * (-(g_c * wm.R11) + wm.R21 + g_a * wm.R31) - s_y * (-(g_c * wm.R12) + wm.R22 + g_a * wm.R32)
            - s_z * (-(g_c * wm.R13) + wm.R23 + g_a * wm.R33);
        auto common_5 = -a_z - g_z - s_x * (g_b * wm.R11 - g_a * wm.R21 + wm.R31) - s_y * (g_b * wm.R12 - g_a * wm.R22 + wm.R32)
            - s_z * (g_b * wm.R13 - g_a * wm.R23 + wm.R33);

        auto common_10 = pow(common_0, 2) + pow(common_1, 2) + pow(common_2, 2);
        auto common_11 = common_2 * common_5 + common_0 * common_4 + common_1 * common_3;

        auto common_20 = sqrt(common_10);
        auto common_21 = pow(common_10, 3. / 2.);
        auto common_22 = fabs(common_11);

        // dgx
        this->template set_global_derivative<1>() =
            ((-(g_c * wm.R12) - g_b * t_y * wm.R12 + wm.R22 + g_a * t_y * wm.R22 + g_a * wm.R32 - t_y * wm.R32) * (common_11))
            / (common_20 * common_22);

        // dgy
        this->template set_global_derivative<2>() =
            ((-wm.R12 + g_b * t_x * wm.R12 + g_c * wm.R22 - g_a * t_x * wm.R22 - g_b * wm.R32 + t_x * wm.R32) * (common_11))
            / (common_20 * common_22);

        // dgz
        this->template set_global_derivative<3>() =
            ((g_c * t_x * wm.R12 + t_y * wm.R12 - t_x * wm.R22 - g_c * t_y * wm.R22 - g_a * t_x * wm.R32 + g_b * t_y * wm.R32)
             * (common_11))
            / (common_20 * common_22);

        // dga
        this->template set_global_derivative<4>() =
            (((s_x * wm.R21 + s_y * wm.R22 + s_z * wm.R23) * (common_2) + (common_0) * (-(s_x * wm.R31) - s_y * wm.R32 - s_z * wm.R33)
              + t_x * wm.R32 * (common_5) + t_x * wm.R22 * (common_4) + (-(t_y * wm.R22) - wm.R32) * (common_3))
             * (common_11))
                / (common_20 * common_22)
            - ((2 * t_x * wm.R22 * (common_0) + 2 * (-(t_y * wm.R22) - wm.R32) * (common_1) + 2 * t_x * wm.R32 * (common_2)) * common_22)
                / (2 * common_21);

        // dgb
        this->template set_global_derivative<5>() =
            (((-(s_x * wm.R11) - s_y * wm.R12 - s_z * wm.R13) * (common_2) + (common_1) * (-(s_x * wm.R31) - s_y * wm.R32 - s_z * wm.R33)
              - t_y * wm.R32 * (common_5) + (-(t_x * wm.R12) + wm.R32) * (common_4) + t_y * wm.R12 * (common_3))
             * (common_11))
                / (common_20 * common_22)
            - ((2 * (-(t_x * wm.R12) + wm.R32) * (common_0) + 2 * t_y * wm.R12 * (common_1)-2 * t_y * wm.R32 * (common_2)) * common_22)
                / (2 * common_21);

        // dgc
        this->template set_global_derivative<6>() =
            (((s_x * wm.R11 + s_y * wm.R12 + s_z * wm.R13) * (common_0) + (s_x * wm.R21 + s_y * wm.R22 + s_z * wm.R23) * (common_1)
              + (-(t_x * wm.R12) + t_y * wm.R22) * (common_5)-wm.R22 * (common_4) + wm.R12 * (common_3))
             * (common_11))
                / (common_20 * common_22)
            - ((-2 * wm.R22 * (common_0) + 2 * wm.R12 * (common_1) + 2 * (-(t_x * wm.R12) + t_y * wm.R22) * (common_2)) * common_22)
                / (2 * common_21);

        // dbx
        this->template set_local_derivative<1>() = ((common_1) * (common_11)) / (common_20 * common_22);

        // dby
        this->template set_local_derivative<2>() = ((common_0) * (common_11)) / (common_20 * common_22);

        // dtx
        this->template set_local_derivative<3>() =
            (((-(g_c * wm.R12) + wm.R22 + g_a * wm.R32) * (common_5) + (-(g_b * wm.R12) + g_a * wm.R22 - wm.R32) * (common_4))
             * (common_11))
                / (common_20 * common_22)
            - ((2 * (-(g_b * wm.R12) + g_a * wm.R22 - wm.R32) * (common_0) + 2 * (-(g_c * wm.R12) + wm.R22 + g_a * wm.R32) * (common_2))
               * common_22)
                / (2 * common_21);

        // dty
        this->template set_local_derivative<4>() =
            (((-wm.R12 + g_c * wm.R22 - g_b * wm.R32) * (common_5) + (g_b * wm.R12 - g_a * wm.R22 + wm.R32) * (common_3)) * (common_11))
                / (common_20 * common_22)
            - ((2 * (g_b * wm.R12 - g_a * wm.R22 + wm.R32) * (common_1) + 2 * (-wm.R12 + g_c * wm.R22 - g_b * wm.R32) * (common_2))
               * common_22)
                / (2 * common_21);
    }

    // auto print() const -> void
    // {
    // std::cout << "#1 Gt" << ROOT::Math::XYZVector(g_x, g_y, g_z) << "  Gr" << ROOT::Math::XYZVector(g_a, g_b, g_c) << "  At"
    // << ROOT::Math::XYZVector(a_x, a_y, a_z) << "  WM: " << euler::make_rotation_matrix(wm) << "#2 Sl"
    // << ROOT::Math::XYZVector(s_x, s_y, s_z) << "  Bt" << ROOT::Math::XYZVector(b_x, b_y, 0) << "  Tt"
    // << ROOT::Math::XYZVector(t_x, t_y, 1) << '\n';
    // }
};

}  // namespace mb
