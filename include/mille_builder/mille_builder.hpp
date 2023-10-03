#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
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

#include "Mille.h"

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

namespace bits
{

constexpr auto size_t2int(size_t val) -> int
{
    return (val <= std::numeric_limits<int>::max()) ? static_cast<int>(val) : -1;
}

constexpr auto int2size_t(int val) -> size_t
{
    return (val < 0) ? __SIZE_MAX__ : static_cast<size_t>(val);
}

constexpr auto f2d(float x) -> double
{
    return static_cast<double>(x);
}

constexpr auto d2f(double x) -> double
{
    return static_cast<float>(x);
}

}  // namespace bits

namespace geom
{

inline auto distance(XYZPoint base1, XYZVector dir1, XYZPoint base2, XYZVector dir2) -> XYZPoint::Scalar
{
    return std::abs((base1 - base2).Dot(dir1.Cross(dir2)) / dir1.Cross(dir2).R());
}

}  // namespace geom

template<typename T, template<class> class R>
struct derivatives
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

    // local system straw coordinates
    T s_x {0};
    T s_y {0};
    T s_z {0};

    // track base
    T b_x {0};
    T b_y {0};

    // track direction
    T t_x {0};
    T t_y {0};

    T common_0;
    T common_1;
    T common_2;
    T common_3;
    T common_4;
    T common_5;
    T common_10;
    T common_11;
    T common_20;
    T common_21;
    T common_22;

    mb::euler::euler_base<T> wm;

    derivatives(T gx, T gy, T gz, T ga, T gb, T gc, T ax, T ay, T az, T alpha, T beta, T gamma)
        : g_x(gx)
        , g_y(gy)
        , g_z(gz)
        , g_a(ga)
        , g_b(gb)
        , g_c(gc)
        , a_x(ax)
        , a_y(ay)
        , a_z(az)
        , wm(R<T>(alpha, beta, gamma))
    {
    }

    auto update(T sx_, T sy_, T sz_, T bx_, T by_, T tx_, T ty_) -> void
    {
        s_x = sx_;
        s_y = sy_;
        s_z = sz_;
        b_x = bx_;
        b_y = by_;
        t_x = tx_;
        t_y = ty_;

        // Manually optimized shortcuts for frequently appearing expressions.
        // Touch it on your own risk.
        common_0 = wm.R12 - g_b * t_x * wm.R12 - g_c * wm.R22 + g_a * t_x * wm.R22 + g_b * wm.R32 - t_x * wm.R32;
        common_1 = g_c * wm.R12 + g_b * t_y * wm.R12 - wm.R22 - g_a * t_y * wm.R22 - g_a * wm.R32 + t_y * wm.R32;
        common_2 = -(g_c * t_x * wm.R12) - t_y * wm.R12 + t_x * wm.R22 + g_c * t_y * wm.R22 + g_a * t_x * wm.R32 - g_b * t_y * wm.R32;
        common_3 = -a_x + b_x - g_x - s_x * (wm.R11 - g_c * wm.R21 + g_b * wm.R31) - s_y * (wm.R12 - g_c * wm.R22 + g_b * wm.R32)
            - s_z * (wm.R13 - g_c * wm.R23 + g_b * wm.R33);
        common_4 = -a_y + b_y - g_y - s_x * (-(g_c * wm.R11) + wm.R21 + g_a * wm.R31) - s_y * (-(g_c * wm.R12) + wm.R22 + g_a * wm.R32)
            - s_z * (-(g_c * wm.R13) + wm.R23 + g_a * wm.R33);
        common_5 = -a_z - g_z - s_x * (g_b * wm.R11 - g_a * wm.R21 + wm.R31) - s_y * (g_b * wm.R12 - g_a * wm.R22 + wm.R32)
            - s_z * (g_b * wm.R13 - g_a * wm.R23 + wm.R33);

        common_10 = pow(common_0, 2) + pow(common_1, 2) + pow(common_2, 2);
        common_11 = common_2 * common_5 + common_0 * common_4 + common_1 * common_3;

        common_20 = sqrt(common_10);
        common_21 = pow(common_10, 3. / 2.);
        common_22 = fabs(common_11);
    }

    constexpr auto dr_dgx() const -> T
    {
        return ((-(g_c * wm.R12) - g_b * t_y * wm.R12 + wm.R22 + g_a * t_y * wm.R22 + g_a * wm.R32 - t_y * wm.R32) * (common_11))
            / (common_20 * common_22);
    }

    constexpr auto dr_dgy() const -> T
    {
        return ((-wm.R12 + g_b * t_x * wm.R12 + g_c * wm.R22 - g_a * t_x * wm.R22 - g_b * wm.R32 + t_x * wm.R32) * (common_11))
            / (common_20 * common_22);
    }

    constexpr auto dr_dgz() const -> T
    {
        return ((g_c * t_x * wm.R12 + t_y * wm.R12 - t_x * wm.R22 - g_c * t_y * wm.R22 - g_a * t_x * wm.R32 + g_b * t_y * wm.R32)
                * (common_11))
            / (common_20 * common_22);
    }

    constexpr auto dr_dga() const -> T
    {
        return (((s_x * wm.R21 + s_y * wm.R22 + s_z * wm.R23) * (common_2) + (common_0) * (-(s_x * wm.R31) - s_y * wm.R32 - s_z * wm.R33)
                 + t_x * wm.R32 * (common_5) + t_x * wm.R22 * (common_4) + (-(t_y * wm.R22) - wm.R32) * (common_3))
                * (common_11))
            / (common_20 * common_22)
            - ((2 * t_x * wm.R22 * (common_0) + 2 * (-(t_y * wm.R22) - wm.R32) * (common_1) + 2 * t_x * wm.R32 * (common_2)) * common_22)
            / (2 * common_21);
    }

    constexpr auto dr_dgb() const -> T
    {
        return (((-(s_x * wm.R11) - s_y * wm.R12 - s_z * wm.R13) * (common_2) + (common_1) * (-(s_x * wm.R31) - s_y * wm.R32 - s_z * wm.R33)
                 - t_y * wm.R32 * (common_5) + (-(t_x * wm.R12) + wm.R32) * (common_4) + t_y * wm.R12 * (common_3))
                * (common_11))
            / (common_20 * common_22)
            - ((2 * (-(t_x * wm.R12) + wm.R32) * (common_0) + 2 * t_y * wm.R12 * (common_1)-2 * t_y * wm.R32 * (common_2)) * common_22)
            / (2 * common_21);
    }

    constexpr auto dr_dgc() const -> T
    {
        return (((s_x * wm.R11 + s_y * wm.R12 + s_z * wm.R13) * (common_0) + (s_x * wm.R21 + s_y * wm.R22 + s_z * wm.R23) * (common_1)
                 + (-(t_x * wm.R12) + t_y * wm.R22) * (common_5)-wm.R22 * (common_4) + wm.R12 * (common_3))
                * (common_11))
            / (common_20 * common_22)
            - ((-2 * wm.R22 * (common_0) + 2 * wm.R12 * (common_1) + 2 * (-(t_x * wm.R12) + t_y * wm.R22) * (common_2)) * common_22)
            / (2 * common_21);
    }

    constexpr auto dr_dbx() const -> T { return ((common_1) * (common_11)) / (common_20 * common_22); }

    constexpr auto dr_dby() const -> T { return ((common_0) * (common_11)) / (common_20 * common_22); }

    constexpr auto dr_dtx() const -> T
    {
        return (((-(g_c * wm.R12) + wm.R22 + g_a * wm.R32) * (common_5) + (-(g_b * wm.R12) + g_a * wm.R22 - wm.R32) * (common_4))
                * (common_11))
            / (common_20 * common_22)
            - ((2 * (-(g_b * wm.R12) + g_a * wm.R22 - wm.R32) * (common_0) + 2 * (-(g_c * wm.R12) + wm.R22 + g_a * wm.R32) * (common_2))
               * common_22)
            / (2 * common_21);
    }

    constexpr auto dr_dty() const -> T
    {
        return (((-wm.R12 + g_c * wm.R22 - g_b * wm.R32) * (common_5) + (g_b * wm.R12 - g_a * wm.R22 + wm.R32) * (common_3)) * (common_11))
            / (common_20 * common_22)
            - ((2 * (g_b * wm.R12 - g_a * wm.R22 + wm.R32) * (common_1) + 2 * (-wm.R12 + g_c * wm.R22 - g_b * wm.R32) * (common_2))
               * common_22)
            / (2 * common_21);
    }

    auto print() const -> void
    {
        std::cout << "#1 Gt" << ROOT::Math::XYZVector(g_x, g_y, g_z) << "  Gr" << ROOT::Math::XYZVector(g_a, g_b, g_c) << "  At"
                  << ROOT::Math::XYZVector(a_x, a_y, a_z) << "  WM: " << euler::make_rotation_matrix(wm) << "#2 Sl"
                  << ROOT::Math::XYZVector(s_x, s_y, s_z) << "  Bt" << ROOT::Math::XYZVector(b_x, b_y, 0) << "  Tt"
                  << ROOT::Math::XYZVector(t_x, t_y, 1) << '\n';
    }
};

enum class Kind
{
    FIXED,
    FREE
};

inline auto to_kind(int v) -> Kind
{
    return v ? Kind::FREE : Kind::FIXED;
}

inline auto operator<<(std::ostream& ofs, Kind rhs) -> std::ostream&
{
    return ofs << std::setw(2) << (rhs == mb::Kind::FREE ? '~' : 'x');
}

template<typename T>
struct parameter
{
    T value;
    Kind kind;

    auto is_free() const -> bool { return kind == Kind::FREE; }
    operator T() const { return value; }
};

template<typename T>
auto operator<<(std::ostream& ofs, const mb::parameter<T>& rhs) -> std::ostream&
{
    return ofs << rhs.kind << std::setw(12) << rhs.value;
}

template<typename T>
struct global_parameters
{
    global_parameters(parameter<T> gx, parameter<T> gy, parameter<T> gz, parameter<T> ga, parameter<T> gb, parameter<T> gc)
    {
        params.reserve(6);
        params.push_back(gx);
        params.push_back(gy);
        params.push_back(gz);
        params.push_back(ga);
        params.push_back(gb);
        params.push_back(gc);
    }

    std::vector<parameter<T>> params;

    auto gx() -> parameter<T>& { return params[0]; }
    auto gy() -> parameter<T>& { return params[1]; }
    auto gz() -> parameter<T>& { return params[2]; }
    auto ga() -> parameter<T>& { return params[3]; }
    auto gb() -> parameter<T>& { return params[4]; }
    auto gc() -> parameter<T>& { return params[5]; }

    auto gx() const -> const parameter<T>& { return params[0]; }
    auto gy() const -> const parameter<T>& { return params[1]; }
    auto gz() const -> const parameter<T>& { return params[2]; }
    auto ga() const -> const parameter<T>& { return params[3]; }
    auto gb() const -> const parameter<T>& { return params[4]; }
    auto gc() const -> const parameter<T>& { return params[5]; }

    auto dump_pede_param(int layer, std::ostream& ofs) -> void
    {
        for (int i = 0; i < 6; ++i) {
            const auto& p = params[i];
            if (p.is_free()) {
                ofs << std::setw(10) << layer * 10 + i + 1 << std::setw(20) << p.value << std::setw(10) << 0.0 << '\n';
            }
        }
    }
};

struct local_parameters
{
    local_parameters(Kind lx0, Kind ly0, Kind ltx, Kind lty)
        : params {lx0, ly0, ltx, lty}
    {
    }

    std::array<Kind, 4> params;

    auto lx0() -> Kind& { return params[0]; }
    auto ly0() -> Kind& { return params[1]; }
    auto ltx() -> Kind& { return params[2]; }
    auto lty() -> Kind& { return params[3]; }

    auto lx0() const -> const Kind& { return params[0]; }
    auto ly0() const -> const Kind& { return params[1]; }
    auto ltx() const -> const Kind& { return params[2]; }
    auto lty() const -> const Kind& { return params[3]; }

    auto dump_pede_param(size_t layer, std::ostream& ofs) -> void
    {
        for (decltype(params)::size_type i = 0; i < 4; ++i) {
            const auto& p = params[i];
            ofs << std::setw(10) << layer * 10 + i + 1 << std::setw(20) << (p == Kind::FREE ? '~' : 'x') << std::setw(10) << 0.0 << '\n';
        }
    }
};

template<template<class> class R>
class mille_builder
{
  public:
    mille_builder(const char* prefix, const char* outFileName, bool asBinary = true, bool writeZero = false)
        : mille_prefix(prefix)
        , mille(outFileName, asBinary, writeZero)
    {
    }
    ~mille_builder() {}

    /**
     * Add tracking planes alignment definitions. Alignment is defined as fixed rotational and transformational elements, and rotational and
     * transformational corrections.
     *
     * @param gx transformational global correction for x
     * @param gy transformational global correction for y
     * @param gz transformational global correction for z
     * @param ga rotational global correction d_alpha
     * @param gb rotational global correction d_beta
     * @param gc rotational global correction d_alpha
     * @param ax alignment x-component
     * @param ay alignment y-component
     * @param az alignment y-component
     * @param alpha rotationa 1st angle
     * @param beta rotationa 2nd angle
     * @param gamma rotationa 3rd angle
     * @param lx0 is base-x relevant for this plane
     * @param ly0 is base-x relevant for this plane
     * @param ltx is base-x relevant for this plane
     * @param lty is base-x relevant for this plane
     */
    auto add_planes_globals(parameter<float> gx,
                            parameter<float> gy,
                            parameter<float> gz,
                            parameter<float> ga,
                            parameter<float> gb,
                            parameter<float> gc,
                            float ax,
                            float ay,
                            float az,
                            float alpha,
                            float beta,
                            float gamma,
                            Kind lx0,
                            Kind ly0,
                            Kind ltx,
                            Kind lty) -> void
    {
        layers_global_pars.emplace_back(gx, gy, gz, ga, gb, gc);
        layers_derivatives.emplace_back(gx.value, gy.value, gz.value, ga.value, gb.value, gc.value, ax, ay, az, alpha, beta, gamma);
        alignment_rotation_corrections.emplace_back(1, -gc, gb, gc, 1, ga, -gb, ga, 1);
        layers_local_pars.emplace_back(lx0, ly0, ltx, lty);
    }

    /**
     * Add set of local variables for given layer of straws.
     *
     * @param layer layer number
     * @param bx base vector x-coordinate
     * @param by base vector y-coordinate
     * @param bz base vector z-coordinate
     * @param tx direction vector x-component
     * @param ty direction vector y-component
     * @param s_x wire position x
     * @param s_y wire position y
     * @param s_z wire position z
     */
    auto add_local(size_t layer, float bx, float by, float tx, float ty, float sx, float sy, float sz, float dr, float sigma) -> void
    {
        assert(layer < layers_global_pars.size());

        auto& derivs = layers_derivatives[layer];
        derivs.update(sx, sy, sz, bx, by, tx, ty);

        std::array<float, 4> local_derivatives;
        int local_cnt = 0;

        const auto& current_layer_locals = layers_local_pars[layer];
        if (current_layer_locals.lx0() == Kind::FREE) {
            local_derivatives[bits::int2size_t(local_cnt)] = derivs.dr_dbx();
            local_cnt++;
        }
        if (current_layer_locals.ly0() == Kind::FREE) {
            local_derivatives[bits::int2size_t(local_cnt)] = derivs.dr_dby();
            local_cnt++;
        }
        if (current_layer_locals.ltx() == Kind::FREE) {
            local_derivatives[bits::int2size_t(local_cnt)] = derivs.dr_dtx();
            local_cnt++;
        }
        if (current_layer_locals.lty() == Kind::FREE) {
            local_derivatives[bits::int2size_t(local_cnt)] = derivs.dr_dty();
            local_cnt++;
        }

        std::array<float, 6> global_derivatives;
        std::array<int, 6> global_deriv_index;
        int global_cnt = 0;
        const auto param_idx_offset = bits::size_t2int(layer) * 10;

        const auto& current_layer_globals = layers_global_pars[layer];
        if (current_layer_globals.gx().is_free()) {
            global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgx();
            global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 1;
            global_cnt++;
        }

        if (current_layer_globals.gy().is_free()) {
            global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgy();
            global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 2;
            global_cnt++;
        }

        if (current_layer_globals.gz().is_free()) {
            global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgz();
            global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 3;
            global_cnt++;
        }

        if (current_layer_globals.ga().is_free()) {
            global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dga();
            global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 4;
            global_cnt++;
        }

        if (current_layer_globals.gb().is_free()) {
            global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgb();
            global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 5;
            global_cnt++;
        }

        if (current_layer_globals.gc().is_free()) {
            global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgc();
            global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 6;
            global_cnt++;
        }

        auto alignment_rotation = euler::make_rotation_matrix(derivs.wm);
        auto sloc = XYZVector(sx, sy, sz);
        auto alignment_translation = XYZVector(derivs.a_x, derivs.a_y, derivs.a_z);
        auto alignment_translation_correction =
            XYZVector(current_layer_globals.gx(), current_layer_globals.gy(), current_layer_globals.gz());
        auto alignment_rotation_correction = alignment_rotation_corrections[layer];
        auto slab = alignment_rotation_correction * alignment_rotation * sloc + alignment_translation + alignment_translation_correction;
        auto srot = alignment_rotation_correction * alignment_rotation * ROOT::Math::XYZVector(0, 1, 0);
        auto res = abs(bits::d2f(geom::distance({bx, by, 0}, {tx, ty, 1.0}, {slab.x(), slab.y(), slab.z()}, srot)) - dr);

        if (verbose) {
            std::cout << " --- LAYER " << layer << " ---\n";
        }
        if (verbose > 1) {
            derivs.print();
            std::cout << "# slab=" << slab << '\n';
        }
        if (verbose) {
            static const size_t value_width = 12;
            std::cout << "+ " << std::setw(3) << layer << "   res=" << std::setw(value_width) << res
                      << "   sigma=" << std::setw(value_width) << sigma << "   NLC=" << local_derivatives.size();
            for (const auto& ld : local_derivatives)
                std::cout << "  " << std::setw(value_width) << ld;

            std::cout << "   NGL=" << global_derivatives.size();
            for (const auto& gd : global_derivatives)
                std::cout << "  " << std::setw(value_width) << gd;
            std::cout << '\n';
        }

        mille.mille(
            local_cnt, local_derivatives.data(), global_cnt, global_derivatives.data(), global_deriv_index.data(), (float)res, sigma);
    }

    /* Call Mille::end()
     */
    auto end() -> void
    {
        if (verbose) {
            std::cout << "--------------------\n";
        }
        mille.end();
    }

    /* Call Mille::kill();
     */
    auto kill() -> void
    {
        if (verbose) {
            std::cout << " KILL  KILL  KILL  KILL\n";
        }
        mille.kill();
    }

    auto write_param_file() -> void
    {
        std::ofstream param_file(mille_prefix + std::string("params.txt"));
        if (!param_file)
            return;

        param_file << "Parameter\n";
        auto max_globals = layers_global_pars.size();
        for (decltype(max_globals) i = 0; i < max_globals; ++i) {
            // if (layers_global_pars[i].is_free())  {}
            layers_global_pars[i].dump_pede_param(i, param_file);
        }
    }
    /**
     * Read straw plane (global paremeters from file.
     *
     * File must start with the keyword describing format of the input. Allowed formats are:
     * * EULERD - 6 parameters: 3 transformation coordinates + 3 Euler angles psi, theta, phi in degrees
     * * EULERR - 6 parameters: 3 transformation coordinates + 3 Euler angles psi, theta, phi in radians
     * * MATRIX - 12 parameters: 3 transformation coordinates + 9 rotational matrix components r11, r12 .. r33
     *
     * * After the parameters, the 6 configuration bits interpreted as transaltion and Euler angles: 0 - free, 0 - fixed;
     * The first line is fllowed by the input data, one set per line.
     */
    auto read_global_data(const std::string& file_name) -> void
    {
        std::ifstream infile(file_name, std::ios::in);

        if (!infile) {
            std::cerr << "File " << file_name << " does not exists.\n";
            abort();
        }

        std::string config_word;
        infile >> config_word;

        if (config_word == "EULERD") {
            float translation_x {}, translation_y {}, translation_z {};
            float rotation_alpha {}, rotation_beta {}, rotation_gamma {};
            int flag_x {}, flag_y {}, flag_z {}, flag_alpha {}, flag_beta {}, flag_gamma {};
            while (true) {
                infile >> translation_x >> translation_y >> translation_z;
                infile >> rotation_alpha >> rotation_beta >> rotation_gamma;
                infile >> flag_x >> flag_y >> flag_z >> flag_alpha >> flag_beta >> flag_gamma;
                if (infile.eof()) {
                    break;
                }

                add_planes_globals({0, to_kind(flag_x)},
                                   {0, to_kind(flag_y)},
                                   {0, to_kind(flag_z)},
                                   {0, to_kind(flag_alpha)},
                                   {0, to_kind(flag_beta)},
                                   {0, to_kind(flag_gamma)},
                                   translation_x,
                                   translation_y,
                                   translation_z,
                                   rotation_alpha * TMath::DegToRad(),
                                   rotation_beta * TMath::DegToRad(),
                                   rotation_gamma * TMath::DegToRad(),
                                   mb::Kind::FREE,
                                   mb::Kind::FREE,
                                   mb::Kind::FREE,
                                   mb::Kind::FREE);
            }
        } /*else if (config_word == "MATRIX") { FIXME need to specify exact rotation, like MATRIX-ZYZ, or similar
            double translation_x {}, translation_y {}, translation_z {};
            double rotation_11 {}, rotation_12 {}, rotation_13 {}, rotation_21 {}, rotation_22 {}, rotation_23 {}, rotation_31 {},
                rotation_32 {}, rotation_33 {};
            int flag_x {}, flag_y {}, flag_z {}, flag_alpha {}, flag_beta {}, flag_gamma {};
            while (true) {
                infile >> translation_x >> translation_y >> translation_z;
                infile >> rotation_11 >> rotation_12 >> rotation_13;
                infile >> rotation_21 >> rotation_22 >> rotation_23;
                infile >> rotation_31 >> rotation_32 >> rotation_33;
                infile >> flag_x >> flag_y >> flag_z >> flag_alpha >> flag_beta >> flag_gamma;
                if (infile.eof())
                    break;

                const auto rotation = Rotation3D(
                    rotation_11, rotation_12, rotation_13, rotation_21, rotation_22, rotation_23, rotation_31, rotation_32, rotation_33);
                const RotationZYX rzyx(rotation);  // FIXME wrong rotation

                add_planes_globals({0, to_kind(flag_x)},
                                   {0, to_kind(flag_y)},
                                   {0, to_kind(flag_z)},
                                   {0, to_kind(flag_alpha)},
                                   {0, to_kind(flag_beta)},
                                   {0, to_kind(flag_gamma)},
                                   translation_x,
                                   translation_y,
                                   translation_z,
                                   rzyx.Psi(),
                                   rzyx.Theta(),
                                   rzyx.Phi(),
                                   mb::Kind::FREE,
                                   mb::Kind::FREE,
                                   mb::Kind::FREE,
                                   mb::Kind::FREE);
            }
        }*/
        else
        {
            abort();
        }
    }

    auto print(bool use_deg) const -> void
    {
        constexpr auto width = 14;
        auto bar = std::string(width * 7, '=');
        std::cout << bar << '\n';
        std::cout << std::left << std::setw(width) << "Layer" << std::setw(width) << "Gx" << std::setw(width) << "Gy" << std::setw(width)
                  << "GzZ" << std::setw(width) << "Ga" << std::setw(width) << "Gb" << std::setw(width) << "Gc" << '\n';
        std::cout << bar << '\n';

        auto cnt = 1;
        for (const auto& global_par : layers_global_pars) {
            std::cout << std::setw(width) << cnt++ << global_par.gx() << global_par.gy() << global_par.gz() << global_par.ga()
                      << global_par.gb() << global_par.gc() << '\n';
        }

        std::cout << bar << '\n';

        cnt = 1;
        for (const auto& local_par : layers_local_pars) {
            std::cout << std::setw(width) << cnt++ << local_par.lx0() << local_par.ly0() << local_par.ltx() << local_par.lty() << '\n';
            ;
        }

        std::cout << bar << '\n';
        std::cout << std::right;
    }

    auto set_verbose(int make_verbose) -> void { verbose = make_verbose; }

  private:
    std::string mille_prefix;
    std::vector<global_parameters<float>> layers_global_pars;
    std::vector<derivatives<float, R>> layers_derivatives;
    std::vector<Rotation3D> alignment_rotation_corrections;
    std::vector<local_parameters> layers_local_pars;

    Mille mille;
    int verbose {0};
};

}  // namespace mb
