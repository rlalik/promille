#pragma once

#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include <Math/Point3Dfwd.h>
#include <Math/Rotation3D.h>
#include <Math/Vector3D.h>

#include "Mille.h"

/**
 * @brief Return the name of this header-only library
 */
inline auto name() -> std::string
{
    return "StrawAlignment";
}

/*
 * Rotation matrix, see https://mathworld.wolfram.com/EulerAngles.html (48..50)
 * for definitions.
 */

using ROOT::Math::Rotation3D;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

namespace SA
{

namespace helper
{
// clang-format off
template<typename T> constexpr auto R11(T psi, T theta, T phi) -> T { return cos(theta)*cos(psi); }
template<typename T> constexpr auto R12(T psi, T theta, T phi) -> T { return cos(theta)*sin(psi); }
template<typename T> constexpr auto R13(T psi, T theta, T phi) -> T { return -sin(theta); }

template<typename T> constexpr auto R21(T psi, T theta, T phi) -> T { return sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi); }
template<typename T> constexpr auto R22(T psi, T theta, T phi) -> T { return sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi); }
template<typename T> constexpr auto R23(T psi, T theta, T phi) -> T { return sin(phi)*cos(theta); }

template<typename T> constexpr auto R31(T psi, T theta, T phi) -> T { return cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi); }
template<typename T> constexpr auto R32(T psi, T theta, T phi) -> T { return cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi); }
template<typename T> constexpr auto R33(T psi, T theta, T phi) -> T { return cos(phi)*cos(theta); }
// clang-format on

template<typename T>
constexpr auto sign(T value) -> T
{
    return std::signbit(value) ? -1 : 1;
}
}  // namespace helper

template<typename T>
constexpr auto make_rotation_matrix(T psi, T theta, T phi) -> Rotation3D
{
    return Rotation3D(helper::R11(psi, theta, phi),
                      helper::R12(psi, theta, phi),
                      helper::R13(psi, theta, phi),
                      helper::R21(psi, theta, phi),
                      helper::R22(psi, theta, phi),
                      helper::R23(psi, theta, phi),
                      helper::R31(psi, theta, phi),
                      helper::R32(psi, theta, phi),
                      helper::R33(psi, theta, phi));
}

inline auto rotate(XYZVector v, Rotation3D R) -> XYZVector
{
    return R * v;
}

template<typename T>
constexpr auto make_point(T x0, T y0, T z0) -> XYZPoint
{
    return XYZPoint(x0, y0, z0);
}

template<typename T>
constexpr auto make_vector(T x0, T y0, T z0) -> XYZVector
{
    return XYZVector(x0, y0, z0);
}

inline auto distance(XYZPoint base1, XYZVector dir1, XYZPoint base2, XYZVector dir2) -> double
{
    return std::abs((base1 - base2).Dot(dir1.Cross(dir2)) / dir1.Cross(dir2).R());
}

template<typename T>
struct derivatives
{
    T cos_psi;
    T sin_psi;

    T cos_theta;
    T sin_theta;

    T cos_phi;
    T sin_phi;

    T sin_theta_cos_phi;
    T sin_theta_sin_phi;
    // T cos_phi_sin_phi;
    // T sin_theta_cos_theta;
    T cos_psi_cos_phi;
    // T cos_psi_cos_psi;
    T cos_psi_cos_theta;
    T cos_psi_sin_phi;
    T sin_psi_cos_phi;
    // T sin_psi_cos_psi;
    // T cos_psi_sin_psi;
    T sin_psi_cos_theta;
    T sin_psi_sin_phi;
    T cos_theta_cos_phi;
    // T cos_theta_cos_theta;
    T cos_theta_sin_phi;
    // T sin_psi_sin_psi;
    // T sin_theta_sin_theta;
    // T sin_phi_cos_phi;
    // T sin_phi_sin_phi;
    T cos_psi_sin_theta;
    T sin_psi_sin_theta;
    // T cos_theta_sin_theta;
    // T cos_phi_cos_phi;
    //
    // T sin_psi_sin_theta_sin_theta;
    // T cos_psi_sin_psi_cos_phi;
    // T cos_psi_sin_psi_cos_psi;
    // T cos_phi_cos_phi_sin_phi;
    // T cos_psi_sin_theta_sin_theta;
    T sin_psi_cos_theta_sin_phi;
    // T cos_psi_sin_phi_cos_phi;
    // T sin_theta_sin_theta_sin_theta;
    // T cos_theta_sin_phi_cos_phi;
    // T sin_phi_sin_phi_sin_phi;
    // T cos_psi_cos_psi_cos_phi;
    // T cos_theta_cos_theta_sin_theta;
    // T cos_psi_cos_psi_cos_psi;
    // T sin_theta_cos_theta_cos_theta;
    T cos_psi_sin_theta_sin_phi;
    // T cos_psi_cos_theta_cos_theta;
    // T sin_theta_sin_theta_sin_phi;
    // T sin_psi_sin_psi_cos_psi;
    // T cos_psi_cos_phi_cos_phi;
    // T sin_phi_cos_phi_sin_phi;
    // T cos_theta_cos_theta_sin_phi;
    // T sin_theta_cos_phi_cos_phi;
    // T cos_theta_sin_theta_sin_phi;
    T sin_psi_sin_theta_cos_phi;
    // T sin_phi_sin_phi_cos_phi;
    // T sin_psi_cos_theta_sin_theta;
    T cos_psi_sin_theta_cos_phi;
    // T sin_theta_sin_theta_cos_phi;
    // T cos_psi_cos_theta_sin_theta;
    // T sin_theta_sin_phi_cos_phi;
    // T sin_phi_cos_phi_cos_phi;
    // T cos_theta_cos_theta_cos_phi;
    // T cos_psi_sin_psi_cos_theta;
    // T sin_psi_cos_psi_cos_theta;
    // T cos_theta_sin_theta_sin_theta;
    T cos_psi_cos_theta_sin_phi;
    // T cos_phi_sin_phi_sin_phi;
    // T cos_phi_cos_phi_cos_phi;
    // T cos_psi_cos_psi_cos_theta;
    T sin_psi_cos_theta_cos_phi;
    // T sin_psi_sin_psi_cos_theta;
    T cos_psi_cos_theta_cos_phi;
    // T sin_psi_cos_psi_sin_theta;
    // T sin_theta_cos_theta_sin_theta;
    // T cos_theta_cos_phi_sin_phi;
    // T sin_psi_sin_theta_cos_theta;
    // T cos_phi_sin_phi_cos_phi;
    // T cos_theta_sin_theta_cos_phi;
    // T sin_psi_cos_psi_sin_phi;
    // T sin_theta_cos_theta_sin_phi;
    // T sin_psi_cos_psi_sin_psi;
    // T sin_psi_sin_psi_sin_theta;
    // T cos_theta_sin_phi_sin_phi;
    // T cos_psi_cos_psi_sin_phi;
    // T cos_psi_cos_psi_sin_psi;
    // T sin_psi_cos_phi_sin_phi;
    // T sin_psi_sin_psi_sin_phi;
    // T sin_psi_sin_psi_sin_psi;
    // T cos_psi_cos_phi_sin_phi;
    // T cos_psi_sin_psi_sin_theta;
    // T sin_psi_cos_psi_cos_phi;
    // T sin_psi_cos_psi_cos_psi;
    // T sin_theta_cos_theta_cos_phi;
    // T sin_theta_cos_phi_sin_phi;
    T sin_psi_sin_theta_sin_phi;
    // T sin_psi_cos_theta_cos_theta;
    // T sin_psi_sin_phi_sin_phi;
    // T sin_psi_cos_phi_cos_phi;
    // T cos_psi_sin_psi_sin_phi;
    // T cos_psi_sin_psi_sin_psi;
    // T cos_psi_cos_psi_sin_theta;
    // T cos_psi_sin_theta_cos_theta;
    // T sin_psi_sin_psi_cos_phi;
    // T sin_theta_sin_theta_cos_theta;
    // T sin_theta_sin_phi_sin_phi;
    // T cos_psi_sin_phi_sin_phi;
    // T cos_theta_cos_theta_cos_theta;
    // T cos_theta_sin_theta_cos_theta;
    // T cos_theta_cos_phi_cos_phi;
    // T sin_psi_sin_phi_cos_phi;

    T U;
    T V;
    T Z;
    T Xa;
    T Ya;
    T Za;
    T bx;
    T by;
    T bz;
    T tx;
    T ty;

    T corr_U;
    T corr_V;
    T corr_Z;

    T common_0;
    T common_1;
    T common_2;
    T common_3;
    T common_4;
    T common_5;
    T common_6;
    T common_7;
    T common_8;
    T common_9;

    T common_10;
    T common_11;
    T common_12;
    T common_13;
    T common_14;
    T common_16;
    T common_17;
    T common_18;
    T common_19;

    derivatives(T psi, T theta, T phi, T U, T V, T Z, T Xa, T Ya, T Za, T bx, T by, T bz, T tx, T ty)
        : U(U)
        , V(V)
        , Z(Z)
        , Xa(Xa)
        , Ya(Ya)
        , Za(Za)
        , bx(bx)
        , by(by)
        , bz(bz)
        , tx(tx)
        , ty(ty)
    {
        cos_psi = cos(psi);
        sin_psi = sin(psi);

        cos_theta = cos(theta);
        sin_theta = sin(theta);

        cos_phi = cos(phi);
        sin_phi = sin(phi);

        sin_theta_cos_phi = sin_theta * cos_phi;
        sin_theta_sin_phi = sin_theta * sin_phi;
        // cos_phi_sin_phi = cos_phi * sin_phi;
        // sin_theta_cos_theta = sin_theta * cos_theta;
        cos_psi_cos_phi = cos_psi * cos_phi;
        // cos_psi_cos_psi = cos_psi * cos_psi;
        cos_psi_cos_theta = cos_psi * cos_theta;
        cos_psi_sin_phi = cos_psi * sin_phi;
        sin_psi_cos_phi = sin_psi * cos_phi;
        // sin_psi_cos_psi = sin_psi * cos_psi;
        // cos_psi_sin_psi = cos_psi * sin_psi;
        sin_psi_cos_theta = sin_psi * cos_theta;
        sin_psi_sin_phi = sin_psi * sin_phi;
        cos_theta_cos_phi = cos_theta * cos_phi;
        // cos_theta_cos_theta = cos_theta * cos_theta;
        cos_theta_sin_phi = cos_theta * sin_phi;
        // sin_psi_sin_psi = sin_psi * sin_psi;
        // sin_theta_sin_theta = sin_theta * sin_theta;
        // sin_phi_cos_phi = sin_phi * cos_phi;
        // sin_phi_sin_phi = sin_phi * sin_phi;
        cos_psi_sin_theta = cos_psi * sin_theta;
        sin_psi_sin_theta = sin_psi * sin_theta;
        // cos_theta_sin_theta = cos_theta * sin_theta;
        // cos_phi_cos_phi = cos_phi * cos_phi;
        //
        // sin_psi_sin_theta_sin_theta = sin_psi * sin_theta * sin_theta;
        // cos_psi_sin_psi_cos_phi = cos_psi * sin_psi * cos_phi;
        // cos_psi_sin_psi_cos_psi = cos_psi * sin_psi * cos_psi;
        // cos_phi_cos_phi_sin_phi = cos_phi * cos_phi * sin_phi;
        // cos_psi_sin_theta_sin_theta = cos_psi * sin_theta * sin_theta;
        sin_psi_cos_theta_sin_phi = sin_psi * cos_theta * sin_phi;
        // cos_psi_sin_phi_cos_phi = cos_psi * sin_phi * cos_phi;
        // sin_theta_sin_theta_sin_theta = sin_theta * sin_theta * sin_theta;
        // cos_theta_sin_phi_cos_phi = cos_theta * sin_phi * cos_phi;
        // sin_phi_sin_phi_sin_phi = sin_phi * sin_phi * sin_phi;
        // cos_psi_cos_psi_cos_phi = cos_psi * cos_psi * cos_phi;
        // cos_theta_cos_theta_sin_theta = cos_theta * cos_theta * sin_theta;
        // cos_psi_cos_psi_cos_psi = cos_psi * cos_psi * cos_psi;
        // sin_theta_cos_theta_cos_theta = sin_theta * cos_theta * cos_theta;
        cos_psi_sin_theta_sin_phi = cos_psi * sin_theta * sin_phi;
        // cos_psi_cos_theta_cos_theta = cos_psi * cos_theta * cos_theta;
        // sin_theta_sin_theta_sin_phi = sin_theta * sin_theta * sin_phi;
        // sin_psi_sin_psi_cos_psi = sin_psi * sin_psi * cos_psi;
        // cos_psi_cos_phi_cos_phi = cos_psi * cos_phi * cos_phi;
        // sin_phi_cos_phi_sin_phi = sin_phi * cos_phi * sin_phi;
        // cos_theta_cos_theta_sin_phi = cos_theta * cos_theta * sin_phi;
        // sin_theta_cos_phi_cos_phi = sin_theta * cos_phi * cos_phi;
        // cos_theta_sin_theta_sin_phi = cos_theta * sin_theta * sin_phi;
        sin_psi_sin_theta_cos_phi = sin_psi * sin_theta * cos_phi;
        // sin_phi_sin_phi_cos_phi = sin_phi * sin_phi * cos_phi;
        // sin_psi_cos_theta_sin_theta = sin_psi * cos_theta * sin_theta;
        cos_psi_sin_theta_cos_phi = cos_psi * sin_theta * cos_phi;
        // sin_theta_sin_theta_cos_phi = sin_theta * sin_theta * cos_phi;
        // cos_psi_cos_theta_sin_theta = cos_psi * cos_theta * sin_theta;
        // sin_theta_sin_phi_cos_phi = sin_theta * sin_phi * cos_phi;
        // sin_phi_cos_phi_cos_phi = sin_phi * cos_phi * cos_phi;
        // cos_theta_cos_theta_cos_phi = cos_theta * cos_theta * cos_phi;
        // cos_psi_sin_psi_cos_theta = cos_psi * sin_psi * cos_theta;
        // sin_psi_cos_psi_cos_theta = sin_psi * cos_psi * cos_theta;
        // cos_theta_sin_theta_sin_theta = cos_theta * sin_theta * sin_theta;
        cos_psi_cos_theta_sin_phi = cos_psi * cos_theta * sin_phi;
        // cos_phi_sin_phi_sin_phi = cos_phi * sin_phi * sin_phi;
        // cos_phi_cos_phi_cos_phi = cos_phi * cos_phi * cos_phi;
        // cos_psi_cos_psi_cos_theta = cos_psi * cos_psi * cos_theta;
        sin_psi_cos_theta_cos_phi = sin_psi * cos_theta * cos_phi;
        // sin_psi_sin_psi_cos_theta = sin_psi * sin_psi * cos_theta;
        cos_psi_cos_theta_cos_phi = cos_psi * cos_theta * cos_phi;
        // sin_psi_cos_psi_sin_theta = sin_psi * cos_psi * sin_theta;
        // sin_theta_cos_theta_sin_theta = sin_theta * cos_theta * sin_theta;
        // cos_theta_cos_phi_sin_phi = cos_theta * cos_phi * sin_phi;
        // sin_psi_sin_theta_cos_theta = sin_psi * sin_theta * cos_theta;
        // cos_phi_sin_phi_cos_phi = cos_phi * sin_phi * cos_phi;
        // cos_theta_sin_theta_cos_phi = cos_theta * sin_theta * cos_phi;
        // sin_psi_cos_psi_sin_phi = sin_psi * cos_psi * sin_phi;
        // sin_theta_cos_theta_sin_phi = sin_theta * cos_theta * sin_phi;
        // sin_psi_cos_psi_sin_psi = sin_psi * cos_psi * sin_psi;
        // sin_psi_sin_psi_sin_theta = sin_psi * sin_psi * sin_theta;
        // cos_theta_sin_phi_sin_phi = cos_theta * sin_phi * sin_phi;
        // cos_psi_cos_psi_sin_phi = cos_psi * cos_psi * sin_phi;
        // cos_psi_cos_psi_sin_psi = cos_psi * cos_psi * sin_psi;
        // sin_psi_cos_phi_sin_phi = sin_psi * cos_phi * sin_phi;
        // sin_psi_sin_psi_sin_phi = sin_psi * sin_psi * sin_phi;
        // sin_psi_sin_psi_sin_psi = sin_psi * sin_psi * sin_psi;
        // cos_psi_cos_phi_sin_phi = cos_psi * cos_phi * sin_phi;
        // cos_psi_sin_psi_sin_theta = cos_psi * sin_psi * sin_theta;
        // sin_psi_cos_psi_cos_phi = sin_psi * cos_psi * cos_phi;
        // sin_psi_cos_psi_cos_psi = sin_psi * cos_psi * cos_psi;
        // sin_theta_cos_theta_cos_phi = sin_theta * cos_theta * cos_phi;
        // sin_theta_cos_phi_sin_phi = sin_theta * cos_phi * sin_phi;
        sin_psi_sin_theta_sin_phi = sin_psi * sin_theta * sin_phi;
        // sin_psi_cos_theta_cos_theta = sin_psi * cos_theta * cos_theta;
        // sin_psi_sin_phi_sin_phi = sin_psi * sin_phi * sin_phi;
        // sin_psi_cos_phi_cos_phi = sin_psi * cos_phi * cos_phi;
        // cos_psi_sin_psi_sin_phi = cos_psi * sin_psi * sin_phi;
        // cos_psi_sin_psi_sin_psi = cos_psi * sin_psi * sin_psi;
        // cos_psi_cos_psi_sin_theta = cos_psi * cos_psi * sin_theta;
        // cos_psi_sin_theta_cos_theta = cos_psi * sin_theta * cos_theta;
        // sin_psi_sin_psi_cos_phi = sin_psi * sin_psi * cos_phi;
        // sin_theta_sin_theta_cos_theta = sin_theta * sin_theta * cos_theta;
        // sin_theta_sin_phi_sin_phi = sin_theta * sin_phi * sin_phi;
        // cos_psi_sin_phi_sin_phi = cos_psi * sin_phi * sin_phi;
        // cos_theta_cos_theta_cos_theta = cos_theta * cos_theta * cos_theta;
        // cos_theta_sin_theta_cos_theta = cos_theta * sin_theta * cos_theta;
        // cos_theta_cos_phi_cos_phi = cos_theta * cos_phi * cos_phi;
        // sin_psi_sin_phi_cos_phi = sin_psi * sin_phi * cos_phi;

        // Manually optimized shortcuts for frequently appearing expressions.
        // Touch it on your own risk.
        corr_U = U + Xa;
        corr_V = V + Ya;
        corr_Z = Z + Za;

        common_0 = corr_Z * sin_psi_cos_theta - corr_U * (cos_psi_sin_phi + sin_psi_sin_theta_cos_phi)
            - corr_V * (cos_psi_cos_phi - sin_psi_sin_theta_sin_phi);
        common_1 = -(cos_theta_sin_phi)-tx * sin_psi_cos_phi - tx * cos_psi_sin_theta_sin_phi;
        common_2 = bz - corr_Z * cos_psi_cos_theta - corr_U * (sin_psi_sin_phi - cos_psi_sin_theta_cos_phi)
            - corr_V * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi);
        common_3 = tx * cos_psi_cos_phi + ty * cos_theta_sin_phi - tx * sin_psi_sin_theta_sin_phi;
        common_4 = -(cos_psi_cos_phi) + ty * sin_psi_cos_phi + ty * cos_psi_sin_theta_sin_phi + sin_psi_sin_theta_sin_phi;
        common_5 = bx - corr_U * cos_theta_cos_phi + corr_V * cos_theta_sin_phi - corr_Z * sin_theta;
        common_6 = pow(common_1, 2.0) + pow(common_4, 2.0) + pow(common_3, 2.0);
        common_7 = (2.0 * pow(common_6, 3.0 / 2.0));

        common_10 = common_5 * common_4 + common_3 * common_2 + common_1 * (by + common_0);
        common_11 = abs(common_10);
        common_12 = sqrt(common_6) * common_11;
    }

    // Manually optimized shortcuts for frequently appearing expressions.
    // Touch it on your own risk.

    constexpr auto dr_dpsi() -> T
    {
        return -(((-(tx * sin_psi_cos_phi) - tx * cos_psi_sin_theta_sin_phi) * common_3
                  + common_1 * (-(tx * cos_psi_cos_phi) + tx * sin_psi_sin_theta_sin_phi)
                  + common_4 * (ty * cos_psi_cos_phi + sin_psi_cos_phi + cos_psi_sin_theta_sin_phi - ty * sin_psi_sin_theta_sin_phi))
                 * common_11)
            / pow(common_6, 3.0 / 2.0)
            + (common_10
               * (common_5 * (ty * cos_psi_cos_phi + sin_psi_cos_phi + cos_psi_sin_theta_sin_phi - ty * sin_psi_sin_theta_sin_phi)
                  + common_1
                      * ((corr_Z)*cos_psi_cos_theta - corr_U * (-(sin_psi_sin_phi) + cos_psi_sin_theta_cos_phi)
                         - corr_V * (-(sin_psi_cos_phi)-cos_psi_sin_theta_sin_phi))
                  + (-(tx * sin_psi_cos_phi) - tx * cos_psi_sin_theta_sin_phi) * common_2 + common_3 * common_0
                  + (-(tx * cos_psi_cos_phi) + tx * sin_psi_sin_theta_sin_phi) * (by + common_0)))
            / common_12;
    }

    constexpr auto dr_dtheta() -> T
    {
        return (((ty * cos_psi_cos_theta_sin_phi + sin_psi_cos_theta_sin_phi) * common_5
                 + common_1 * (-corr_U * sin_psi_cos_theta_cos_phi + corr_V * sin_psi_cos_theta_sin_phi - (corr_Z)*sin_psi_sin_theta)
                 + (-corr_Z * cos_theta + corr_U * sin_theta_cos_phi - corr_V * sin_theta_sin_phi) * common_4
                 + (-(-corr_U * cos_psi_cos_theta_cos_phi) - corr_V * cos_psi_cos_theta_sin_phi + corr_Z * cos_psi_sin_theta) * common_3
                 + (-(tx * sin_psi_cos_theta_sin_phi) - ty * sin_theta_sin_phi) * common_2
                 + (-(tx * cos_psi_cos_theta_sin_phi) + sin_theta_sin_phi) * (by + common_0))
                * common_10)
            / common_12
            - ((2.0 * (-(tx * cos_psi_cos_theta_sin_phi) + sin_theta_sin_phi) * common_1
                + 2.0 * (ty * cos_psi_cos_theta_sin_phi + sin_psi_cos_theta_sin_phi) * common_4
                + 2.0 * (-(tx * sin_psi_cos_theta_sin_phi) - ty * sin_theta_sin_phi) * common_3)
               * common_11)
            / common_7;
    }

    constexpr auto dr_dphi() -> T
    {
        return ((common_5 * (cos_psi_sin_phi - ty * sin_psi_sin_phi + ty * cos_psi_sin_theta_cos_phi + sin_psi_sin_theta_cos_phi)
                 + (corr_V * cos_theta_cos_phi + corr_U * cos_theta_sin_phi) * common_4
                 + common_3
                     * (-corr_V * (-(sin_psi_sin_phi) + cos_psi_sin_theta_cos_phi) - corr_U * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi))
                 + (ty * cos_theta_cos_phi - tx * cos_psi_sin_phi - tx * sin_psi_sin_theta_cos_phi) * common_2
                 + common_1
                     * (-corr_V * (-(cos_psi_sin_phi)-sin_psi_sin_theta_cos_phi) - corr_U * (cos_psi_cos_phi - sin_psi_sin_theta_sin_phi))
                 + (-(cos_theta_cos_phi) + tx * sin_psi_sin_phi - tx * cos_psi_sin_theta_cos_phi) * (by + common_0))
                * common_10)
            / common_12
            - ((2.0 * (-(cos_theta_cos_phi) + tx * sin_psi_sin_phi - tx * cos_psi_sin_theta_cos_phi) * common_1
                + 2.0 * (cos_psi_sin_phi - ty * sin_psi_sin_phi + ty * cos_psi_sin_theta_cos_phi + sin_psi_sin_theta_cos_phi) * common_4
                + 2.0 * (ty * cos_theta_cos_phi - tx * cos_psi_sin_phi - tx * sin_psi_sin_theta_cos_phi) * common_3)
               * common_11)
            / common_7;
    }

    constexpr auto dr_dXa() -> T
    {
        return ((common_1 * (-(cos_psi_sin_phi)-sin_psi_sin_theta_cos_phi) - cos_theta_cos_phi * common_4
                 + (-(sin_psi_sin_phi) + cos_psi_sin_theta_cos_phi) * common_3)
                * common_10)
            / common_12;
    }

    constexpr auto dr_dYa() -> T
    {
        return ((common_1 * (-(cos_psi_cos_phi) + sin_psi_sin_theta_sin_phi) + cos_theta_sin_phi * common_4
                 + (-(sin_psi_cos_phi)-cos_psi_sin_theta_sin_phi) * common_3)
                * common_10)
            / common_12;
    }

    constexpr auto dr_dZa() -> T
    {
        return ((sin_psi_cos_theta * common_1 - sin_theta * common_4 - cos_psi_cos_theta * common_3) * common_10) / common_12;
    }

    constexpr auto dr_dbx() -> T { return (common_4 * common_10) / common_12; }

    constexpr auto dr_dby() -> T { return (common_1 * common_10) / common_12; }

    constexpr auto dr_dbz() -> T { return (common_3 * common_10) / common_12; }

    constexpr auto dr_dtx() -> T
    {
        return (((cos_psi_cos_phi - sin_psi_sin_theta_sin_phi) * common_2
                 + (-(sin_psi_cos_phi)-cos_psi_sin_theta_sin_phi) * (by + common_0))
                * common_10)
            / common_12
            - ((2.0 * (-(sin_psi_cos_phi)-cos_psi_sin_theta_sin_phi) * common_1
                + 2.0 * (cos_psi_cos_phi - sin_psi_sin_theta_sin_phi) * common_3)
               * common_11)
            / common_7;
    }

    constexpr auto dr_dty() -> T
    {
        return ((common_5 * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi) + cos_theta_sin_phi * common_2) * common_10) / common_12
            - ((2.0 * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi) * common_4 + 2.0 * cos_theta_sin_phi * common_3) * common_11)
            / common_7;
    }
};

enum class Kind
{
    FIXED,
    FREE
};

template<typename T>
struct parameter
{
    T value;
    Kind kind;

    operator bool() const { return kind == Kind::FREE; }
    explicit operator T() const { return value; }
};

template<typename T>
struct global_parameters
{
    global_parameters(parameter<T> Xa, parameter<T> Ya, parameter<T> Za, parameter<T> psi, parameter<T> phi, parameter<T> theta)
        : Xa(Xa)
        , Ya(Ya)
        , Za(Za)
        , psi(psi)
        , theta(theta)
        , phi(phi)
    {
        wire_rotation = make_rotation_matrix<T>(psi, theta, phi) * make_vector<T>(0, 1, 0);
    }

    parameter<T> Xa;
    parameter<T> Ya;
    parameter<T> Za;

    parameter<T> psi;
    parameter<T> phi;
    parameter<T> theta;

    XYZVector wire_rotation;
};

class MilleBuilder
{
  public:
    MilleBuilder(const char* outFileName, bool asBinary = true, bool writeZero = false)
        : mille(outFileName, asBinary, writeZero)
    {
    }
    ~MilleBuilder() {};

    auto add_planes_global(
        parameter<float> Xa, parameter<float> Ya, parameter<float> Za, parameter<float> psi, parameter<float> theta, parameter<float> phi)
    {
        global_parameters<float> gp = {Xa, Ya, Za, psi, theta, phi};
        layers_global_pars.push_back(std::move(gp));
    }

    /**
     * Add set of local variables for given layer of straws.
     *
     * @param layer layer number
     * @param x0 base vector x-coordinate
     * @param y0 base vector y-coordinate
     * @param z0 base vector z-coordinate
     * @param tx direction vector x-component
     * @param ty direction vector y-component
     * @param wx wire position x
     * @param wy wire position y
     * @param wz wire position z
     * @param rx wire rotation around x
     * @param ry wire rotation around y
     * @param rz wire rotation around z
     */
    auto add_local(int layer, float bx, float by, float bz, float tx, float ty, float U, float V, float Z, float dr, float sigma)
    {
        assert(layer < layers_global_pars.size());

        const auto param_idx_offset = layer * 10;
        const auto& current_layer_pars = layers_global_pars[layer];

        auto derivs = derivatives<float>(current_layer_pars.psi.value,
                                         current_layer_pars.theta.value,
                                         current_layer_pars.phi.value,
                                         U,
                                         V,
                                         Z,
                                         current_layer_pars.Xa.value,
                                         current_layer_pars.Ya.value,
                                         current_layer_pars.Za.value,
                                         bx,
                                         by,
                                         bz,
                                         tx,
                                         ty);

        std::vector<float> global_derivatives(10);
        std::vector<int> global_deriv_index(10);

        int global_params_count = 0;

        if (!current_layer_pars.Xa) {
            global_params_count++;

            global_derivatives.push_back(derivs.dr_dXa());
            global_deriv_index.push_back(param_idx_offset + 0);
        }

        if (!current_layer_pars.Ya) {
            global_params_count++;

            global_derivatives.push_back(derivs.dr_dYa());
            global_deriv_index.push_back(param_idx_offset + 1);
        }

        if (!current_layer_pars.Za) {
            global_params_count++;

            global_derivatives.push_back(derivs.dr_dZa());
            global_deriv_index.push_back(param_idx_offset + 2);
        }

        if (!current_layer_pars.phi) {
            global_params_count++;

            global_derivatives.push_back(derivs.dr_dphi());
            global_deriv_index.push_back(param_idx_offset + 3);
        }

        if (!current_layer_pars.theta) {
            global_params_count++;

            global_derivatives.push_back(derivs.dr_dtheta());
            global_deriv_index.push_back(param_idx_offset + 4);
        }

        if (!current_layer_pars.psi) {
            global_params_count++;

            global_derivatives.push_back(derivs.dr_dpsi());
            global_deriv_index.push_back(param_idx_offset + 5);
        }

        std::vector<float> local_derivatives(10);

        local_derivatives.push_back(derivs.dr_dbx());
        local_derivatives.push_back(derivs.dr_dby());
        local_derivatives.push_back(derivs.dr_dbz());
        local_derivatives.push_back(derivs.dr_dtx());
        local_derivatives.push_back(derivs.dr_dty());

        mille.mille(5,
                    local_derivatives.data(),
                    global_params_count,
                    global_derivatives.data(),
                    global_deriv_index.data(),
                    distance({bx, by, bz},
                             {tx, ty, 1.0},
                             {U + current_layer_pars.Xa, V + current_layer_pars.Ya, Z + current_layer_pars.Za},
                             current_layer_pars.wire_rotation)
                        - dr,
                    sigma);
    }

    /* Call Mille::end()
     */
    auto end() -> void { mille.end(); }

    /* Call Mille::kill();
     */
    auto kill() -> void { mille.kill(); }

  private:
    std::vector<global_parameters<float>> layers_global_pars;
    Mille mille;
};

}  // namespace SA
