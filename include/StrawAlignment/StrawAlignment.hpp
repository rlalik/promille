#pragma once

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

#include "Mille.h"

/*
 * Rotation matrix, see https://mathworld.wolfram.com/EulerAngles.html (48..50)
 * for definitions.
 */

using ROOT::Math::Rotation3D;
using ROOT::Math::RotationZYX;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

namespace SA
{

namespace geom
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

template<typename T>
constexpr auto make_rotation_matrix(T psi, T theta, T phi) -> Rotation3D
{
    return Rotation3D(R11(psi, theta, phi),
                      R12(psi, theta, phi),
                      R13(psi, theta, phi),
                      R21(psi, theta, phi),
                      R22(psi, theta, phi),
                      R23(psi, theta, phi),
                      R31(psi, theta, phi),
                      R32(psi, theta, phi),
                      R33(psi, theta, phi));
}

template<typename T>
constexpr auto make_rotation_matrix(T r11, T r12, T r13, T r21, T r22, T r23, T r31, T r32, T r33) -> Rotation3D
{
    return Rotation3D(r11, r12, r13, r21, r22, r23, r31, r32, r33);
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

}  // namespace geom

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

    T sx;
    T sy;
    T sz;

    T common_0;
    T common_1;
    T common_2;
    T common_3;
    T common_4;
    T common_5;
    T common_6;
    T common_7;

    T common_10;
    T common_11;
    T common_12;
    T common_13;

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
        sx = U + Xa;
        sy = V + Ya;
        sz = Z + Za;

#define VER2

#ifdef VER2
        common_0 = by + sz * sin_psi_cos_theta - sx * (cos_psi_sin_phi + sin_psi_sin_theta_cos_phi)
            - sy * (cos_psi_cos_phi - sin_psi_sin_theta_sin_phi);
        common_1 = -(cos_theta_sin_phi)-tx * sin_psi_cos_phi - tx * cos_psi_sin_theta_sin_phi;
        common_2 = bz - sz * cos_psi_cos_theta - sx * (sin_psi_sin_phi - cos_psi_sin_theta_cos_phi)
            - sy * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi);
        common_3 = tx * cos_psi_cos_phi + ty * cos_theta_sin_phi - tx * sin_psi_sin_theta_sin_phi;
        common_4 = -(cos_psi_cos_phi) + ty * sin_psi_cos_phi + ty * cos_psi_sin_theta_sin_phi + sin_psi_sin_theta_sin_phi;
        common_5 = bx - sx * cos_theta_cos_phi + sy * cos_theta_sin_phi - sz * sin_theta;
        common_6 = pow(common_1, 2.0) + pow(common_4, 2.0) + pow(common_3, 2.0);
        common_7 = 2.0 * pow(common_6, 3.0 / 2.0);

        common_10 = common_5 * common_4 + common_3 * common_2 + common_1 * common_0;
        common_11 = abs(common_10);
        common_12 = sqrt(common_6) * common_11;
        common_13 = common_10 / common_12;
#else
        common_0 = sz * sin_psi_cos_theta - sx * (cos_psi_sin_phi + sin_psi_sin_theta_cos_phi)
            - sy * (cos_psi_cos_phi - sin_psi_sin_theta_sin_phi);
        common_1 = -(cos_theta_sin_phi)-tx * sin_psi_cos_phi - tx * cos_psi_sin_theta_sin_phi;
        common_2 = bz - sz * cos_psi_cos_theta - sx * (sin_psi_sin_phi - cos_psi_sin_theta_cos_phi)
            - sy * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi);
        common_3 = tx * cos_psi_cos_phi + ty * cos_theta_sin_phi - tx * sin_psi_sin_theta_sin_phi;
        common_4 = -(cos_psi_cos_phi) + ty * sin_psi_cos_phi + ty * cos_psi_sin_theta_sin_phi + sin_psi_sin_theta_sin_phi;
        common_5 = bx - sx * cos_theta_cos_phi + sy * cos_theta_sin_phi - sz * sin_theta;
        common_6 = pow(common_1, 2.0) + pow(common_4, 2.0) + pow(common_3, 2.0);
        common_7 = 2.0 * pow(common_6, 3.0 / 2.0);

        common_10 = common_5 * common_4 + common_3 * common_2 + common_1 * (by + common_0);
        common_11 = abscommon_10;
        common_12 = sqrt(common_6) * common_11;
#endif
    }

    // Manually optimized shortcuts for frequently appearing expressions.
    // Touch it on your own risk.

#ifdef VER2
    // #    include "StrawAlignment/functions_bodies.h"

    constexpr auto dr_dpsi() -> T
    {
        return -1 / 2.0
            * ((2.0 * (-(tx * sin_psi_cos_phi) - tx * cos_psi_sin_theta_sin_phi) * common_3
                + 2.0 * common_1 * (-(tx * cos_psi_cos_phi) + tx * sin_psi_sin_theta_sin_phi)
                + 2.0 * common_4 * (ty * cos_psi_cos_phi + sin_psi_cos_phi + cos_psi_sin_theta_sin_phi - ty * sin_psi_sin_theta_sin_phi))
               * common_11)
            / pow(common_6, 3.0 / 2.0)
            + (common_10
               * (common_5 * (ty * cos_psi_cos_phi + sin_psi_cos_phi + cos_psi_sin_theta_sin_phi - ty * sin_psi_sin_theta_sin_phi)
                  + common_1
                      * (sz * cos_psi_cos_theta + sy * sin_psi_cos_phi + sy * cos_psi_sin_theta_sin_phi
                         - sx * (-(sin_psi_sin_phi) + cos_psi_sin_theta_cos_phi))
                  + (-(tx * sin_psi_cos_phi) - tx * cos_psi_sin_theta_sin_phi) * common_2
                  + (-(tx * cos_psi_cos_phi) + tx * sin_psi_sin_theta_sin_phi) * common_0
                  + common_3
                      * (-sx * (cos_psi_sin_phi + sin_psi_sin_theta_cos_phi) - sy * (cos_psi_cos_phi - sin_psi_sin_theta_sin_phi)
                         + sz * sin_psi_cos_theta)))
            / common_12;
    }

    constexpr auto dr_dtheta() -> T
    {
        return ((ty * cos_psi_cos_theta_sin_phi + sin_psi_cos_theta_sin_phi) * common_5
                + common_1 * (-(sx * sin_psi_cos_theta_cos_phi) + sy * sin_psi_cos_theta_sin_phi - sz * sin_psi_sin_theta)
                + (-(sz * cos_theta) + sx * sin_theta_cos_phi - sy * sin_theta_sin_phi) * common_4
                + (sx * cos_psi_cos_theta_cos_phi - sy * cos_psi_cos_theta_sin_phi + sz * cos_psi_sin_theta) * common_3
                + (-(tx * sin_psi_cos_theta_sin_phi) - ty * sin_theta_sin_phi) * common_2
                + (-(tx * cos_psi_cos_theta_sin_phi) + sin_theta_sin_phi) * common_0)
            * common_13
            - ((2.0 * (-(tx * cos_psi_cos_theta_sin_phi) + sin_theta_sin_phi) * common_1
                + 2.0 * (ty * cos_psi_cos_theta_sin_phi + sin_psi_cos_theta_sin_phi) * common_4
                + 2.0 * (-(tx * sin_psi_cos_theta_sin_phi) - ty * sin_theta_sin_phi) * common_3)
               * common_11)
            / common_7;
    }

    constexpr auto dr_dphi() -> T
    {
        return -1 / 2.0
            * ((2.0 * (-(cos_theta_cos_phi) + tx * sin_psi_sin_phi - tx * cos_psi_sin_theta_cos_phi) * common_1
                + 2.0 * (cos_psi_sin_phi - ty * sin_psi_sin_phi + ty * cos_psi_sin_theta_cos_phi + sin_psi_sin_theta_cos_phi) * common_4
                + 2.0 * (ty * cos_theta_cos_phi - tx * cos_psi_sin_phi - tx * sin_psi_sin_theta_cos_phi) * common_3)
               * common_11)
            / pow(common_6, 3.0 / 2.0)
            + (common_10
               * (common_5 * (cos_psi_sin_phi - ty * sin_psi_sin_phi + ty * cos_psi_sin_theta_cos_phi + sin_psi_sin_theta_cos_phi)
                  + (sy * cos_theta_cos_phi + sx * cos_theta_sin_phi) * common_4
                  + common_3 * (-sx * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi) - sy * (-(sin_psi_sin_phi) + cos_psi_sin_theta_cos_phi))
                  + (ty * cos_theta_cos_phi - tx * cos_psi_sin_phi - tx * sin_psi_sin_theta_cos_phi) * common_2
                  + (-(cos_theta_cos_phi) + tx * sin_psi_sin_phi - tx * cos_psi_sin_theta_cos_phi) * common_0
                  + common_1
                      * (sy * cos_psi_sin_phi + sy * sin_psi_sin_theta_cos_phi - sx * (cos_psi_cos_phi - sin_psi_sin_theta_sin_phi))))
            / common_12;
    }

    constexpr auto dr_dXa() -> T
    {
        return (common_1 * (-(cos_psi_sin_phi)-sin_psi_sin_theta_cos_phi) - cos_theta_cos_phi * common_4
                + (-(sin_psi_sin_phi) + cos_psi_sin_theta_cos_phi) * common_3)
            * common_13;
    }

    constexpr auto dr_dYa() -> T
    {
        return (common_1 * (-(cos_psi_cos_phi) + sin_psi_sin_theta_sin_phi) + cos_theta_sin_phi * common_4
                + (-(sin_psi_cos_phi)-cos_psi_sin_theta_sin_phi) * common_3)
            * common_13;
    }

    constexpr auto dr_dZa() -> T
    {
        return (sin_psi_cos_theta * common_1 - sin_theta * common_4 - cos_psi_cos_theta * common_3) * common_13;
    }

    constexpr auto dr_dbx() -> T { return common_4 * common_13; }

    constexpr auto dr_dby() -> T { return common_1 * common_13; }

    constexpr auto dr_dbz() -> T { return common_3 * common_13; }

    constexpr auto dr_dtx() -> T
    {
        return ((cos_psi_cos_phi - sin_psi_sin_theta_sin_phi) * common_2 + (-(sin_psi_cos_phi)-cos_psi_sin_theta_sin_phi) * common_0)
            * common_13
            - ((2.0 * (-(sin_psi_cos_phi)-cos_psi_sin_theta_sin_phi) * common_1
                + 2.0 * (cos_psi_cos_phi - sin_psi_sin_theta_sin_phi) * common_3)
               * common_11)
            / common_7;
    }

    constexpr auto dr_dty() -> T
    {
        return (common_5 * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi) + cos_theta_sin_phi * common_2) * common_13
            - ((2.0 * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi) * common_4 + 2.0 * cos_theta_sin_phi * common_3) * common_11)
            / common_7;
    }

#else
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
                      * (sz * cos_psi_cos_theta - sx * (-(sin_psi_sin_phi) + cos_psi_sin_theta_cos_phi)
                         - sy * (-(sin_psi_cos_phi)-cos_psi_sin_theta_sin_phi))
                  + (-(tx * sin_psi_cos_phi) - tx * cos_psi_sin_theta_sin_phi) * common_2 + common_3 * common_0
                  + (-(tx * cos_psi_cos_phi) + tx * sin_psi_sin_theta_sin_phi) * (by + common_0)))
            / common_12;
    }

    constexpr auto dr_dtheta() -> T
    {
        return (((ty * cos_psi_cos_theta_sin_phi + sin_psi_cos_theta_sin_phi) * common_5
                 + common_1 * (-sx * sin_psi_cos_theta_cos_phi + sy * sin_psi_cos_theta_sin_phi - (sz)*sin_psi_sin_theta)
                 + (-sz * cos_theta + sx * sin_theta_cos_phi - sy * sin_theta_sin_phi) * common_4
                 + (-(-sx * cos_psi_cos_theta_cos_phi) - sy * cos_psi_cos_theta_sin_phi + sz * cos_psi_sin_theta) * common_3
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
                 + (sy * cos_theta_cos_phi + sx * cos_theta_sin_phi) * common_4
                 + common_3 * (-sy * (-(sin_psi_sin_phi) + cos_psi_sin_theta_cos_phi) - sx * (sin_psi_cos_phi + cos_psi_sin_theta_sin_phi))
                 + (ty * cos_theta_cos_phi - tx * cos_psi_sin_phi - tx * sin_psi_sin_theta_cos_phi) * common_2
                 + common_1 * (-sy * (-(cos_psi_sin_phi)-sin_psi_sin_theta_cos_phi) - sx * (cos_psi_cos_phi - sin_psi_sin_theta_sin_phi))
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
#endif
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

template<typename T>
struct parameter
{
    T value;
    Kind kind;

    explicit operator bool() const { return kind == Kind::FREE; }
    operator T() const { return value; }

    auto operator<<(std::ofstream& ofs) { ofs << value << ' ' << (kind == Kind::FREE ? '~' : ' '); }
};

template<typename T>
struct global_parameters
{
    global_parameters(parameter<T> Xa, parameter<T> Ya, parameter<T> Za, parameter<T> psi, parameter<T> theta, parameter<T> phi)
        : Xa(Xa)
        , Ya(Ya)
        , Za(Za)
        , psi(psi)
        , theta(theta)
        , phi(phi)
    {
        wire_rotation = geom::make_rotation_matrix<T>(psi, theta, phi) * geom::make_vector<T>(0, 1, 0);
    }

    parameter<T> Xa;
    parameter<T> Ya;
    parameter<T> Za;

    parameter<T> psi;
    parameter<T> theta;
    parameter<T> phi;

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

        std::vector<float> local_derivatives;
        local_derivatives.reserve(5);

        local_derivatives.push_back(derivs.dr_dbx());
        local_derivatives.push_back(derivs.dr_dby());
        local_derivatives.push_back(derivs.dr_dbz());
        local_derivatives.push_back(derivs.dr_dtx());
        local_derivatives.push_back(derivs.dr_dty());

        std::vector<float> global_derivatives;
        global_derivatives.reserve(6);
        std::vector<int> global_deriv_index;
        global_deriv_index.reserve(6);

        if (current_layer_pars.Xa) {
            global_derivatives.push_back(derivs.dr_dXa());
            global_deriv_index.push_back(param_idx_offset + 1);
        }

        if (current_layer_pars.Ya) {
            global_derivatives.push_back(derivs.dr_dYa());
            global_deriv_index.push_back(param_idx_offset + 2);
        }

        if (current_layer_pars.Za) {
            global_derivatives.push_back(derivs.dr_dZa());
            global_deriv_index.push_back(param_idx_offset + 3);
        }

        if (current_layer_pars.phi) {
            global_derivatives.push_back(derivs.dr_dphi());
            global_deriv_index.push_back(param_idx_offset + 4);
        }

        if (current_layer_pars.theta) {
            global_derivatives.push_back(derivs.dr_dtheta());
            global_deriv_index.push_back(param_idx_offset + 5);
        }

        if (current_layer_pars.psi) {
            global_derivatives.push_back(derivs.dr_dpsi());
            global_deriv_index.push_back(param_idx_offset + 6);
        }

        auto res = geom::distance({bx, by, bz},
                                  {tx, ty, 1.0},
                                  {U + current_layer_pars.Xa, V + current_layer_pars.Ya, Z + current_layer_pars.Za},
                                  current_layer_pars.wire_rotation)
            - dr;

        if (verbose) {
            static const auto value_width = 12;
            std::cout << "+ " << std::setw(3) << layer << "   res=" << std::setw(value_width) << res
                      << "   sigma=" << std::setw(value_width) << sigma << "   NLC=" << local_derivatives.size();
            for (const auto& ld : local_derivatives)
                std::cout << "  " << std::setw(value_width) << ld;

            std::cout << "   NGL=" << global_derivatives.size();
            for (const auto& gd : global_derivatives)
                std::cout << "  " << std::setw(value_width) << gd;
            std::cout << '\n';
        }

        mille.mille(local_derivatives.size(),
                    local_derivatives.data(),
                    global_derivatives.size(),
                    global_derivatives.data(),
                    global_deriv_index.data(),
                    res,
                    sigma);
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
            double x, y, z, psi, theta, phi;
            int bx, by, bz, bpsi, btheta, bphi;
            while (true) {
                infile >> x >> y >> z >> psi >> theta >> phi;
                infile >> bx >> by >> bz >> bpsi >> btheta >> bphi;
                if (infile.eof())
                    break;

                add_planes_global({x, to_kind(bx)},
                                  {y, to_kind(by)},
                                  {z, to_kind(bz)},
                                  {psi * TMath::DegToRad(), to_kind(bpsi)},
                                  {theta * TMath::DegToRad(), to_kind(btheta)},
                                  {phi * TMath::DegToRad(), to_kind(bpsi)});
            }
        } else if (config_word == "MATRIX") {
            double x, y, z, r11, r12, r13, r21, r22, r23, r31, r32, r33;
            int bx, by, bz, bpsi, btheta, bphi;
            while (true) {
                infile >> x >> y >> z >> r11 >> r12 >> r13 >> r21 >> r22 >> r23 >> r31 >> r32 >> r33;
                infile >> bx >> by >> bz >> bpsi >> btheta >> bphi;
                if (infile.eof())
                    break;

                auto r = geom::make_rotation_matrix(r11, r12, r13, r21, r22, r23, r31, r32, r33);
                RotationZYX ra(r);

                add_planes_global({x, to_kind(bx)},
                                  {y, to_kind(by)},
                                  {z, to_kind(bz)},
                                  {ra.Psi(), to_kind(bpsi)},
                                  {ra.Theta(), to_kind(btheta)},
                                  {ra.Phi(), to_kind(bpsi)});
            }
        } else {
            abort();
        }
    }

    auto print(bool use_deg) const -> void
    {
        auto rad_deg = use_deg ? TMath::RadToDeg() : 1.0;

        constexpr auto width = 14;
        auto bar = std::string(width * 7, '=');
        std::cout << bar << '\n';
        std::cout << std::left << std::setw(width) << "Layer" << std::setw(width) << "X" << std::setw(width) << "Y" << std::setw(width)
                  << "Z" << std::setw(width) << "Psi" << std::setw(width) << "Theta" << std::setw(width) << "Phi" << '\n';
        std::cout << bar << '\n';

        auto i = 1;
        for (const auto& gl : layers_global_pars) {
            std::cout << std::setw(width) << i++ << std::setw(width) << gl.Xa << std::setw(width) << gl.Ya << std::setw(width) << gl.Za
                      << std::setw(width) << gl.psi * rad_deg << std::setw(width) << gl.theta * rad_deg << std::setw(width)
                      << gl.phi * rad_deg << '\n';
        }

        std::cout << bar << '\n';
        std::cout << std::right;
    }

    auto set_verbose(bool make_verbose) { verbose = make_verbose; }

  private:
    std::vector<global_parameters<float>> layers_global_pars;
    Mille mille;
    bool verbose {false};
};

}  // namespace SA
