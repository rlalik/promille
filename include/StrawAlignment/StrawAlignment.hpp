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
#include <StrawAlignment/EulerAngles.hpp>
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

template<typename T>
constexpr auto sign(T value) -> T
{
    return std::signbit(value) ? -1 : 1;
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

template<typename T, template<class> class R>
struct derivatives
{
    // gloal translational corrections
    T gx {0};
    T gy {0};
    T gz {0};

    // gloal rotational corrections
    T ga {0};
    T gb {0};
    T gc {0};

    // translational alignment of straws
    T ax {0};
    T ay {0};
    T az {0};

    // local system straw coordinates
    T sx {0};
    T sy {0};
    T sz {0};

    // track base
    T bx {0};
    T by {0};
    T bz {0};
    // track direction
    T tx {0};
    T ty {0};

    T common_0;
    T common_1;
    T common_2;
    T common_3;
    T common_4;

    SA::euler::euler_base<T> wm;

    derivatives(T gx, T gy, T gz, T ga, T gb, T gc, T ax, T ay, T az, T alpha, T beta, T gamma)
        : gx(gx)
        , gy(gy)
        , gz(gz)
        , ga(ga)
        , gb(gb)
        , gc(gc)
        , ax(ax)
        , ay(ay)
        , az(az)
        , wm(R<T>(alpha, beta, gamma))
    {
        // std::cout << "Gt" << ROOT::Math::XYZVector(gx, gy, gz) << "\nGr" << ROOT::Math::XYZVector(ga, gb, gc) << "\nAt"
        // << ROOT::Math::XYZVector(ax, ay, az) << "\nWM: " << euler::make_rotation_matrix(wm) << '\n';
    }

    auto update(T sx_, T sy_, T sz_, T bx_, T by_, T bz_, T tx_, T ty_) -> void
    {
        sx = sx_;
        sy = sy_;
        sz = sz_;
        bx = bx_;
        by = by_;
        bz = bz_;
        tx = tx_;
        ty = ty_;

        // std::cout << "Sl" << ROOT::Math::XYZVector(sx, sy, sz) << "\nBt" << ROOT::Math::XYZVector(bx, by, bz) << "\nTt"
        // << ROOT::Math::XYZVector(tx, ty, 1) << '\n';

        // Manually optimized shortcuts for frequently appearing expressions.
        // Touch it on your own risk.
        /*
                common_0 = abs((-ay + by - sx * wm.R21 - sy * wm.R22 - sz * wm.R23) * (wm.R12 - tx * wm.R32)
                               + (-ax + bx - sx * wm.R11 - sy * wm.R12 - sz * wm.R13) * (-wm.R22 + ty * wm.R32)
                               + (-(ty * wm.R12) + tx * wm.R22) * (-az + bz - sx * wm.R31 - sy * wm.R32 - sz * wm.R33));
                common_1 = sqrt(pow(-(ty * wm.R12) + tx * wm.R22, 2) + pow(wm.R12 - tx * wm.R32, 2) + pow(-wm.R22 + ty * wm.R32, 2));
                common_2 = pow(pow(-(ty * wm.R12) + tx * wm.R22, 2) + pow(wm.R12 - tx * wm.R32, 2) + pow(-wm.R22 + ty * wm.R32, 2), 3.0
           / 2.0); common_3 = ((-ay + by - sx * wm.R21 - sy * wm.R22 - sz * wm.R23) * (wm.R12 - tx * wm.R32)
                            + (-ax + bx - sx * wm.R11 - sy * wm.R12 - sz * wm.R13) * (-wm.R22 + ty * wm.R32)
                            + (-(ty * wm.R12) + tx * wm.R22) * (-az + bz - sx * wm.R31 - sy * wm.R32 - sz * wm.R33));*/

        common_0 = pow(wm.R12 - gb * tx * wm.R12 - gc * wm.R22 + ga * tx * wm.R22 + gb * wm.R32 - tx * wm.R32, 2)
            + pow(gc * wm.R12 + gb * ty * wm.R12 - wm.R22 - ga * ty * wm.R22 - ga * wm.R32 + ty * wm.R32, 2)
            + pow(-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32, 2);
        common_1 = pow(common_0, 3. / 2.);
        common_2 = ((-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32)
                        * (-2 * az + bz - sx * (gb * wm.R11 - ga * wm.R21 + wm.R31) - sy * (gb * wm.R12 - ga * wm.R22 + wm.R32)
                           - sz * (gb * wm.R13 - ga * wm.R23 + wm.R33))
                    + (wm.R12 - gb * tx * wm.R12 - gc * wm.R22 + ga * tx * wm.R22 + gb * wm.R32 - tx * wm.R32)
                        * (-2 * ay + by - sx * (-(gc * wm.R11) + wm.R21 + ga * wm.R31) - sy * (-(gc * wm.R12) + wm.R22 + ga * wm.R32)
                           - sz * (-(gc * wm.R13) + wm.R23 + ga * wm.R33))
                    + (gc * wm.R12 + gb * ty * wm.R12 - wm.R22 - ga * ty * wm.R22 - ga * wm.R32 + ty * wm.R32)
                        * (-2 * ax + bx - sx * (wm.R11 - gc * wm.R21 + gb * wm.R31) - sy * (wm.R12 - gc * wm.R22 + gb * wm.R32)
                           - sz * (wm.R13 - gc * wm.R23 + gb * wm.R33)));
        common_3 = fabs(common_2);
        common_4 = sqrt(common_0);
    }

    // Manually optimized shortcuts for frequently appearing expressions.
    // Touch it on your own risk.

    constexpr auto dr_dgx() -> T { return 0; }

    constexpr auto dr_dgy() -> T { return 0; }

    constexpr auto dr_dgz() -> T { return 0; }

    constexpr auto dr_dga() -> T
    {
        return (((sx * wm.R21 + sy * wm.R22 + sz * wm.R23)
                     * (-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32)
                 + (wm.R12 - gb * tx * wm.R12 - gc * wm.R22 + ga * tx * wm.R22 + gb * wm.R32 - tx * wm.R32)
                     * (-(sx * wm.R31) - sy * wm.R32 - sz * wm.R33)
                 + tx * wm.R32
                     * (-2 * az + bz - sx * (gb * wm.R11 - ga * wm.R21 + wm.R31) - sy * (gb * wm.R12 - ga * wm.R22 + wm.R32)
                        - sz * (gb * wm.R13 - ga * wm.R23 + wm.R33))
                 + tx * wm.R22
                     * (-2 * ay + by - sx * (-(gc * wm.R11) + wm.R21 + ga * wm.R31) - sy * (-(gc * wm.R12) + wm.R22 + ga * wm.R32)
                        - sz * (-(gc * wm.R13) + wm.R23 + ga * wm.R33))
                 + (-(ty * wm.R22) - wm.R32)
                     * (-2 * ax + bx - sx * (wm.R11 - gc * wm.R21 + gb * wm.R31) - sy * (wm.R12 - gc * wm.R22 + gb * wm.R32)
                        - sz * (wm.R13 - gc * wm.R23 + gb * wm.R33)))
                * common_2)
            / (common_4 * common_3)
            - ((2 * tx * wm.R22 * (wm.R12 - gb * tx * wm.R12 - gc * wm.R22 + ga * tx * wm.R22 + gb * wm.R32 - tx * wm.R32)
                + 2 * (-(ty * wm.R22) - wm.R32) * (gc * wm.R12 + gb * ty * wm.R12 - wm.R22 - ga * ty * wm.R22 - ga * wm.R32 + ty * wm.R32)
                + 2 * tx * wm.R32
                    * (-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32))
               * common_3)
            / (2 * common_1);
    }

    constexpr auto dr_dgb() -> T
    {
        return (((-(sx * wm.R11) - sy * wm.R12 - sz * wm.R13)
                     * (-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32)
                 + (gc * wm.R12 + gb * ty * wm.R12 - wm.R22 - ga * ty * wm.R22 - ga * wm.R32 + ty * wm.R32)
                     * (-(sx * wm.R31) - sy * wm.R32 - sz * wm.R33)
                 - ty * wm.R32
                     * (-2 * az + bz - sx * (gb * wm.R11 - ga * wm.R21 + wm.R31) - sy * (gb * wm.R12 - ga * wm.R22 + wm.R32)
                        - sz * (gb * wm.R13 - ga * wm.R23 + wm.R33))
                 + (-(tx * wm.R12) + wm.R32)
                     * (-2 * ay + by - sx * (-(gc * wm.R11) + wm.R21 + ga * wm.R31) - sy * (-(gc * wm.R12) + wm.R22 + ga * wm.R32)
                        - sz * (-(gc * wm.R13) + wm.R23 + ga * wm.R33))
                 + ty * wm.R12
                     * (-2 * ax + bx - sx * (wm.R11 - gc * wm.R21 + gb * wm.R31) - sy * (wm.R12 - gc * wm.R22 + gb * wm.R32)
                        - sz * (wm.R13 - gc * wm.R23 + gb * wm.R33)))
                * common_2)
            / (common_4 * common_3)
            - ((2 * (-(tx * wm.R12) + wm.R32) * (wm.R12 - gb * tx * wm.R12 - gc * wm.R22 + ga * tx * wm.R22 + gb * wm.R32 - tx * wm.R32)
                + 2 * ty * wm.R12 * (gc * wm.R12 + gb * ty * wm.R12 - wm.R22 - ga * ty * wm.R22 - ga * wm.R32 + ty * wm.R32)
                - 2 * ty * wm.R32
                    * (-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32))
               * common_3)
            / (2 * common_1);
    }

    constexpr auto dr_dgc() -> T
    {
        return (((sx * wm.R11 + sy * wm.R12 + sz * wm.R13)
                     * (wm.R12 - gb * tx * wm.R12 - gc * wm.R22 + ga * tx * wm.R22 + gb * wm.R32 - tx * wm.R32)
                 + (sx * wm.R21 + sy * wm.R22 + sz * wm.R23)
                     * (gc * wm.R12 + gb * ty * wm.R12 - wm.R22 - ga * ty * wm.R22 - ga * wm.R32 + ty * wm.R32)
                 + (-(tx * wm.R12) + ty * wm.R22)
                     * (-2 * az + bz - sx * (gb * wm.R11 - ga * wm.R21 + wm.R31) - sy * (gb * wm.R12 - ga * wm.R22 + wm.R32)
                        - sz * (gb * wm.R13 - ga * wm.R23 + wm.R33))
                 - wm.R22
                     * (-2 * ay + by - sx * (-(gc * wm.R11) + wm.R21 + ga * wm.R31) - sy * (-(gc * wm.R12) + wm.R22 + ga * wm.R32)
                        - sz * (-(gc * wm.R13) + wm.R23 + ga * wm.R33))
                 + wm.R12
                     * (-2 * ax + bx - sx * (wm.R11 - gc * wm.R21 + gb * wm.R31) - sy * (wm.R12 - gc * wm.R22 + gb * wm.R32)
                        - sz * (wm.R13 - gc * wm.R23 + gb * wm.R33)))
                * common_2)
            / (common_4 * common_3)
            - ((-2 * wm.R22 * (wm.R12 - gb * tx * wm.R12 - gc * wm.R22 + ga * tx * wm.R22 + gb * wm.R32 - tx * wm.R32)
                + 2 * wm.R12 * (gc * wm.R12 + gb * ty * wm.R12 - wm.R22 - ga * ty * wm.R22 - ga * wm.R32 + ty * wm.R32)
                + 2 * (-(tx * wm.R12) + ty * wm.R22)
                    * (-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32))
               * common_3)
            / (2 * common_1);
    }

    constexpr auto dr_dbx() -> T
    {
        return ((gc * wm.R12 + gb * ty * wm.R12 - wm.R22 - ga * ty * wm.R22 - ga * wm.R32 + ty * wm.R32) * common_2)
            / (common_4 * common_3);
    }

    constexpr auto dr_dby() -> T
    {
        return ((wm.R12 - gb * tx * wm.R12 - gc * wm.R22 + ga * tx * wm.R22 + gb * wm.R32 - tx * wm.R32) * common_2)
            / (common_4 * common_3);
    }

    constexpr auto dr_dbz() -> T
    {
        return ((-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32) * common_2)
            / (common_4 * common_3);
    }

    constexpr auto dr_dtx() -> T
    {
        return (((-(gc * wm.R12) + wm.R22 + ga * wm.R32)
                     * (-2 * az + bz - sx * (gb * wm.R11 - ga * wm.R21 + wm.R31) - sy * (gb * wm.R12 - ga * wm.R22 + wm.R32)
                        - sz * (gb * wm.R13 - ga * wm.R23 + wm.R33))
                 + (-(gb * wm.R12) + ga * wm.R22 - wm.R32)
                     * (-2 * ay + by - sx * (-(gc * wm.R11) + wm.R21 + ga * wm.R31) - sy * (-(gc * wm.R12) + wm.R22 + ga * wm.R32)
                        - sz * (-(gc * wm.R13) + wm.R23 + ga * wm.R33)))
                * common_2)
            / (common_4 * common_3)
            - ((2 * (-(gb * wm.R12) + ga * wm.R22 - wm.R32)
                    * (wm.R12 - gb * tx * wm.R12 - gc * wm.R22 + ga * tx * wm.R22 + gb * wm.R32 - tx * wm.R32)
                + 2 * (-(gc * wm.R12) + wm.R22 + ga * wm.R32)
                    * (-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32))
               * common_3)
            / (2 * common_1);
    }

    constexpr auto dr_dty() -> T
    {
        return (((-wm.R12 + gc * wm.R22 - gb * wm.R32)
                     * (-2 * az + bz - sx * (gb * wm.R11 - ga * wm.R21 + wm.R31) - sy * (gb * wm.R12 - ga * wm.R22 + wm.R32)
                        - sz * (gb * wm.R13 - ga * wm.R23 + wm.R33))
                 + (gb * wm.R12 - ga * wm.R22 + wm.R32)
                     * (-2 * ax + bx - sx * (wm.R11 - gc * wm.R21 + gb * wm.R31) - sy * (wm.R12 - gc * wm.R22 + gb * wm.R32)
                        - sz * (wm.R13 - gc * wm.R23 + gb * wm.R33)))
                * common_2)
            / (common_4 * common_3)
            - ((2 * (gb * wm.R12 - ga * wm.R22 + wm.R32)
                    * (gc * wm.R12 + gb * ty * wm.R12 - wm.R22 - ga * ty * wm.R22 - ga * wm.R32 + ty * wm.R32)
                + 2 * (-wm.R12 + gc * wm.R22 - gb * wm.R32)
                    * (-(gc * tx * wm.R12) - ty * wm.R12 + tx * wm.R22 + gc * ty * wm.R22 + ga * tx * wm.R32 - gb * ty * wm.R32))
               * common_3)
            / (2 * common_1);
    }

    // constexpr auto dr_dpsi() -> T
    // {
    //     return -0.5
    //         * (common_0
    //            * (2 * (-(ty * wm.R12) + tx * wm.R22) * (-(ty * wm.dr1_R12) + tx * wm.dr1_R22)
    //               + 2 * (wm.R12 - tx * wm.R32) * (wm.dr1_R12 - tx * wm.dr1_R32)
    //               + 2 * (-wm.R22 + ty * wm.R32) * (-wm.dr1_R22 + ty * wm.dr1_R32)))
    //         / common_2
    //         + (common_3
    //            * ((-wm.R22 + ty * wm.R32) * (-(sx * wm.dr1_R11) - sy * wm.dr1_R12 - sz * wm.dr1_R13)
    //               + (-az + bz - sx * wm.R31 - sy * wm.R32 - sz * wm.R33) * (-(ty * wm.dr1_R12) + tx * wm.dr1_R22)
    //               + (wm.R12 - tx * wm.R32) * (-(sx * wm.dr1_R21) - sy * wm.dr1_R22 - sz * wm.dr1_R23)
    //               + (-ay + by - sx * wm.R21 - sy * wm.R22 - sz * wm.R23) * (wm.dr1_R12 - tx * wm.dr1_R32)
    //               + (-ax + bx - sx * wm.R11 - sy * wm.R12 - sz * wm.R13) * (-wm.dr1_R22 + ty * wm.dr1_R32)
    //               + (-(ty * wm.R12) + tx * wm.R22) * (-(sx * wm.dr1_R31) - sy * wm.dr1_R32 - sz * wm.dr1_R33)))
    //         / (common_1 * common_0);
    // }
    //
    // constexpr auto dr_dtheta() -> T
    // {
    //     return -0.5
    //         * (common_0
    //            * (2 * (-(ty * wm.R12) + tx * wm.R22) * (-(ty * wm.dr2_R12) + tx * wm.dr2_R22)
    //               + 2 * (wm.R12 - tx * wm.R32) * (wm.dr2_R12 - tx * wm.dr2_R32)
    //               + 2 * (-wm.R22 + ty * wm.R32) * (-wm.dr2_R22 + ty * wm.dr2_R32)))
    //         / common_2
    //         + (common_3
    //            * ((-wm.R22 + ty * wm.R32) * (-(sx * wm.dr2_R11) - sy * wm.dr2_R12 - sz * wm.dr2_R13)
    //               + (-az + bz - sx * wm.R31 - sy * wm.R32 - sz * wm.R33) * (-(ty * wm.dr2_R12) + tx * wm.dr2_R22)
    //               + (wm.R12 - tx * wm.R32) * (-(sx * wm.dr2_R21) - sy * wm.dr2_R22 - sz * wm.dr2_R23)
    //               + (-ay + by - sx * wm.R21 - sy * wm.R22 - sz * wm.R23) * (wm.dr2_R12 - tx * wm.dr2_R32)
    //               + (-ax + bx - sx * wm.R11 - sy * wm.R12 - sz * wm.R13) * (-wm.dr2_R22 + ty * wm.dr2_R32)
    //               + (-(ty * wm.R12) + tx * wm.R22) * (-(sx * wm.dr2_R31) - sy * wm.dr2_R32 - sz * wm.dr2_R33)))
    //         / (common_1 * common_0);
    // }
    //
    // constexpr auto dr_dphi() -> T
    // {
    //     return -0.5
    //         * (common_0
    //            * (2 * (-(ty * wm.R12) + tx * wm.R22) * (-(ty * wm.dr3_R12) + tx * wm.dr3_R22)
    //               + 2 * (wm.R12 - tx * wm.R32) * (wm.dr3_R12 - tx * wm.dr3_R32)
    //               + 2 * (-wm.R22 + ty * wm.R32) * (-wm.dr3_R22 + ty * wm.dr3_R32)))
    //         / common_2
    //         + (common_3
    //            * ((-wm.R22 + ty * wm.R32) * (-(sx * wm.dr3_R11) - sy * wm.dr3_R12 - sz * wm.dr3_R13)
    //               + (-az + bz - sx * wm.R31 - sy * wm.R32 - sz * wm.R33) * (-(ty * wm.dr3_R12) + tx * wm.dr3_R22)
    //               + (wm.R12 - tx * wm.R32) * (-(sx * wm.dr3_R21) - sy * wm.dr3_R22 - sz * wm.dr3_R23)
    //               + (-ay + by - sx * wm.R21 - sy * wm.R22 - sz * wm.R23) * (wm.dr3_R12 - tx * wm.dr3_R32)
    //               + (-ax + bx - sx * wm.R11 - sy * wm.R12 - sz * wm.R13) * (-wm.dr3_R22 + ty * wm.dr3_R32)
    //               + (-(ty * wm.R12) + tx * wm.R22) * (-(sx * wm.dr3_R31) - sy * wm.dr3_R32 - sz * wm.dr3_R33)))
    //         / (common_1 * common_0);
    // }
    //
    // constexpr auto dr_dax() -> T { return ((wm.R22 - ty * wm.R32) * common_3) / (common_1 * common_0); }
    //
    // constexpr auto dr_day() -> T { return ((-wm.R12 + tx * wm.R32) * common_3) / (common_1 * common_0); }
    //
    // constexpr auto dr_daz() -> T { return ((ty * wm.R12 - tx * wm.R22) * common_3) / (common_1 * common_0); }
    //
    // constexpr auto dr_dbx() -> T { return ((-wm.R22 + ty * wm.R32) * common_3) / (common_1 * common_0); }
    //
    // constexpr auto dr_dby() -> T { return ((wm.R12 - tx * wm.R32) * common_3) / (common_1 * common_0); }
    //
    // constexpr auto dr_dbz() -> T { return ((-(ty * wm.R12) + tx * wm.R22) * common_3) / (common_1 * common_0); }
    //
    // constexpr auto dr_dtx() -> T
    // {
    //     return ((-((-ay + by - sx * wm.R21 - sy * wm.R22 - sz * wm.R23) * wm.R32)
    //              + wm.R22 * (-az + bz - sx * wm.R31 - sy * wm.R32 - sz * wm.R33))
    //             * common_3)
    //         / (common_1 * common_0)
    //         - ((2 * wm.R22 * (-(ty * wm.R12) + tx * wm.R22) - 2 * wm.R32 * (wm.R12 - tx * wm.R32)) * common_0) / (2 * common_2);
    // }
    //
    // constexpr auto dr_dty() -> T
    // {
    //     return (((-ax + bx - sx * wm.R11 - sy * wm.R12 - sz * wm.R13) * wm.R32
    //              - wm.R12 * (-az + bz - sx * wm.R31 - sy * wm.R32 - sz * wm.R33))
    //             * common_3)
    //         / (common_1 * common_0)
    //         - ((-2 * wm.R12 * (-(ty * wm.R12) + tx * wm.R22) + 2 * wm.R32 * (-wm.R22 + ty * wm.R32)) * common_0) / (2 * common_2);
    // }
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

    auto is_free() const -> bool { return kind == Kind::FREE; }
    operator T() const { return value; }
};

template<typename T>
auto operator<<(std::ostream& ofs, const SA::parameter<T>& rhs) -> std::ostream&
{
    return ofs << std::setw(2) << (rhs.kind == SA::Kind::FREE ? '~' : 'x') << std::setw(12) << rhs.value;
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

    auto dump_pede_param(int layer, std::ostream& ofs)
    {
        for (int i = 0; i < 6; ++i) {
            const auto& p = params[i];
            if (p.is_free()) {
                ofs << std::setw(10) << layer * 10 + i + 1 << std::setw(20) << p.value << std::setw(10) << 0.0 << '\n';
            }
        }
    }
};

template<template<class> class R>
class MilleBuilder
{
  public:
    MilleBuilder(const char* prefix, const char* outFileName, bool asBinary = true, bool writeZero = false)
        : prefix(prefix)
        , mille(outFileName, asBinary, writeZero)
    {
    }
    ~MilleBuilder() {};

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
                            float gamma)
    {
        layers_global_pars.emplace_back(gx, gy, gz, ga, gb, gc);
        layers_derivatives.emplace_back(gx.value, gy.value, gz.value, ga.value, gb.value, gc.value, ax, ay, az, alpha, beta, gamma);
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
     * @param sx wire position x
     * @param sy wire position y
     * @param sz wire position z
     */
    auto add_local(int layer, float bx, float by, float bz, float tx, float ty, float sx, float sy, float sz, float dr, float sigma)
    {
        assert(layer < layers_global_pars.size());

        const auto param_idx_offset = layer * 10;
        const auto& current_layer = layers_global_pars[layer];

        auto& derivs = layers_derivatives[layer];
        derivs.update(sx, sy, sz, bx, by, bz, tx, ty);

        std::array<float, 5> local_derivatives;

        local_derivatives[0] = derivs.dr_dbx();
        local_derivatives[1] = derivs.dr_dby();
        local_derivatives[2] = derivs.dr_dbz();
        local_derivatives[3] = derivs.dr_dtx();
        local_derivatives[4] = derivs.dr_dty();

        std::array<float, 6> global_derivatives;
        std::array<int, 6> global_deriv_index;
        int global_cnt = 0;

        if (current_layer.gx().is_free()) {
            global_derivatives[global_cnt] = derivs.dr_dgx();
            global_deriv_index[global_cnt] = param_idx_offset + 1;
            global_cnt++;
        }

        if (current_layer.gy().is_free()) {
            global_derivatives[global_cnt] = derivs.dr_dgy();
            global_deriv_index[global_cnt] = param_idx_offset + 2;
            global_cnt++;
        }

        if (current_layer.gz().is_free()) {
            global_derivatives[global_cnt] = derivs.dr_dgz();
            global_deriv_index[global_cnt] = param_idx_offset + 3;
            global_cnt++;
        }

        if (current_layer.ga().is_free()) {
            global_derivatives[global_cnt] = derivs.dr_dga();
            global_deriv_index[global_cnt] = param_idx_offset + 4;
            global_cnt++;
        }

        if (current_layer.gb().is_free()) {
            global_derivatives[global_cnt] = derivs.dr_dgb();
            global_deriv_index[global_cnt] = param_idx_offset + 5;
            global_cnt++;
        }

        if (current_layer.gc().is_free()) {
            global_derivatives[global_cnt] = derivs.dr_dgc();
            global_deriv_index[global_cnt] = param_idx_offset + 6;
            global_cnt++;
        }

        auto alignment_rotation = euler::make_rotation_matrix(derivs.wm);
        auto sloc = XYZVector(sx, sy, sz);
        auto alignment_translation = XYZVector(derivs.ax, derivs.ay, derivs.az);
        auto alignment_translation_correction = XYZVector(current_layer.gx(), current_layer.gy(), current_layer.gz());
        auto alignment_rotation_correction = Rotation3D(1,
                                                        -current_layer.gc(),
                                                        current_layer.gb(),
                                                        current_layer.gc(),
                                                        1,
                                                        current_layer.ga(),
                                                        -current_layer.gb(),
                                                        current_layer.ga(),
                                                        1);
        auto slab = alignment_rotation_correction * alignment_rotation * sloc + alignment_translation + alignment_translation_correction;
        auto srot = alignment_rotation_correction * alignment_rotation * ROOT::Math::XYZVector(0, 1, 0);
        auto res = abs(geom::distance({bx, by, bz}, {tx, ty, 1.0}, {slab.x(), slab.y(), slab.z()}, srot) - dr);

        if (verbose > 1) {
            std::cout << "Rotation: " << alignment_rotation;
            std::cout << "% sloc=" << sloc << "   slab=" << slab << '\n';
        }
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

        mille.mille(5, local_derivatives.data(), global_cnt, global_derivatives.data(), global_deriv_index.data(), res, sigma);
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

    auto write_param_file()
    {
        std::ofstream param_file(prefix + std::string("params.txt"));
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
            double x, y, z, psi, theta, phi;
            int bx, by, bz, bpsi, btheta, bphi;
            while (true) {
                infile >> x >> y >> z >> psi >> theta >> phi;
                infile >> bx >> by >> bz >> bpsi >> btheta >> bphi;
                if (infile.eof())
                    break;

                add_planes_globals({0, to_kind(bx)},
                                   {0, to_kind(by)},
                                   {0, to_kind(bz)},
                                   {0, to_kind(bpsi)},
                                   {0, to_kind(btheta)},
                                   {0, to_kind(bpsi)},
                                   x,
                                   y,
                                   z,
                                   bpsi * TMath::DegToRad(),
                                   btheta * TMath::DegToRad(),
                                   bphi * TMath::DegToRad());
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

                add_planes_globals({0, to_kind(bx)},
                                   {0, to_kind(by)},
                                   {0, to_kind(bz)},
                                   {0, to_kind(bpsi)},
                                   {0, to_kind(btheta)},
                                   {0, to_kind(bpsi)},
                                   x,
                                   y,
                                   z,
                                   ra.Psi(),
                                   ra.Theta(),
                                   ra.Phi());
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
        std::cout << std::left << std::setw(width) << "Layer" << std::setw(width) << "Gx" << std::setw(width) << "Gy" << std::setw(width)
                  << "GzZ" << std::setw(width) << "Ga" << std::setw(width) << "Gb" << std::setw(width) << "Gc" << '\n';
        std::cout << bar << '\n';

        auto i = 1;
        for (const auto& gl : layers_global_pars) {
            std::cout << std::setw(width) << i++ << gl.gx() << gl.gy() << gl.gz() << gl.ga() << gl.gb() << gl.gc() << '\n';
        }

        std::cout << bar << '\n';
        std::cout << std::right;
    }

    auto set_verbose(int make_verbose) { verbose = make_verbose; }

  private:
    std::string prefix;
    std::vector<global_parameters<float>> layers_global_pars;
    std::vector<derivatives<float, R>> layers_derivatives;
    Mille mille;
    int verbose {0};
};

}  // namespace SA
