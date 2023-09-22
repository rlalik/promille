#pragma once

#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include <Math/Point3Dfwd.h>
#include <Math/Rotation3D.h>
#include <Math/Vector3D.h>

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
template<typename T> constexpr auto R11(T phi, T theta, T ypsilon) -> T { return cos(theta)*cos(phi); }
template<typename T> constexpr auto R12(T phi, T theta, T ypsilon) -> T { return cos(theta)*sin(phi); }
template<typename T> constexpr auto R13(T phi, T theta, T ypsilon) -> T { return -sin(theta); }

template<typename T> constexpr auto R21(T phi, T theta, T ypsilon) -> T { return sin(ypsilon)*sin(theta)*cos(phi) - cos(ypsilon)*sin(phi); }
template<typename T> constexpr auto R22(T phi, T theta, T ypsilon) -> T { return sin(ypsilon)*sin(theta)*sin(phi) + cos(ypsilon)*cos(phi); }
template<typename T> constexpr auto R23(T phi, T theta, T ypsilon) -> T { return sin(ypsilon)*cos(theta); }

template<typename T> constexpr auto R31(T phi, T theta, T ypsilon) -> T { return cos(ypsilon)*sin(theta)*cos(phi) + sin(ypsilon)*sin(phi); }
template<typename T> constexpr auto R32(T phi, T theta, T ypsilon) -> T { return cos(ypsilon)*sin(theta)*sin(phi) - sin(ypsilon)*cos(phi); }
template<typename T> constexpr auto R33(T phi, T theta, T ypsilon) -> T { return cos(ypsilon)*cos(theta); }
// clang-format on

}  // namespace helper

template<typename T>
constexpr auto make_rotation_matrix(T phi, T theta, T ypsilon) -> Rotation3D
{
    return Rotation3D(helper::R11(phi, theta, ypsilon),
                      helper::R12(phi, theta, ypsilon),
                      helper::R13(phi, theta, ypsilon),
                      helper::R21(phi, theta, ypsilon),
                      helper::R22(phi, theta, ypsilon),
                      helper::R23(phi, theta, ypsilon),
                      helper::R31(phi, theta, ypsilon),
                      helper::R32(phi, theta, ypsilon),
                      helper::R33(phi, theta, ypsilon));
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

inline auto distance(XYZPoint base1,
                     XYZVector dir1,
                     XYZPoint base2,
                     XYZVector dir2) -> double
{
    return std::abs((base1 - base2).Dot(dir1.Cross(dir2))
                    / dir1.Cross(dir2).R());
}

template<typename T>
struct global_parameters
{
    std::optional<T> X_coarse;
    std::optional<T> Y_coarse;
    std::optional<T> Z_coarse;

    std::optional<T> X_fine;
    std::optional<T> Y_fine;
    std::optional<T> Z_fine;

    std::optional<T> phi;
    std::optional<T> theta;
    std::optional<T> ypsilon;
};

struct MilleBuilder
{
    template<typename T>
    auto add_planes_global(std::optional<T> X_coarse,
                           std::optional<T> Y_coarse,
                           std::optional<T> Z_coarse,
                           std::optional<T> X_fine,
                           std::optional<T> Y_fine,
                           std::optional<T> Z_fine,
                           std::optional<T> phi,
                           std::optional<T> theta,
                           std::optional<T> ypsilon)
    {
        globals.emplace_back(X_coarse,
                             Y_coarse,
                             Z_coarse,
                             X_fine,
                             Y_fine,
                             Z_fine,
                             phi,
                             theta,
                             ypsilon);
    }

    std::vector<global_parameters<double>> globals;
};

}  // namespace SA
