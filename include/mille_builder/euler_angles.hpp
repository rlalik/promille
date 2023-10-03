#pragma once

#include <Math/Rotation3D.h>

namespace mb::euler
{

template<typename EulerStruct>
auto make_rotation_matrix(EulerStruct es) -> ROOT::Math::Rotation3D
{
    return ROOT::Math::Rotation3D(es.R11, es.R12, es.R13, es.R21, es.R22, es.R23, es.R31, es.R32, es.R33);
}

template<typename T>
struct euler_base
{
    T c1, s1, c2, s2, c3, s3;
    T R11 {0}, R12 {0}, R13 {0}, R21 {0}, R22 {0}, R23 {0}, R31 {0}, R32 {0}, R33 {0};

    euler_base(T a1, T a2, T a3)
        : c1(cos(a1))
        , s1(sin(a1))
        , c2(cos(a2))
        , s2(sin(a2))
        , c3(cos(a3))
        , s3(sin(a3))
    {
    }
};

template<typename T>
struct zyz : euler_base<T>
{
    using euler_base<T>::c1;
    using euler_base<T>::s1;
    using euler_base<T>::c2;
    using euler_base<T>::s2;
    using euler_base<T>::c3;
    using euler_base<T>::s3;

    zyz(T a1, T a2, T a3)
        : euler_base<T>(a1, a2, a3)
    {
        this->R11 = c1 * c2 * c3 - s1 * s3;
        this->R12 = -(c3 * s1) - c1 * c2 * s3;
        this->R13 = c1 * s2;
        this->R21 = c2 * c3 * s1 + c1 * s3;
        this->R22 = c1 * c3 - c2 * s1 * s3;
        this->R23 = s1 * s2;
        this->R31 = -(c3 * s2);
        this->R32 = s2 * s3;
        this->R33 = c2;
    }
};

};  // namespace mb::euler
