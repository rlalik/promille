#pragma once

#include <Math/Rotation3D.h>

namespace promille::euler
{

template<typename EulerStruct>
auto make_rotation_matrix(EulerStruct es) -> ROOT::Math::Rotation3D
{
    return ROOT::Math::Rotation3D(es.R11, es.R12, es.R13, es.R21, es.R22, es.R23, es.R31, es.R32, es.R33);
}

template<typename T>
struct euler_base
{
    T a1, a2, a3;  // angles 1,2,3 of rotation
    T c1, s1, c2, s2, c3, s3;  // cosinus and sinus of angles 1,2,3
    T R11 {0}, R12 {0}, R13 {0}, R21 {0}, R22 {0}, R23 {0}, R31 {0}, R32 {0}, R33 {0};  // rotation matrix elements

    auto determinat() const -> T
    {
        return R11 * R22 * R33 + R12 * R23 * R31 + R13 * R21 * R32 - R13 * R22 * R31 - R12 * R21 * R33 - R11 * R23 * R32;
    }

    auto print() const -> void
    {
        printf("E-ANGLES: %f %f %f\n", a1, a2, a3);
        printf("E-TRIG: %f %f   %f  %f   %f %f\n", s1, c1, s2, c2, s3, c3);
        printf("E-MATRIX: %f %f %f %f %f %f %f %f %f\n", R11, R12, R13, R21, R22, R23, R31, R32, R33);
    }

  protected:
    euler_base(T a1, T a2, T a3)
        : a1(a1)
        , a2(a2)
        , a3(a3)
    {
        init_sin_cos();
    }

    euler_base(T r11, T r12, T r13, T r21, T r22, T r23, T r31, T r32, T r33)
        : R11(r11)
        , R12(r12)
        , R13(r13)
        , R21(r21)
        , R22(r22)
        , R23(r23)
        , R31(r31)
        , R32(r32)
        , R33(r33)
    {
    }

    auto init_sin_cos() -> void
    {
        c1 = cos(a1);
        s1 = sin(a1);
        c2 = cos(a2);
        s2 = sin(a2);
        c3 = cos(a3);
        s3 = sin(a3);
    }
};

template<typename T>
struct zyz : euler_base<T>
{
    using euler_base<T>::a1;
    using euler_base<T>::a2;
    using euler_base<T>::a3;

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

    zyz(T r11, T r12, T r13, T r21, T r22, T r23, T r31, T r32, T r33)
        : euler_base<T>(r11, r12, r13, r21, r22, r23, r31, r32, r33)
    {
        if (abs(r33) != 1.0) {
            this->a1 = atan2(r23, r13);
            this->a2 = atan2(sqrt(1.0 - r33 * r33), r33);
            this->a3 = atan2(r32, -r31);
        } else {
            this->a1 = atan2(r21, r11);
            this->a2 = atan2(sqrt(1.0 - r33 * r33), r33);
            this->a3 = 0;
        }

        this->init_sin_cos();
    }
};

};  // namespace promille::euler
