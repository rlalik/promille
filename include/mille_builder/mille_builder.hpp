#pragma once

// internal
#include <mille_builder/alignment_model.hpp>
#include <mille_builder/euler_angles.hpp>

// millepede
#include "Mille.h"

// ROOT
#include <Math/EulerAngles.h>
#include <Math/Point3D.h>
#include <Math/Rotation3D.h>
#include <Math/RotationZYX.h>
#include <Math/Vector3D.h>
#include <TMath.h>

// system
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

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
struct parameter_state;

template<typename T>
struct global_parameter
{
    int id;
    T value;
    std::string description;

    global_parameter(int id, T value, std::string description = "")
        : id(id)
        , value(value)
        , description(std::move(description))
    {
    }
    operator T() const { return value; }

    auto make_state_from_parameter(Kind kind = Kind::FREE) -> parameter_state<global_parameter<T>>
    {
        return parameter_state<global_parameter<T>>(this, kind);
    }
};

template<typename T>
auto operator<<(std::ostream& ofs, const mb::global_parameter<T>& rhs) -> std::ostream&
{
    return ofs << std::setw(6) << rhs.id << std::setw(12) << rhs.value << "  " << rhs.description;
}

template<typename T>
struct local_parameter
{
    int id;
    std::string description;

    local_parameter(int id, std::string description = "")
        : id(id)
        , description(std::move(description))
    {
    }

    auto make_state_from_parameter(Kind kind = Kind::FREE) -> parameter_state<local_parameter<T>>
    {
        return parameter_state<local_parameter<T>>(this, kind);
    }
};

template<typename T>
auto operator<<(std::ostream& ofs, const mb::local_parameter<T>& rhs) -> std::ostream&
{
    return ofs << std::setw(6) << rhs.id << "  " << rhs.description;
}

template<typename T>
struct parameter_state
{
    T* ptr;
    Kind kind;

    explicit constexpr parameter_state(T* ptr, Kind kind = Kind::FREE)
        : ptr(ptr)
        , kind(kind)
    {
    }

    auto is_free() const -> bool { return kind == Kind::FREE; }

    auto operator=(Kind rhs) -> void { kind = rhs; }
};

template<typename T>
auto operator<<(std::ostream& ofs, const mb::parameter_state<T>& rhs) -> std::ostream&
{
    return ofs << std::setw(6) << rhs.ptr->id << rhs.kind << std::setw(12) << rhs.ptr->value << "  " << rhs.ptr->description;
}

template<typename T, std::size_t Nglobal, std::size_t Nlocal>
struct measurement_plane
{
    std::array<parameter_state<global_parameter<T>>, Nglobal> globals;
    std::array<parameter_state<local_parameter<T>>, Nlocal> locals;

    measurement_plane(std::array<parameter_state<global_parameter<T>>, Nglobal>&& g,
                      std::array<parameter_state<local_parameter<T>>, Nlocal>&& l)
        : globals(std::move(g))
        , locals(std::move(l))
    {
    }

    template<typename... StateKind>
    auto set_globals_configuration(StateKind... kinds) -> measurement_plane<T, Nglobal, Nlocal>&
    {
        // parameter_state<float> ps(nullptr);
        // globals = {kinds...};
        set_params_configuration(globals, std::make_index_sequence<sizeof...(StateKind)> {}, kinds...);
        return *this;
    }

    template<typename... StateKind>
    auto set_locals_configuration(StateKind... kinds) -> measurement_plane<T, Nglobal, Nlocal>&
    {
        set_params_configuration(locals, std::make_index_sequence<sizeof...(StateKind)> {}, kinds...);
        return *this;
    }

  private:
    template<typename Param, typename... StateKind, std::size_t... Is>
    auto set_params_configuration(Param& p, std::index_sequence<Is...> const&, StateKind... kinds) -> measurement_plane<T, Nglobal, Nlocal>&
    {
        ((set_param_kind(p, Is, kinds)), ...);
        return *this;
    }

    template<typename Param>
    auto set_param_kind(Param& p, size_t idx, Kind kind)
    {
        p.at(idx).kind = kind;
    }
};

template<typename ResiduaModel>
class mille_builder
{
  public:
    using main_type = float;

    using global_param_type_t = global_parameter<main_type>;
    using local_param_type_t = local_parameter<main_type>;

    using measurement_plane_t = measurement_plane<main_type, ResiduaModel::n_globals, ResiduaModel::n_locals>;
    using measurement_plane_array_t = std::vector<measurement_plane_t>;

    using local_parameter_array_t = std::array<std::unique_ptr<local_param_type_t>, ResiduaModel::n_locals>;

  public:
    mille_builder(const char* prefix, const char* outFileName, bool asBinary = true, bool writeZero = false)
        : mille_prefix(prefix)
        , mille(outFileName, asBinary, writeZero)
        , local_parameters_list(make_locals(std::make_integer_sequence<decltype(ResiduaModel::n_locals), ResiduaModel::n_locals> {}))
    {
    }
    ~mille_builder() {}

    auto add_global_parameter(int id, float value, std::string description = "") -> void
    {
        global_parameters_list.push_back(std::make_unique<global_param_type_t>(id, value, std::move(description)));
        global_parameters_map[id] = global_parameters_list.back().get();
    }

    auto set_local_parameter(int id, std::string description) -> void
    {
        local_parameters_list.at(id)->description = std::move(description);
    }

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

    template<typename... GlobalIds>
    auto add_plane(int plane_id, GlobalIds... gids) -> measurement_plane_t&
    {
        // static_assert(sizeof...(gids) == ResiduaModel::n_globals, "Wrong number of global parameters ids");

        // printf("SIZES %ld %ld\n", sizeof...(gids), sizeof...(pd));

        // {global_parameters_map.at(gids)->make_state_from_parameter()...},

        plane_globals_map[plane_id] = plane_configuration.size();
        plane_configuration.emplace_back(make_plane_globals_state(gids...),
                                         make_plane_locals_state<decltype(ResiduaModel::n_locals)>(
                                             std::make_integer_sequence<decltype(ResiduaModel::n_locals), ResiduaModel::n_locals> {}));

        // plane_configuration.back().globals = {global_parameters_map.at(gids)->make_state_from_parameter()...};
        // plane_configuration.back().globals = { global_parameters_list[global_parameters_map.at(gids)].get() ...}; FIXME DO I NEED IT?

        // plane_configuration.back().locals = detail::create_array(Kind::FREE, std::make_index_sequence<ResiduaModel::n_locals>()); FIXME
        // !!!
        global_planes.push_back(ResiduaModel(global_parameters_map.at(gids)->value...));

        return plane_configuration.back();

        // TODO
        // plane_id2idx_mapping[globals_id] = layers_global_pars.size();

        // layers_global_pars.emplace_back(gx, gy, gz, ga, gb, gc);
        // layers_derivatives.emplace_back(gx.value, gy.value, gz.value, ga.value, gb.value, gc.value, ax, ay, az, alpha, beta, gamma);
        // alignment_rotation_corrections.emplace_back(1, -gc, gb, gc, 1, ga, -gb, ga, 1);
        // layers_local_pars.emplace_back(lx0, ly0, ltx, lty);
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
    template<typename... LocalPars>
    auto add_measurement(int plane_id, float dr, float sigma, LocalPars... args) -> void
    {
        auto& residua_model = global_planes[plane_globals_map.at(plane_id)];

        auto distance = residua_model.get_distance(args...);

        std::array<float, ResiduaModel::n_globals> global_derivatives;
        std::array<int, ResiduaModel::n_globals> global_deriv_index;
        int global_cnt = 0;

        const auto param_idx_offset = plane_id * encoder;

        for (int i = 0; i < ResiduaModel::n_globals; ++i) {
            const auto& param_state = plane_configuration[plane_globals_map.at(plane_id)].globals[i];

            if (param_state.is_free()) {
                global_derivatives[global_cnt] = residua_model.global_derivative(i + 1);
                global_deriv_index[global_cnt] = param_state.ptr->id;
                global_cnt++;
            }
        }

        residua_model.print();

        auto plane_idx = plane_globals_map.at(plane_id);
        for (const auto& global_par : plane_configuration[plane_idx].globals) {
            const auto& global_parameter = *global_par.ptr;

            printf("Parameter [%d] = %f for plane %d\n", plane_idx, global_parameter.value, plane_id);

            if (global_par.is_free()) {
            }
            //     global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgx();
            //     global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 1;
            //     global_cnt++;
        }

        // const auto it_idx = plane_id2idx_mapping.find(globals_id);
        // assert(it_idx != plane_id2idx_mapping.cend());
        // auto idx = it_idx->second;
        //
        // auto& derivs = layers_derivatives[idx];
        // derivs.update(sx, sy, sz, bx, by, tx, ty);
        //
        // std::array<float, 4> local_derivatives;
        // int local_cnt = 0;
        //
        // const auto& current_plane_locals = layers_local_pars[idx];
        // if (current_plane_locals.lx0() == Kind::FREE) {
        //     local_derivatives[bits::int2size_t(local_cnt)] = derivs.dr_dbx();
        //     local_cnt++;
        // }
        // if (current_plane_locals.ly0() == Kind::FREE) {
        //     local_derivatives[bits::int2size_t(local_cnt)] = derivs.dr_dby();
        //     local_cnt++;
        // }
        // if (current_plane_locals.ltx() == Kind::FREE) {
        //     local_derivatives[bits::int2size_t(local_cnt)] = derivs.dr_dtx();
        //     local_cnt++;
        // }
        // if (current_plane_locals.lty() == Kind::FREE) {
        //     local_derivatives[bits::int2size_t(local_cnt)] = derivs.dr_dty();
        //     local_cnt++;
        // }
        //
        // std::array<float, 6> global_derivatives;
        // std::array<int, 6> global_deriv_index;
        // int global_cnt = 0;
        // const auto param_idx_offset = globals_id * 10;
        //
        // const auto& current_plane_globals = layers_global_pars[idx];
        // if (current_plane_globals.gx().is_free()) {
        //     global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgx();
        //     global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 1;
        //     global_cnt++;
        // }
        //
        // if (current_plane_globals.gy().is_free()) {
        //     global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgy();
        //     global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 2;
        //     global_cnt++;
        // }
        //
        // if (current_plane_globals.gz().is_free()) {
        //     global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgz();
        //     global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 3;
        //     global_cnt++;
        // }
        //
        // if (current_plane_globals.ga().is_free()) {
        //     global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dga();
        //     global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 4;
        //     global_cnt++;
        // }
        //
        // if (current_plane_globals.gb().is_free()) {
        //     global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgb();
        //     global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 5;
        //     global_cnt++;
        // }
        //
        // if (current_plane_globals.gc().is_free()) {
        //     global_derivatives[bits::int2size_t(global_cnt)] = derivs.dr_dgc();
        //     global_deriv_index[bits::int2size_t(global_cnt)] = param_idx_offset + 6;
        //     global_cnt++;
        // }
        //
        // auto alignment_rotation = euler::make_rotation_matrix(derivs.wm);
        // auto sloc = XYZVector(sx, sy, sz);
        // auto alignment_translation = XYZVector(derivs.a_x, derivs.a_y, derivs.a_z);
        // auto alignment_translation_correction =
        //     XYZVector(current_plane_globals.gx(), current_plane_globals.gy(), current_plane_globals.gz());
        // auto alignment_rotation_correction = alignment_rotation_corrections[idx];
        // auto slab = alignment_rotation_correction * alignment_rotation * sloc + alignment_translation + alignment_translation_correction;
        // auto srot = alignment_rotation_correction * alignment_rotation * ROOT::Math::XYZVector(0, 1, 0);
        // auto res = abs(bits::d2f(geom::distance({bx, by, 0}, {tx, ty, 1.0}, {slab.x(), slab.y(), slab.z()}, srot)) - dr);

        if (verbose) {
            std::cout << " --- GLOBALS ID " << plane_id << " ---\n";
        }
        if (verbose > 1) {
            // derivs.print();
            // std::cout << "# slab=" << slab << '\n';
        }
        if (verbose) {
            static const size_t value_width = 12;
            // std::cout << "+ " << std::setw(3) << globals_id << "   res=" << std::setw(value_width) << res
            // << "   sigma=" << std::setw(value_width) << sigma << "   NLC=" << local_derivatives.size();
            // for (const auto& ld : local_derivatives)
            // std::cout << "  " << std::setw(value_width) << ld;

            std::cout << "   NGL=" << global_cnt;
            for (int i = 0; i < global_cnt; ++i)
                std::cout << "  " << std::setw(value_width) << global_derivatives[i];
            std::cout << '\n';
        }
        //
        // mille.mille(
        //     local_cnt, local_derivatives.data(), global_cnt, global_derivatives.data(), global_deriv_index.data(), (float)res, sigma);
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
        // TODO
        // param_file << "Parameter\n";
        // auto max_globals = layers_global_pars.size();
        // for (decltype(max_globals) i = 0; i < max_globals; ++i) {
        //     // if (layers_global_pars[i].is_free())  {}
        //     layers_global_pars[i].dump_pede_param(i, param_file);
        // }
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
            int cnt = 0;
            while (true) {
                infile >> translation_x >> translation_y >> translation_z;
                infile >> rotation_alpha >> rotation_beta >> rotation_gamma;
                infile >> flag_x >> flag_y >> flag_z >> flag_alpha >> flag_beta >> flag_gamma;
                if (infile.eof()) {
                    break;
                }

                add_plane(cnt++,
                          {1, 0, to_kind(flag_x)},  // FIXME param ids are wrong
                          {2, 0, to_kind(flag_y)},
                          {3, 0, to_kind(flag_z)},
                          {4, 0, to_kind(flag_alpha)},
                          {5, 0, to_kind(flag_beta)},
                          {6, 0, to_kind(flag_gamma)},
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

    auto print(bool /*use_deg*/) const -> void
    {
        constexpr auto width = 14;
        const auto bar = std::string(width * 7, '=');
        std::cout << bar << '\n';

        std::cout << "Global parameters" << '\n';
        std::cout << bar << '\n';
        for (const auto& par : global_parameters_list) {
            std::cout << *par << '\n';
        }
        std::cout << bar << '\n';

        std::cout << "Local parameters" << '\n';
        std::cout << bar << '\n';
        for (const auto& par : local_parameters_list) {
            std::cout << *par << '\n';
        }
        std::cout << bar << '\n';

        std::cout << "Measurement planes" << '\n';
        std::cout << bar << '\n';
        std::cout << std::right;
        for (const auto& plane_ids_it : plane_globals_map) {
            std::cout << '[' << std::setw(5) << plane_ids_it.first << ']';

            for (const auto id : plane_configuration[plane_ids_it.second].globals) {
                std::cout << std::setw(5) << id.ptr->id << id.kind;
            }

            std::cout << "  |";

            for (const auto id : plane_configuration[plane_ids_it.second].locals) {
                std::cout << std::setw(5) << id.ptr->id << id.kind;
            }

            // for (const auto id : plane_configuration[plane_ids_it.second].locals) {
            //     std::cout << std::setw(5) << local_parameters_list[local_parameters_map.at(id)].name << " (" << std::setw(14) << 0 <<
            //     ')';
            // }

            std::cout << '\n';
        }
        std::cout << std::left;
        std::cout << bar << '\n';

        // constexpr auto width = 14;
        // auto bar = std::string(width * 7, '=');
        // std::cout << bar << '\n';
        // std::cout << std::left << std::setw(width) << "Layer" << std::setw(width) << "Gx" << std::setw(width) << "Gy" << std::setw(width)
        //           << "GzZ" << std::setw(width) << "Ga" << std::setw(width) << "Gb" << std::setw(width) << "Gc" << '\n';
        // std::cout << bar << '\n';
        //
        // auto cnt = 1;
        // for (const auto& global_par : layers_global_pars) {
        //     std::cout << std::setw(width) << cnt++ << global_par.gx() << global_par.gy() << global_par.gz() << global_par.ga()
        //               << global_par.gb() << global_par.gc() << '\n';
        // }
        //
        // std::cout << bar << '\n';
        //
        // cnt = 1;
        // for (const auto& local_par : layers_local_pars) {
        //     std::cout << std::setw(width) << cnt++ << local_par.lx0() << local_par.ly0() << local_par.ltx() << local_par.lty() << '\n';
        //     ;
        // }
        //
        // std::cout << bar << '\n';
        // std::cout << std::right;
    }

    template<typename... LocalsKind>
    auto set_locals_configuration(LocalsKind... kind) -> void
    {
    }

    auto set_verbose(int make_verbose) -> void { verbose = make_verbose; }

  private:
    // template<typename Array, std::size_t... I>
    // auto a2t_impl(const Array& a, std::index_sequence<I...>)
    // {
    //     return std::make_tuple(a[I]...);
    // }

    template<typename... Is>
    auto make_plane_globals_state(Is... ints) -> decltype(measurement_plane_t::globals)
    {
        return {global_parameters_map.at(ints)->make_state_from_parameter()...};
    }

    template<typename Is, Is... ints>
    auto make_plane_locals_state(std::integer_sequence<Is, ints...> int_seq) -> decltype(measurement_plane_t::locals)
    {
        return {local_parameters_list.at(ints)->make_state_from_parameter()...};
    }

    template<typename Is, Is... ints>
    auto make_locals(std::integer_sequence<Is, ints...>) -> local_parameter_array_t
    {
        return {std::make_unique<local_parameter<main_type>>(ints)...};
    }

  private:
    std::string mille_prefix;

    // std::map<int, size_t> plane_id2idx_mapping;
    // std::vector<global_parameters<float>> layers_global_pars;
    // // std::vector<derivatives<float, R>> layers_derivatives;
    // std::vector<Rotation3D> alignment_rotation_corrections;
    // std::vector<local_parameters> layers_local_pars;

    std::vector<std::unique_ptr<global_param_type_t>> global_parameters_list;
    std::map<int, global_param_type_t*> global_parameters_map;

    local_parameter_array_t local_parameters_list;

    std::vector<ResiduaModel> global_planes;
    std::map<int, int> plane_globals_map;

    measurement_plane_array_t plane_configuration;

    Mille mille;
    int encoder {100};
    int verbose {0};
};

}  // namespace mb
