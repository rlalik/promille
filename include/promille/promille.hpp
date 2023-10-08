#pragma once

// millepede
#include "Mille.h"

// system
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace promille
{

template<typename T, std::size_t Nglobal, std::size_t Nlocal, typename PointType, typename VectorType, typename... ExtraArgs>
struct residual_model_base
{
    static const size_t n_globals = Nglobal;
    static const size_t n_locals = Nlocal;

    auto get_residual(PointType base, VectorType track) -> T
    {
        calc_derivatives(base, track);
        return calc_residual(base, track);
    }

    auto get_residual(PointType base, VectorType track, ExtraArgs... args) -> T
    {
        update_extras(args...);
        calc_derivatives(base, track);
        return calc_residual(base, track);
    }

    auto global_derivative(size_t variable_number) const -> T
    {
        assert(variable_number > 0 and variable_number <= Nglobal);
        return global_derivatives.at(variable_number - 1);
    }

    auto local_derivative(size_t variable_number) const -> T
    {
        assert(variable_number > 0 and variable_number <= Nlocal);
        return local_derivatives.at(variable_number - 1);
    }

    auto print() const -> void
    {
        for (size_t i = 0; i < Nglobal; ++i) {
            std::cout << " Global derivative " << i + 1 << " = " << global_derivative(i + 1) << '\n';
        }

        for (size_t i = 0; i < Nlocal; ++i) {
            std::cout << " Local derivative " << i + 1 << " = " << local_derivative(i + 1) << '\n';
        }
    }

  protected:
    template<size_t variable_number>
    auto set_global_derivative() -> T&
    {
        static_assert(variable_number > 0 and variable_number <= Nglobal, "Global parameter index out of bounds");
        return global_derivatives.at(variable_number - 1);
    }

    template<size_t variable_number>
    auto set_local_derivative() -> T&
    {
        static_assert(variable_number > 0 and variable_number <= Nlocal, "Local parameter index out of bounds");
        return local_derivatives.at(variable_number - 1);
    }

  protected:
    constexpr residual_model_base() = default;

  private:
    std::array<T, Nglobal> global_derivatives {0};
    std::array<T, Nlocal> local_derivatives {0};

    virtual auto update_extras(ExtraArgs... args) -> void = 0;
    virtual auto calc_residual(PointType base, VectorType track) -> T = 0;
    virtual auto calc_derivatives(PointType base, VectorType track) -> void = 0;
};

enum class Kind
{
    FIXED,
    FREE
};

inline auto operator<<(std::ostream& ofs, Kind rhs) -> std::ostream&
{
    return ofs << std::setw(2) << (rhs == promille::Kind::FREE ? '~' : 'x');
}

template<typename T>
struct parameter_state;

template<typename T>
struct global_parameter
{
    size_t id;
    T value;
    std::string description;

    global_parameter(size_t id, T value, std::string description = "")
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
auto operator<<(std::ostream& ofs, const promille::global_parameter<T>& rhs) -> std::ostream&
{
    return ofs << std::setw(6) << rhs.id << std::setw(12) << rhs.value << "  " << rhs.description;
}

template<typename T>
struct local_parameter
{
    size_t id;
    std::string description;

    local_parameter(size_t id, std::string description = "")
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
auto operator<<(std::ostream& ofs, const promille::local_parameter<T>& rhs) -> std::ostream&
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
auto operator<<(std::ostream& ofs, const promille::parameter_state<T>& rhs) -> std::ostream&
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
        set_params_configuration(globals, std::make_index_sequence<sizeof...(StateKind)> {}, kinds...);
        return *this;
    }

    auto set_globals_configuration(const std::array<Kind, Nglobal>& kinds) -> measurement_plane<T, Nglobal, Nlocal>&
    {
        for (size_t i = 0; i < kinds.size(); ++i)
            globals[i].kind = kinds[i];
        return *this;
    }

    template<typename... StateKind>
    auto set_locals_configuration(StateKind... kinds) -> measurement_plane<T, Nglobal, Nlocal>&
    {
        set_params_configuration(locals, std::make_index_sequence<sizeof...(StateKind)> {}, kinds...);
        return *this;
    }

    auto set_locals_configuration(const std::array<Kind, Nlocal>& kinds) -> measurement_plane<T, Nglobal, Nlocal>&
    {
        for (size_t i = 0; i < kinds.size(); ++i)
            locals[i].kind = kinds[i];
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
    auto set_param_kind(Param& p, size_t idx, Kind kind) -> void
    {
        p.at(idx).kind = kind;
    }
};

template<typename ResiduaModel>
class promille
{
  public:
    using main_type = float;

    using global_param_type_t = global_parameter<main_type>;
    using local_param_type_t = local_parameter<main_type>;

    using measurement_plane_t = measurement_plane<main_type, ResiduaModel::n_globals, ResiduaModel::n_locals>;
    using measurement_plane_array_t = std::vector<measurement_plane_t>;

    using local_parameter_array_t = std::array<std::unique_ptr<local_param_type_t>, ResiduaModel::n_locals>;

  public:
    /** The Mille Wrapper class.
     *
     * @param prefix
     * @param outFileName the filename to store MillePede data
     * @param asBianry the output file format, see MillePede documentation
     * @param writeZero see MillePede documentation
     */
    promille(const char* prefix, const char* outFileName, bool asBinary = true, bool writeZero = false)
        : mille_prefix(prefix)
        , mille(outFileName, asBinary, writeZero)
        , local_parameters_list(make_locals(std::make_index_sequence<ResiduaModel::n_locals> {}))
    {
    }

    /** Add global apramater with given id and initial value, optionally also description of the global parameter.
     * @param gid global parameter id, must be positive integer larger than 0
     * @param value inttial parameter value
     * @param description the description
     * @return the added id
     */
    auto add_global_parameter(size_t gid, float value, std::string description = "") -> size_t
    {
        global_parameters_list.push_back(std::make_unique<global_param_type_t>(gid, value, std::move(description)));
        global_parameters_map[gid] = global_parameters_list.back().get();
        return gid;
    }

    /** Set local parameter description.
     * @param lid local parameter id, should be in range 1 to N where N is number of the local params in the residual model
     * @param description the param description
     */
    auto set_local_parameter(size_t lid, std::string description) -> void
    {
        local_parameters_list.at(lid)->description = std::move(description);
    }

    /**
     * Add measurement plane global parameters configuration.
     *
     * @param plane_id the measurement plane id
     * @param gids list of global parameters ids associated to this measurement plane
     * @return measurement plane object
     */

    template<typename... GlobalIds>
    auto add_plane(size_t plane_id, GlobalIds... gids) -> measurement_plane_t&
    {
        // static_assert(sizeof...(gids) == ResiduaModel::n_globals, "Wrong number of global parameters ids");

        plane_globals_map[plane_id] = plane_configuration.size();
        plane_configuration.emplace_back(make_plane_globals_state(gids...),
                                         make_plane_locals_state(std::make_index_sequence<ResiduaModel::n_locals> {}));

        global_planes.push_back(ResiduaModel(global_parameters_map.at(gids)->value...));

        return plane_configuration.back();
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
    template<typename... LocalExtraPars>
    auto add_measurement(size_t plane_id, float sigma, LocalExtraPars... args) -> void
    {
        auto& residua_model = global_planes[plane_globals_map.at(plane_id)];

        auto residuum = residua_model.get_residual(args...);

        std::array<float, ResiduaModel::n_globals> global_derivatives;
        std::array<int, ResiduaModel::n_globals> global_deriv_index;
        int global_cnt = 0;

        for (size_t i = 0; i < ResiduaModel::n_globals; ++i) {
            const auto& param_state = plane_configuration[plane_globals_map.at(plane_id)].globals[i];

            if (param_state.is_free()) {
                global_derivatives[global_cnt] = residua_model.global_derivative(i + 1);
                global_deriv_index[global_cnt] = param_state.ptr->id;
                global_cnt++;
            }
        }

        std::array<float, ResiduaModel::n_locals> local_derivatives;
        int local_cnt = 0;

        for (size_t i = 0; i < ResiduaModel::n_locals; ++i) {
            const auto& param_state = plane_configuration[plane_globals_map.at(plane_id)].locals[i];

            if (param_state.is_free()) {
                local_derivatives[local_cnt] = residua_model.local_derivative(i + 1);
                local_cnt++;
            }
        }

        if (verbose) {
            std::cout << " --- GLOBALS ID " << plane_id << " ---\n";
        }
        if (verbose > 1) {
            residua_model.print();
        }
        if (verbose) {
            static const size_t value_width = 12;
            std::cout << "   res=" << std::setw(value_width) << residuum << "   sigma=" << std::setw(value_width) << sigma
                      << "   NLC=" << local_derivatives.size();
            for (const auto& ld : local_derivatives)
                std::cout << "  " << std::setw(value_width) << ld;

            std::cout << "   NGL=" << global_cnt;
            for (size_t i = 0; i < global_cnt; ++i)
                std::cout << std::setw(5) << global_deriv_index[i] << ":" << std::setw(value_width) << global_derivatives[i];
            std::cout << '\n';
        }
        //
        mille.mille(
            local_cnt, local_derivatives.data(), global_cnt, global_derivatives.data(), global_deriv_index.data(), (float)residuum, sigma);
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

    auto set_verbose(int make_verbose) -> void { verbose = make_verbose; }

  private:
    template<typename... Is>
    auto make_plane_globals_state(Is... ints) -> decltype(measurement_plane_t::globals)
    {
        return {global_parameters_map.at(ints)->make_state_from_parameter()...};
    }

    template<std::size_t... Is>
    auto make_plane_locals_state(std::index_sequence<Is...> int_seq) -> decltype(measurement_plane_t::locals)
    {
        return {local_parameters_list.at(Is)->make_state_from_parameter()...};
    }

    template<typename Is, Is... ints>
    auto make_locals(std::integer_sequence<Is, ints...> /* unused */) -> local_parameter_array_t
    {
        return {std::make_unique<local_parameter<main_type>>(ints)...};
    }

    std::string mille_prefix;

    // std::map<int, size_t> plane_id2idx_mapping;
    // std::vector<global_parameters<float>> layers_global_pars;
    // // std::vector<derivatives<float, R>> layers_derivatives;
    // std::vector<Rotation3D> alignment_rotation_corrections;
    // std::vector<local_parameters> layers_local_pars;

    std::vector<std::unique_ptr<global_param_type_t>> global_parameters_list;
    std::map<size_t, global_param_type_t*> global_parameters_map;

    local_parameter_array_t local_parameters_list;

    std::vector<ResiduaModel> global_planes;
    std::map<size_t, size_t> plane_globals_map;

    measurement_plane_array_t plane_configuration;

    Mille mille;
    int verbose {0};
};

}  // namespace promille
