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

template<typename T, size_t Nlocals, size_t Nglobals>
struct residual_model_base
{
    static const size_t n_globals = Nglobals;
    static const size_t n_locals = Nlocals;

    virtual auto residual() const -> T = 0;

    constexpr auto global_derivative(size_t variable_number) const -> T
    {
        assert(variable_number < Nglobals);
        return global_derivatives.at(variable_number);
    }

    constexpr auto local_derivative(size_t variable_number) const -> T
    {
        assert(variable_number < Nlocals);
        return local_derivatives.at(variable_number);
    }

    virtual auto print() const -> void
    {
        std::cout << "MODEL: Ng = " << n_globals << "  Nl = " << n_locals << '\n';
        std::cout << "drd_globals: ";
        for (size_t i = 0; i < n_globals; ++i) {
            std::cout << std::setw(12) << global_derivative(i);
        }
        std::cout << "\ndrd_locals : ";
        for (size_t i = 0; i < n_locals; ++i) {
            std::cout << std::setw(12) << local_derivative(i);
        }
        std::cout << '\n';
    }

    virtual ~residual_model_base() = default;

  protected:
    template<size_t variable_number>
    auto set_global_derivative() -> T&
    {
        static_assert(variable_number < n_globals, "Global parameter index out of bounds");
        return global_derivatives.at(variable_number);
    }

    template<size_t variable_number>
    auto set_local_derivative() -> T&
    {
        static_assert(variable_number < n_locals, "Local parameter index out of bounds");
        return local_derivatives.at(variable_number);
    }

    residual_model_base() = default;

  private:
    std::array<T, Nglobals> global_derivatives;
    std::array<T, Nlocals> local_derivatives;
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
struct parameter_kind;

template<typename T>
struct global_parameter
{
    using type = T;

    size_t id;
    T value;
    std::string description;
    bool is_anywere_free {false};

    global_parameter(size_t id, T value, std::string description = "")
        : id(id)
        , value(value)
        , description(std::move(description))
    {
    }

    global_parameter(const global_parameter&) = default;
    global_parameter(global_parameter&&) = default;

    auto operator=(const global_parameter&) = delete;
    auto operator=(global_parameter&&) = delete;

    operator T() const { return value; }

    auto make_parameter_kind(Kind kind = Kind::FIXED) -> parameter_kind<T> { return parameter_kind<T>(this, kind); }

    auto dump_pede_param(std::ostream& ofs) -> void
    {
        ofs << std::setw(10) << id << std::setw(20) << value << std::setw(10) << 0.0 << '\n';
    }
};

template<typename T>
auto operator<<(std::ostream& ofs, const promille::global_parameter<T>& rhs) -> std::ostream&
{
    return ofs << std::setw(6) << rhs.id << std::setw(12) << rhs.value << "  " << rhs.description;
}

template<typename T>
struct parameter_kind
{
  private:
    global_parameter<T>* ptr;
    Kind kind;

  public:
    explicit constexpr parameter_kind(global_parameter<T>* param_ptr, Kind param_kind = Kind::FREE)
        : ptr(param_ptr)
        , kind(param_kind)
    {
    }

    parameter_kind(const parameter_kind&) = default;
    parameter_kind(parameter_kind&&) = default;

    auto operator=(const parameter_kind&) = delete;
    auto operator=(parameter_kind&&) = delete;

    auto get() -> global_parameter<T>& { return *ptr; }

    auto get() const -> const global_parameter<T>& { return *ptr; }

    auto get_kind() -> Kind { return kind; }

    auto get_kind() const -> Kind { return kind; }

    auto set_kind(Kind value) -> void
    {
        kind = value;
        if (kind == Kind::FREE) {
            ptr->is_anywere_free = true;
        }
    }

    auto is_free() const -> bool { return kind == Kind::FREE; }

    auto operator=(Kind rhs) -> void { kind = rhs; }

    template<typename _T>
    friend auto operator<<(std::ostream& ofs, const promille::parameter_kind<_T>& rhs) -> std::ostream&;
};

template<typename _T>
auto operator<<(std::ostream& ofs, const promille::parameter_kind<_T>& rhs) -> std::ostream&
{
    return ofs << std::setw(6) << rhs.ptr->id << rhs.kind << std::setw(12) << rhs.ptr->value << "  " << rhs.ptr->description;
}

template<typename T, typename ResidualModel>
struct measurement_plane
{
    ResidualModel residual_model;

    std::array<parameter_kind<T>, ResidualModel::n_globals> globals_kind;
    std::array<Kind, ResidualModel::n_locals> locals_kind;

    template<typename... GlobalParameters>
    measurement_plane(Mille* mille_ptr, GlobalParameters&... param)
        : mille(mille_ptr)
        , residual_model(param.value...)
        , globals_kind({param.make_parameter_kind()...})
    {
    }

    auto model() -> ResidualModel& { return residual_model; }

    template<typename... StateKind>
    auto set_globals_configuration(StateKind... kinds) -> measurement_plane<T, ResidualModel>&
    {
        static_assert(sizeof...(kinds) == ResidualModel::n_globals, "Number of globals kinds do not match model definition");

        set_params_configuration(globals_kind, std::make_index_sequence<sizeof...(StateKind)> {}, kinds...);

        return *this;
    }

    auto set_globals_configuration(const std::array<Kind, ResidualModel::n_globals>& kinds) -> measurement_plane<T, ResidualModel>&
    {
        for (size_t i = 0; i < kinds.size(); ++i) {
            globals_kind[i].set_kind(kinds[i]);
        }
        return *this;
    }

    template<typename... StateKind>
    auto set_locals_configuration(StateKind... kinds) -> measurement_plane<T, ResidualModel>&
    {
        locals_kind = {kinds...};
        return *this;
    }

    auto set_locals_configuration(const std::array<Kind, ResidualModel::n_locals>& kinds) -> measurement_plane<T, ResidualModel>&
    {
        locals_kind = kinds;
        return *this;
    }

    /**
     * Add set of local variables for given layer of straws.
     *
     */
    template<typename... ModelExtraArguments>
    auto add_measurement(float sigma, ModelExtraArguments... extra_args) -> measurement_plane<T, ResidualModel>&
    {
        if constexpr (sizeof...(extra_args)) {
            residual_model.recalculate(extra_args...);
        }

        auto residuum = residual_model.residual();

        std::vector<float> global_derivatives(residual_model.n_globals);
        std::vector<int> global_deriv_index(residual_model.n_globals);

        int global_cnt = 0;

        for (int i = 0; i < residual_model.n_globals; ++i) {
            const auto& param_state = globals_kind[i];

            if (param_state.is_free()) {
                global_derivatives[global_cnt] = residual_model.global_derivative(i);
                global_deriv_index[global_cnt] = param_state.get().id;
                global_cnt++;
            }
        }

        std::vector<float> local_derivatives(residual_model.n_locals);
        int local_cnt = 0;

        for (int i = 0; i < residual_model.n_locals; ++i) {
            const auto& param_state = locals_kind[i];

            if (param_state == Kind::FREE) {
                local_derivatives[local_cnt] = residual_model.local_derivative(i);
                local_cnt++;
            }
        }

        if (verbose_flag) {
            residual_model.print();
            static const size_t value_width = 12;
            std::cout << "   res=" << std::setw(value_width) << residuum << "   sigma=" << std::setw(value_width) << sigma
                      << "   NLC=" << local_cnt << std::right;
            for (size_t i = 0; i < local_cnt; ++i)
                std::cout << std::setw(5) << i << ":" << std::setw(value_width) << local_derivatives[i];

            std::cout << std::left << "   NGL=" << global_cnt << std::right;
            for (size_t i = 0; i < global_cnt; ++i)
                std::cout << std::setw(5) << global_deriv_index[i] << ":" << std::setw(value_width) << global_derivatives[i];
            std::cout << std::left << '\n';
        }

        mille->mille(local_cnt,
                     local_derivatives.data(),
                     global_cnt,
                     global_derivatives.data(),
                     global_deriv_index.data(),
                     static_cast<float>(residuum),
                     sigma);

        return *this;
    }

    auto print() const -> void
    {
        constexpr auto width = 14;
        const auto bar = std::string(width * 7, '=');

        residual_model.print();

        std::cout << bar << '\n' << " GLOBAL:";
        for (const auto id : globals_kind) {
            std::cout << std::setw(5) << id.get().id << id.get_kind();
        }

        std::cout << "  | LOCAL:";
        for (const auto kind : locals_kind) {
            std::cout << std::setw(5) << kind;
        }
        std::cout << '\n' << bar << '\n';
    }

    auto set_verbose(int make_verbose) -> measurement_plane<T, ResidualModel>&
    {
        verbose_flag = make_verbose;
        return *this;
    }

  private:
    template<typename Param, typename... StateKind, std::size_t... Is>
    auto set_params_configuration(Param& p, std::index_sequence<Is...> const&, StateKind... kinds) -> measurement_plane<T, ResidualModel>&
    {
        ((p[Is].set_kind(kinds)), ...);
        return *this;
    }

    bool verbose_flag {false};
    Mille* mille {nullptr};
};

template<typename T>
class promille;

template<typename T, typename ResidualModel>
struct model_planes
{
    using residual_model_t = ResidualModel;
    using measurement_plane_t = measurement_plane<T, ResidualModel>;
    using measurement_plane_array_t = std::vector<measurement_plane_t>;

    model_planes(promille<T>* promille_ptr)
        : pro_mille(promille_ptr)
    {
    }

    /**
     * Add measurement plane global parameters configuration.
     *
     * @param plane_id the measurement plane id
     * @param gids list of global parameters ids associated to this measurement plane
     * @return measurement plane object
     */
    template<typename... GlobalIds>
    auto add_plane(size_t plane_id, GlobalIds... globals_ids) -> measurement_plane_t&
    {
        static_assert(ResidualModel::n_globals == sizeof...(globals_ids), "Wrong number of local parameters in model");

        auto new_plane =
            plane_indexes_map.emplace(plane_id, measurement_plane_t(&pro_mille->get_mille(), pro_mille->global_parameter(globals_ids)...));

        if (!new_plane.second)
            throw;

        return new_plane.first->second;
    }

    auto plane(size_t plane_id) -> measurement_plane_t&
    {
        auto the_plane = plane_indexes_map.find(plane_id);
        if (the_plane != plane_indexes_map.end())
            return the_plane->second;

        throw std::runtime_error("No model plane id=" + std::to_string(plane_id));
    }

    auto print() const -> void
    {
        constexpr auto width = 14;
        const auto bar = std::string(width * 7, '=');

        std::cout << "Measurement planes (" << plane_indexes_map.size() << ")\n";
        std::cout << bar << '\n';
        std::cout << std::right;
        for (const auto& plane_ids_it : plane_indexes_map) {
            std::cout << '[' << std::setw(5) << plane_ids_it.first << ']';
            std::cout << " GLOBAL:";

            for (const auto [key, value] : plane_indexes_map) {
                // std::cout << std::setw(5) << id.get().id << id.get_kind();
                value.print();
            }
        }
        std::cout << bar << '\n';
    }

    auto set_verbose(int make_verbose) -> model_planes<T, ResidualModel>&
    {
        verbose_flag = make_verbose;
        return *this;
    }

  private:
    std::map<size_t, measurement_plane_t> plane_indexes_map;

    promille<T>* pro_mille {nullptr};
    bool verbose_flag {false};
};

template<typename T = float>
class promille
{
  public:
    using global_param_type_t = ::promille::global_parameter<T>;

  public:
    /** The Mille Wrapper class.
     *
     * @param prefix
     * @param outFileName the filename to store MillePede data
     * @param asBianry the output file format, see MillePede documentation
     * @param writeZero see MillePede documentation
     */
    promille(const char* prefix, const char* outFileName, bool asBianry = true, bool writeZero = true)
        : mille_prefix(prefix)
        , mille(Mille(outFileName, asBianry, writeZero))
    {
    }

    /** Add global apramater with given id and initial value, optionally also description of the global parameter.
     * @param parameter_id global parameter id, must be positive integer larger than 0
     * @param value inttial parameter value
     * @param description the description
     * @return the added id
     */
    auto add_global_parameter(size_t parameter_id, float value, std::string description = "") -> size_t
    {
        if (global_parameters_map.find(parameter_id) != global_parameters_map.end()) {
            abort();
        }

        global_parameters_map.emplace(parameter_id, std::make_unique<global_param_type_t>(parameter_id, value, std::move(description)));
        return parameter_id;
    }

    auto global_parameter(size_t parameter_id) const -> const global_param_type_t&
    {
        auto the_parameter = global_parameters_map.find(parameter_id);
        if (the_parameter != global_parameters_map.end())
            return *the_parameter->second;

        throw std::runtime_error("No global parameter id=" + std::to_string(parameter_id));
    }

    auto global_parameter(size_t parameter_id) -> global_param_type_t&
    {
        return (const_cast<global_param_type_t&>(const_cast<const promille*>(this)->global_parameter(parameter_id)));
    }

    template<typename ResidualModel>
    auto make_model_planes() -> model_planes<T, ResidualModel>
    {
        return model_planes<T, ResidualModel>(this);
    }

    /** Call Mille::end()
     */
    auto end() -> void
    {
        if (verbose) {
            std::cout << "--------------------\n";
        }
        mille.end();
    }

    /** Call Mille::kill();
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
        if (!param_file) {
            return;
        }

        param_file << "Parameter\n";
        for (const auto& param : global_parameters_map) {
            if (param.second->is_anywere_free) {
                param.second->dump_pede_param(param_file);
            }
        }
    }

    auto print() const -> void
    {
        constexpr auto width = 14;
        const auto bar = std::string(width * 7, '=');
        std::cout << bar << '\n';

        std::cout << "Global parameters" << '\n';
        std::cout << bar << '\n';
        for (const auto& par : global_parameters_map) {
            std::cout << *par.second << '\n';
        }
        std::cout << bar << '\n';
    }

    auto get_mille() -> Mille& { return mille; }

    auto set_verbose(int make_verbose) -> void { verbose = make_verbose; }

  private:
    std::string mille_prefix;

    std::map<size_t, std::unique_ptr<global_param_type_t>> global_parameters_map;

    std::map<size_t, size_t> plane_indexes_map;

    Mille mille;
    int verbose {0};
};

}  // namespace promille
