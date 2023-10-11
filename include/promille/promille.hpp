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

template<typename T, std::size_t Nglobals, std::size_t Nlocals, typename PointType, typename VectorType, typename... ExtraArgs>
struct residual_model_base
{
    static const size_t n_globals = Nglobals;
    static const size_t n_locals = Nlocals;

    auto residual(PointType base, VectorType track) -> T
    {
        calc_derivatives(base, track);
        return calc_residual(base, track);
    }

    auto residual(PointType base, VectorType track, ExtraArgs... args) -> T
    {
        update_extras(args...);
        calc_derivatives(base, track);
        return calc_residual(base, track);
    }

    auto global_derivative(size_t variable_number) const -> T
    {
        assert(variable_number < Nglobals);
        return global_derivatives.at(variable_number);
    }

    auto local_derivative(size_t variable_number) const -> T
    {
        assert(variable_number < Nlocals);
        return local_derivatives.at(variable_number);
    }

    auto print() const -> void
    {
        for (size_t i = 0; i < Nglobals; ++i) {
            std::cout << " Global derivative " << i + 1 << " = " << global_derivative(i) << '\n';
        }

        for (size_t i = 0; i < Nlocals; ++i) {
            std::cout << " Local derivative " << i + 1 << " = " << local_derivative(i) << '\n';
        }
    }

  protected:
    template<size_t variable_number>
    auto set_global_derivative() -> T&
    {
        static_assert(variable_number < Nglobals, "Global parameter index out of bounds");
        return global_derivatives.at(variable_number);
    }

    template<size_t variable_number>
    auto set_local_derivative() -> T&
    {
        static_assert(variable_number < Nlocals, "Local parameter index out of bounds");
        return local_derivatives.at(variable_number);
    }

  protected:
    constexpr residual_model_base() = default;

  private:
    std::array<T, Nglobals> global_derivatives {0};
    std::array<T, Nlocals> local_derivatives {0};

    virtual auto update_extras(ExtraArgs... args) -> void {};
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

    auto make_state_from_parameter(Kind kind = Kind::FIXED) -> parameter_state<global_parameter<T>>
    {
        return parameter_state<global_parameter<T>>(this, kind);
    }

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
struct local_parameter
{
    using type = T;

    size_t id;
    std::string description;
    bool is_anywere_free {false};

    local_parameter(size_t id, std::string description = "")
        : id(id)
        , description(std::move(description))
    {
    }

    local_parameter(const local_parameter&) = default;
    local_parameter(local_parameter&&) = default;

    auto operator=(const local_parameter&) = delete;
    auto operator=(local_parameter&&) = delete;

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
  private:
    T* ptr;
    Kind kind;

  public:
    explicit constexpr parameter_state(T* ptr, Kind kind = Kind::FREE)
        : ptr(ptr)
        , kind(kind)
    {
    }

    parameter_state(const parameter_state&) = default;
    parameter_state(parameter_state&&) = default;

    auto operator=(const parameter_state&) = delete;
    auto operator=(parameter_state&&) = delete;

    auto get() -> T& { return *ptr; }

    auto get() const -> const T& { return *ptr; }

    auto get_kind() -> Kind { return kind; }

    auto get_kind() const -> const Kind { return kind; }

    auto set_kind(Kind value) -> void
    {
        kind = value;
        if (kind == Kind::FREE)
            ptr->is_anywere_free = true;
    }

    auto is_free() const -> bool { return kind == Kind::FREE; }

    auto operator=(Kind rhs) -> void { kind = rhs; }

    template<typename _T>
    friend auto operator<<(std::ostream& ofs, const promille::parameter_state<_T>& rhs) -> std::ostream&;
};

template<typename _T>
auto operator<<(std::ostream& ofs, const promille::parameter_state<_T>& rhs) -> std::ostream&
{
    return ofs << std::setw(6) << rhs.ptr->id << rhs.kind << std::setw(12) << rhs.ptr->value << "  " << rhs.ptr->description;
}

template<typename T, size_t Nglobals, size_t Nlocals>
struct alignment_plane
{
    static const size_t n_globals = Nglobals;
    static const size_t n_locals = Nlocals;

    std::array<parameter_state<global_parameter<T>>, Nglobals> globals;
    std::array<parameter_state<local_parameter<T>>, Nlocals> locals;

    template<typename... Ts>
    alignment_plane(std::array<parameter_state<global_parameter<T>>, Nglobals>&& g,
                    std::array<parameter_state<local_parameter<T>>, Nlocals>&& l)
        : globals(std::move(g))
        , locals(std::move(l))
    {
    }

    alignment_plane(const alignment_plane&) = default;

    template<typename... StateKind>
    auto set_globals_configuration(StateKind... kinds) -> alignment_plane<T, Nglobals, Nlocals>&
    {
        static_assert(sizeof...(kinds) == Nglobals, "Number of Kinds values do not match global parameters.");

        set_params_configuration(globals, std::make_index_sequence<sizeof...(StateKind)> {}, kinds...);
        return *this;
    }

    auto set_globals_configuration(const std::array<Kind, Nglobals>& kinds) -> alignment_plane<T, Nglobals, Nlocals>&
    {
        for (size_t i = 0; i < kinds.size(); ++i)
            globals[i].set_kind(kinds[i]);
        return *this;
    }

    template<typename... StateKind>
    auto set_locals_configuration(StateKind... kinds) -> alignment_plane<T, Nglobals, Nlocals>&
    {
        set_params_configuration(locals, std::make_index_sequence<sizeof...(StateKind)> {}, kinds...);
        return *this;
    }

    auto set_locals_configuration(const std::array<Kind, Nlocals>& kinds) -> alignment_plane<T, Nglobals, Nlocals>&
    {
        for (size_t i = 0; i < kinds.size(); ++i)
            locals[i].set_kind(kinds[i]);
        return *this;
    }

  private:
    template<typename Param, typename... StateKind, std::size_t... Is>
    auto set_params_configuration(Param& p, std::index_sequence<Is...> const&, StateKind... kinds) -> alignment_plane<T, Nglobals, Nlocals>&
    {
        ((set_param_kind(p, Is, kinds)), ...);
        return *this;
    }

    template<typename Param>
    auto set_param_kind(Param& p, size_t idx, Kind kind) -> void
    {
        p.at(idx).set_kind(kind);
    }
};

template<typename AlignmentPlane, typename ResidualModel>
class measurement_plane
{
  private:
    AlignmentPlane ptr_alignment;
    ResidualModel ptr_model;

  public:
    measurement_plane(AlignmentPlane alignment_ptr, ResidualModel residual_ptr)
        : ptr_alignment(alignment_ptr)
        , ptr_model(residual_ptr)
    {
        static_assert(AlignmentPlane::n_globals == ResidualModel::n_globals, "Number of globals parameters mismatch.");
        static_assert(AlignmentPlane::n_locals == ResidualModel::n_locals, "Number of locals parameters mismatch.");

        // std::cout << "GLOBAL STATE with alignment " << ptr_alignment << " and model " << ptr_model << '\n';
    }

    template<typename... StateKind>
    auto set_globals_configuration(StateKind... kinds) -> measurement_plane<AlignmentPlane, ResidualModel>&
    {
        static_assert(sizeof...(kinds) == AlignmentPlane::n_globals, "Number of Kinds values do not match global parameters.");

        ptr_alignment.set_globals_configuration(kinds...);
        return *this;
    }

    auto set_globals_configuration(const std::array<Kind, AlignmentPlane::n_globals>& kinds)
        -> measurement_plane<AlignmentPlane, ResidualModel>&
    {
        ptr_alignment.set_globals_configuration(kinds);
        return *this;
    }

    template<typename... StateKind>
    auto set_locals_configuration(StateKind... kinds) -> measurement_plane<AlignmentPlane, ResidualModel>&
    {
        ptr_alignment.set_locals_configuration(kinds...);
        return *this;
    }

    auto set_locals_configuration(const std::array<Kind, AlignmentPlane::n_locals>& kinds)
        -> measurement_plane<AlignmentPlane, ResidualModel>&
    {
        ptr_alignment.set_locals_configuration(kinds);
        return *this;
    }

    auto get_alignment() -> AlignmentPlane& { return ptr_alignment; }
    auto get_model() -> ResidualModel& { return ptr_model; }

    auto get_alignment() const -> const AlignmentPlane& { return ptr_alignment; }
    auto get_model() const -> const ResidualModel& { return ptr_model; }
};

template<typename ResidualModel>
class promille
{
  public:
    using main_type = float;

    using global_param_type_t = global_parameter<main_type>;
    using local_param_type_t = local_parameter<main_type>;

    using alignment_plane_t = alignment_plane<main_type, ResidualModel::n_globals, ResidualModel::n_locals>;
    using residual_model_t = ResidualModel;
    using measurement_plane_t = measurement_plane<alignment_plane_t, residual_model_t>;

    using alignment_plane_array_t = std::vector<alignment_plane_t>;
    using residual_model_array_t = std::vector<residual_model_t>;
    using measurement_plane_array_t = std::vector<measurement_plane_t>;

    using local_parameter_array_t = std::array<std::unique_ptr<local_param_type_t>, ResidualModel::n_locals>;

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
        , local_parameters_array(make_locals(std::make_index_sequence<ResidualModel::n_locals> {}))
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
        if (global_parameters_map.find(parameter_id) != global_parameters_map.end())
            abort();

        global_parameters_array.push_back(std::make_unique<global_param_type_t>(parameter_id, value, std::move(description)));
        global_parameters_map[parameter_id] = global_parameters_array.back().get();
        return parameter_id;
    }

    /** Set local parameter description.
     * @param lid local parameter id, should be in range 1 to N where N is number of the local params in the residual model
     * @param description the param description
     */
    auto set_local_parameter(size_t lid, std::string description) -> void
    {
        local_parameters_array.at(lid)->description = std::move(description);
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
        // static_assert(sizeof...(gids) == ResidualModel::n_globals, "Wrong number of global parameters ids");

        plane_indexes_map[plane_id] = measurement_plane_array.size();

        measurement_plane_array.emplace_back(
            alignment_plane_t(make_plane_globals_state(globals_ids...),
                              make_plane_locals_state(std::make_index_sequence<ResidualModel::n_locals> {})),
            residual_model_t(global_parameters_map[globals_ids]->value...));

        return measurement_plane_array.back();
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
        auto& the_plane = measurement_plane_array[plane_indexes_map[plane_id]];

        auto residuum = the_plane.get_model().residual(args...);

        std::array<float, ResidualModel::n_globals> global_derivatives;
        std::array<int, ResidualModel::n_globals> global_deriv_index;
        int global_cnt = 0;

        for (size_t i = 0; i < ResidualModel::n_globals; ++i) {
            const auto& param_state = the_plane.get_alignment().globals[i];

            if (param_state.is_free()) {
                global_derivatives[global_cnt] = the_plane.get_model().global_derivative(i);
                global_deriv_index[global_cnt] = param_state.get().id;
                global_cnt++;
            }
        }

        std::array<float, ResidualModel::n_locals> local_derivatives;
        int local_cnt = 0;

        for (size_t i = 0; i < ResidualModel::n_locals; ++i) {
            const auto& param_state = the_plane.get_alignment().locals[i];

            if (param_state.is_free()) {
                local_derivatives[local_cnt] = the_plane.get_model().local_derivative(i);
                local_cnt++;
            }
        }

        if (verbose) {
            std::cout << " --- GLOBALS ID " << plane_id << " ---\n";
        }
        if (verbose > 1) {
            the_plane.get_model().print();
        }
        if (verbose) {
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

        param_file << "Parameter\n";
        auto max_globals = global_parameters_array.size();
        for (decltype(max_globals) i = 0; i < max_globals; ++i) {
            if (global_parameters_array[i]->is_anywere_free)
                global_parameters_array[i]->dump_pede_param(param_file);
        }
    }

    auto print(bool /*use_deg*/) const -> void
    {
        constexpr auto width = 14;
        const auto bar = std::string(width * 7, '=');
        std::cout << bar << '\n';

        std::cout << "Global parameters" << '\n';
        std::cout << bar << '\n';
        for (const auto& par : global_parameters_array) {
            std::cout << *par << '\n';
        }
        std::cout << bar << '\n';

        std::cout << "Local parameters" << '\n';
        std::cout << bar << '\n';
        for (const auto& par : local_parameters_array) {
            std::cout << *par << '\n';
        }
        std::cout << bar << '\n';

        std::cout << "Measurement planes (" << plane_indexes_map.size() << ")\n";
        std::cout << bar << '\n';
        std::cout << std::right;
        for (const auto& plane_ids_it : plane_indexes_map) {
            std::cout << '[' << std::setw(5) << plane_ids_it.first << ']';
            std::cout << " GLOBAL:";

            for (const auto id : measurement_plane_array[plane_ids_it.second].get_alignment().globals) {
                std::cout << std::setw(5) << id.get().id << id.get_kind();
            }

            std::cout << "  | LOCAL:";

            for (const auto id : measurement_plane_array[plane_ids_it.second].get_alignment().locals) {
                std::cout << std::setw(5) << id.get().id << id.get_kind();
            }
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

    auto get_mille() -> Mille& { return mille; }

    auto set_verbose(int make_verbose) -> void { verbose = make_verbose; }

  private:
    template<typename... Is>
    auto make_plane_globals_state(Is... globals_ids) -> decltype(alignment_plane_t::globals)
    {
        return {global_parameters_map[globals_ids]->make_state_from_parameter()...};
    }

    template<std::size_t... Is>
    auto make_plane_locals_state(std::index_sequence<Is...> int_seq) -> decltype(alignment_plane_t::locals)
    {
        return {local_parameters_array[Is]->make_state_from_parameter()...};
    }

    template<typename Is, Is... ints>
    auto make_locals(std::integer_sequence<Is, ints...> /* unused */) -> local_parameter_array_t
    {
        return {std::make_unique<local_parameter<main_type>>(ints)...};
    }

    std::string mille_prefix;

    std::vector<std::unique_ptr<global_param_type_t>> global_parameters_array;
    std::map<size_t, global_param_type_t*> global_parameters_map;

    local_parameter_array_t local_parameters_array;

    std::map<size_t, size_t> plane_indexes_map;
    measurement_plane_array_t measurement_plane_array;

    Mille mille;
    int verbose {0};
};

}  // namespace promille
