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

template<typename T, size_t Nlocals>
struct residual_model_base
{
    const size_t n_globals;
    static const size_t n_locals = Nlocals;

    virtual auto residual() const -> T = 0;

    auto global_derivative(size_t variable_number) const -> T
    {
        assert(variable_number < n_globals);
        return global_derivatives.at(variable_number);
    }

    auto local_derivative(size_t variable_number) const -> T
    {
        assert(variable_number < n_locals);
        return local_derivatives.at(variable_number);
    }

    auto print() const -> void
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

    residual_model_base() = delete;

    residual_model_base(size_t nglobals)
        : n_globals(nglobals)
        , global_derivatives(n_globals)
    {
    }

  private:
    std::vector<T> global_derivatives;
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

    auto make_state_from_parameter(Kind kind = Kind::FIXED) -> parameter_state<T> { return parameter_state<T>(this, kind); }

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
struct parameter_state
{
  private:
    global_parameter<T>* ptr;
    Kind kind;

  public:
    explicit constexpr parameter_state(global_parameter<T>* parameter_ptr, Kind parameter_kind = Kind::FREE)
        : ptr(parameter_ptr)
        , kind(parameter_kind)
    {
    }

    parameter_state(const parameter_state&) = default;
    parameter_state(parameter_state&&) = default;

    auto operator=(const parameter_state&) = delete;
    auto operator=(parameter_state&&) = delete;

    auto get() -> global_parameter<T>& { return *ptr; }

    auto get() const -> const global_parameter<T>& { return *ptr; }

    auto get_kind() -> Kind { return kind; }

    auto get_kind() const -> Kind { return kind; }

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

template<typename T, size_t Nlocals>
struct alignment_plane
{
    const size_t n_globals;
    const size_t n_locals;

    std::vector<parameter_state<T>> globals;
    std::array<Kind, Nlocals> locals;

    template<typename... GlobalsIds>
    alignment_plane(GlobalsIds... globals_ids)
        : n_globals(sizeof...(globals_ids))
        , n_locals(Nlocals)
        , globals({globals_ids...})
    {
    }

    alignment_plane(const alignment_plane&) = default;

    template<typename... StateKind>
    auto set_globals_configuration(StateKind... kinds) -> alignment_plane<T, Nlocals>&
    {
        // static_assert(sizeof...(kinds) == Nglobals, "Number of Kinds values do not match global parameters.");

        assert(sizeof...(kinds) == n_globals);

        set_params_configuration(globals, std::make_index_sequence<sizeof...(StateKind)> {}, kinds...);
        return *this;
    }

    auto set_globals_configuration(const std::vector<Kind>& kinds) -> alignment_plane<T, Nlocals>&
    {
        assert(kinds.size() == n_globals);

        for (size_t i = 0; i < kinds.size(); ++i)
            globals[i].set_kind(kinds[i]);
        return *this;
    }

    template<typename... StateKind>
    auto set_locals_configuration(StateKind... kinds) -> alignment_plane<T, Nlocals>&
    {
        // set_params_configuration(locals, std::make_index_sequence<sizeof...(StateKind)> {}, kinds...);
        locals = {kinds...};
        return *this;
    }

    auto set_locals_configuration(const std::array<Kind, Nlocals>& kinds) -> alignment_plane<T, Nlocals>&
    {
        locals = kinds;
        return *this;
    }

  private:
    template<typename Param, typename... StateKind, std::size_t... Is>
    auto set_params_configuration(Param& p, std::index_sequence<Is...> const&, StateKind... kinds) -> alignment_plane<T, Nlocals>&
    {
        ((set_param_kind(p, Is, kinds)), ...);
        return *this;
    }

    template<typename Param>
    auto set_param_kind(Param& p, size_t idx, Kind kind) -> void
    {
        p.at(idx) = kind;
    }
};

template<typename T, typename AlignmentPlane, std::size_t Nlocals>
class measurement_plane
{
  private:
    AlignmentPlane ptr_alignment;
    std::unique_ptr<residual_model_base<T, Nlocals>> ptr_model;

  public:
    template<typename ResidualModel>
    measurement_plane(AlignmentPlane alignment_ptr, ResidualModel model)
        : ptr_alignment(alignment_ptr)
        , ptr_model(std::make_unique<ResidualModel>(model))
    {
        assert(alignment_ptr.n_globals == model.n_globals);
    }

    template<typename... StateKind>
    auto set_globals_configuration(StateKind... kinds) -> measurement_plane<T, AlignmentPlane, Nlocals>&
    {
        assert(ptr_alignment.n_globals == sizeof...(kinds));

        ptr_alignment.set_globals_configuration(kinds...);
        return *this;
    }

    template<std::size_t Nglobals>
    auto set_globals_configuration(const std::array<Kind, Nglobals>& kinds) -> measurement_plane<T, AlignmentPlane, Nlocals>&
    {
        assert(Nglobals == ptr_alignment.globals.size());

        for (size_t i = 0; i < Nglobals; ++i)
            ptr_alignment.globals[i] = kinds[i];

        return *this;
    }

    template<typename... StateKind>
    auto set_locals_configuration(StateKind... kinds) -> measurement_plane<T, AlignmentPlane, Nlocals>&
    {
        ptr_alignment.set_locals_configuration(kinds...);
        return *this;
    }

    auto set_locals_configuration(const std::array<Kind, Nlocals>& kinds) -> measurement_plane<T, AlignmentPlane, Nlocals>&
    {
        ptr_alignment.set_locals_configuration(kinds);
        return *this;
    }

    auto get_alignment() -> AlignmentPlane& { return ptr_alignment; }
    auto get_model() -> residual_model_base<T, Nlocals>* { return ptr_model.get(); }

    auto get_alignment() const -> const AlignmentPlane& { return ptr_alignment; }
    auto get_model() const -> const residual_model_base<T, Nlocals>* { return ptr_model.get(); }
};

template<std::size_t Nlocals>
class promille
{
  public:
    using main_type = float;
    using global_param_type_t = global_parameter<main_type>;
    using alignment_plane_t = alignment_plane<main_type, Nlocals>;
    using measurement_plane_t = measurement_plane<main_type, alignment_plane_t, Nlocals>;
    using measurement_plane_array_t = std::vector<measurement_plane_t>;

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
    // auto set_local_parameter(size_t lid, std::string description) -> void
    // {
    //     local_parameters_array.at(lid)->description = std::move(description);
    // }

    /**
     * Add measurement plane global parameters configuration.
     *
     * @param plane_id the measurement plane id
     * @param gids list of global parameters ids associated to this measurement plane
     * @return measurement plane object
     */

    template<typename ResidualModel, typename... GlobalIds>
    auto add_plane(size_t plane_id, GlobalIds... globals_ids) -> measurement_plane_t&
    {
        static_assert(Nlocals == ResidualModel::n_locals, "Wrong number of local parameters in model");

        plane_indexes_map[plane_id] = measurement_plane_array.size();

        measurement_plane_array.emplace_back(alignment_plane_t(global_parameters_map[globals_ids]->make_state_from_parameter()...),
                                             ResidualModel(global_parameters_map[globals_ids]->value...));

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
    auto add_measurement(size_t plane_id, float sigma) -> void
    {
        auto& the_plane = measurement_plane_array[plane_indexes_map[plane_id]];

        auto model = the_plane.get_model();
        auto residuum = model->residual();

        std::vector<float> global_derivatives(model->n_globals);
        std::vector<int> global_deriv_index(model->n_globals);

        int global_cnt = 0;

        for (size_t i = 0; i < model->n_globals; ++i) {
            const auto& param_state = the_plane.get_alignment().globals[i];

            if (param_state.is_free()) {
                global_derivatives[global_cnt] = the_plane.get_model()->global_derivative(i);
                global_deriv_index[global_cnt] = param_state.get().id;
                global_cnt++;
            }
        }

        std::vector<float> local_derivatives(model->n_locals);
        int local_cnt = 0;

        for (size_t i = 0; i < model->n_locals; ++i) {
            const auto& param_state = the_plane.get_alignment().locals[i];

            if (param_state == Kind::FREE) {
                local_derivatives[local_cnt] = the_plane.get_model()->local_derivative(i);
                local_cnt++;
            }
        }

        if (verbose) {
            std::cout << " --- GLOBALS ID " << plane_id << " ---\n";
        }
        if (verbose > 1) {
            the_plane.get_model()->print();
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
        mille.mille(local_cnt,
                    local_derivatives.data(),
                    global_cnt,
                    global_derivatives.data(),
                    global_deriv_index.data(),
                    static_cast<float>(residuum),
                    sigma);
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

            for (const auto kind : measurement_plane_array[plane_ids_it.second].get_alignment().locals) {
                std::cout << std::setw(5) << kind;
            }
            std::cout << '\n';
            measurement_plane_array[plane_ids_it.second].get_model()->print();
        }
        std::cout << std::left;
        std::cout << bar << '\n';
    }

    auto get_mille() -> Mille& { return mille; }

    auto set_verbose(int make_verbose) -> void { verbose = make_verbose; }

  private:
    template<typename... Is>
    auto make_plane_globals_state(Is... globals_ids) -> decltype(alignment_plane_t::globals)
    {
        return {global_parameters_map[globals_ids]->make_state_from_parameter()...};
    }

    std::string mille_prefix;

    std::vector<std::unique_ptr<global_param_type_t>> global_parameters_array;
    std::map<size_t, global_param_type_t*> global_parameters_map;

    std::map<size_t, size_t> plane_indexes_map;
    measurement_plane_array_t measurement_plane_array;

    Mille mille;
    int verbose {0};
};

}  // namespace promille
