#pragma once
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace math {

/******************************************************************************************/

template <class Matrix>
vec<std::size_t> cut_points(Matrix const &L, uint n) {
    boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS > g(n);
    for (auto j : range(n))
        for (auto i : range(j))
            if (L(i, j)) boost::add_edge(i, j, g);

    vec<std::size_t> points;
    boost::articulation_points(g, std::back_inserter(points));
    return points;
}

/******************************************************************************************/

struct ExponentialFitOptions {
    real tolerance = 1e-12;
    std::size_t iters = 1000;
};

// Fit a^T exp(b t) == c, in the logarithm space
template <class O, class T>
O exponential_fit(Col<T> const &a, Col<T> const &b, O const c, ExponentialFitOptions const &ops={}) {
    auto re = [](auto const &t) {
        if constexpr(is_complex<T>::value) return t.real();
        else return t;
    };
    // Initial guess
    O t = re(la::accu(-a / b)), lo = 0, up = t;
    // Double until we can get an upper bound
    for (auto iter : range(ops.iters)) {
        O const num = re(la::dot(a, la::exp(up * b)));
        if (num < c) break;
        lo = up;
        up *= 2;
    }
    // Bisection method
    for (auto iter : range(ops.iters)) {
        t = 0.5 * (lo + up);
        ASSERT(std::isfinite(t));
        O const num = re(la::dot(a, la::exp(t * b)));
        if (std::abs(num - c) < ops.tolerance) {return t;}
        
        if (num > c) lo = t;
        else up = t;
    }

    // "should maybe make hybrid newton/bisection?"
    // Col<T> e;
    // T t = la::accu(-a / b);
    // Col<T> e;
    // for (auto iter : range(ops.iters)) {
    //     e = la::exp(t * b);
    //     T const num = la::dot(a, e);
    //     T const val = std::log(num) - std::log(c);
    //     if (std::abs(val) < ops.tolerance) return t;
    //     t -= val * num / la::dot(a % b, e);
    // }
    ERROR("exponential fit via bisection did not converge", a.t(), b.t(), c, t);
}

/******************************************************************************************/

struct StationaryReferenceSystem {
    Mat<real> Hf; // weighted rate matrix (positive semidefinite)
    Mat<real> Sf; // inverse of weighted rate matrix (positive semidefinite)

    void update(Mat<real> const &H, Col<real> const &p) {
        Hf = H - la::diagmat((H * la::sqrt(p)) / la::sqrt(p));
        real scale = la::trace(Hf) / Hf.n_rows;
        Col<real> b = la::sqrt(p) / la::accu(p);
        Sf = la::inv(Hf + scale * la::outer(b, b)) - la::outer(b, b) / scale;

        // ASSERT(la::eig_sym(ef, Uf, Hf));
        // Vf = la::diagmat(la::sqrt(p)) * Uf;
        // Uf = la::diagmat(1/la::sqrt(p)) * Uf;
    }
};

/******************************************************************************************/

struct DiscreteReferenceSubnetwork {
    Mat<real> R; // rate matrix (negative definite)
    Mat<real> Z; // negative inverse of rate matrix (positive definite)
    Col<real> p; // proportional to stationary probability; held for convenience, not used in any formulas
    std::size_t m_capacity = std::numeric_limits<la::uword>::max();

    REFLECT(DiscreteReferenceSubnetwork, R, Z, p);

    DiscreteReferenceSubnetwork() = default;
    DiscreteReferenceSubnetwork(std::size_t c) : m_capacity(c) {}

    auto length() const {return R.n_cols;}
    auto capacity() const {return m_capacity;}

    void refresh() {Z = -la::inv(R);}

    Mat<real> fundamental_matrix() const {return Z;}
    Mat<real> rate_matrix() const {return R;}

    template <class F>
    void augment(real boltz, real total_rate, F &&rate) {
        REQUIRE(total_rate, >, 0);
        auto const n = length();
        p.resize(n+1);
        p(n) = boltz;
        
        R.resize(n+1, n+1);
        R(n, n) = -total_rate;
        
        for (auto i : range(n)) {
            R(n, i) = rate(n, i);
            R(i, n) = rate(i, n);
        }
        refresh();
    }

    void swap_states(uint i, uint j) {
        R.swap_rows(i, j); R.swap_cols(i, j);
        Z.swap_rows(i, j); Z.swap_cols(i, j);
        p.swap_rows(i, j);
    }
    
    void shed() {
        auto const n = length()-1;
        R.shed_col(n); R.shed_row(n);
        p.shed_row(n);
        refresh();
    }

    // O(N) Unconditional mean escape time
    real escape_time(uint i) const {return la::accu(Z.row(i));}

    // O(N^3) Conditional occupation times in the free chain (source of prototyping bug. not used!)
    Col<real> free_conditional_occupation_times(uint i, uint j) {
        Col<real> s = la::sum(R, 1);
        s(j) = 0;
        return -la::solve(R - la::diagmat(s), Col<real>(length(), la::fill::ones));
    }

    // O(N^2) Unconditional escape probabilities starting from state i
    Col<real> escape_probabilities(uint i) const {
        return -(Z.row(i).t() % la::sum(R, 1));
    }

    // O(N) Occupation times conditioned on starting at state i and escaping at state j
    Col<real> conditional_occupation_times(uint i, uint j) const {
        return Z.row(i).t() % Z.col(j) / Z(i, j);
    }

    real conditional_time(uint i, uint j) const {return la::accu(conditional_occupation_times(i, j));}

    // O(N^2) Variance of the total occupation time
    real conditional_occupation_variance(uint i, uint j) const {
        return 2 * la::as_scalar(Z.row(i) * Z * Z.col(j)) / Z(i, j);
    }

    // O(N^2) How much occupancy time is lost, starting from i and removing state j
    Col<real> cut_time_losses(uint i) const {
        return Z.row(i).t() / Z.diag() % la::sum(Z, 1);
    }

    // O(N^2) probabilities of exiting at node before hitting j
    Col<real> escape_probabilities_before(uint i, uint j) const {
        return (Z(i, j) / Z(j, j) * Z.row(j) - Z.row(i)).t() % la::sum(R, 1);
    }
};

/******************************************************************************************/

struct HittingPopulations {
    real n_incomplete=0, t_incomplete=0;
    real n_complete=0, t_complete=0;

    REFLECT(HittingPopulations, n_incomplete, t_incomplete, n_complete, t_complete);

    HittingPopulations propagate(Mat<real> const &R, Mat<real> const &Z, Mat<real> const &Z2, real sign, uint a, uint b, uint i, uint j, uint end) const {
        if (i != end && j != end) {
            QUICK_REQUIRE(i, !=, j, end);
            real phitj = (Z(a, j) * Z(j, b)) / (Z(a, b) * Z(j, j));
            real nai = sign * Z(a, i) * R(i, i);

            real ni0 = nai * (Z(i, b) * Z(j, j) - Z(i, j) * Z(j, b)) / (Z(a, b) * Z(j, j));
            real ni1 = 1 - phitj;
            real ti0 = nai * (Z2(i, b) / Z(a, b) + (Z(i, j) * Z2(j, j)) / (Z(j, j) * Z(a, j)) * phitj - (Z2(i, j) * Z(j, b) + Z(i, j) * Z2(j, b)) / (Z(a, b) * Z(j, j)));
            real ti1 = (Z2(j, j) / Z(j, j) - Z2(j, b) / Z(j, b) - Z2(a, j) / Z(a, j) + (Z(j, j) * Z2(a, b)) / (Z(a, j) * Z(j, b))) * phitj;
            real ti2 = 1 - phitj;
            real nc0 = nai * Z(i, j) / Z(a, j) * phitj;
            real nc1 = phitj;
            real tc0 = nai * (Z(j, j) * Z2(i, j) - Z(i, j) * Z2(j, j)) / (Z(j, j) * Z(a, j)) * phitj;
            real tc1 = (Z2(a, j) / Z(a, j) - Z2(j, j) / Z(j, j)) * phitj;
            real tc2 = phitj;

            return {
                .n_incomplete=ni0 + ni1 * n_incomplete,
                .t_incomplete=ti0 + ti1 * n_incomplete + ti2 * t_incomplete,
                .n_complete=n_complete + nc0 + nc1 * n_incomplete,
                .t_complete=t_complete + tc0 + tc1 * n_incomplete + tc2 * t_incomplete
            };
        }
        if (i != end) {
            real nai = sign * Z(a, i) * R(i, i);
            return {
                .n_incomplete=n_incomplete + nai * Z(i, b) / Z(a, b),
                .t_incomplete=t_incomplete + Z2(a, b) / Z(a, b) * n_incomplete + nai * Z2(i, b) / Z(a, b),
                .n_complete=n_complete,
                .t_complete=t_complete
            };
        }
        if (j != end) {
            real phitj = (Z(a, j) * Z(j, b)) / (Z(a, b) * Z(j, j));
            return {
                .n_incomplete=(1 - phitj) * n_incomplete,
                .t_incomplete=(1 - phitj) * t_incomplete + (Z2(j, j) / Z(j, j) - Z2(j, b) / Z(j, b) - Z2(a, j) / Z(a, j) + (Z(j, j) * Z2(a, b)) / (Z(a, j) * Z(j, b))) * phitj * n_incomplete,
                .n_complete=n_complete + phitj * n_incomplete,
                .t_complete=t_complete + phitj * t_incomplete + (Z2(a, j) / Z(a, j) - Z2(j, j) / Z(j, j)) * phitj * n_incomplete
            };
        }
        return {
            .n_incomplete=n_incomplete,
            .t_incomplete=t_incomplete + Z2(a, b) / Z(a, b) * n_incomplete,
            .n_complete=n_complete,
            .t_complete=t_complete
        };
    }
};

/******************************************************************************************/


struct DiscreteReferenceSquareSubnetwork : DiscreteReferenceSubnetwork {
    using base_type = DiscreteReferenceSubnetwork;
    Mat<real> Z2;

    DiscreteReferenceSquareSubnetwork() = default;

    explicit DiscreteReferenceSquareSubnetwork(std::size_t capacity) : base_type(capacity) {}

    Mat<real> squared_fundamental_matrix() const {return Z2;}

    void swap_states(uint i, uint j) {
        Z2.swap_cols(i, j); Z2.swap_rows(i, j);
        base_type::swap_states(i, j);
    }

    void refresh() {Z2 = Z * Z;}

    void shed() {base_type::shed(); refresh();}
    
    template <class ...Ts>
    void augment(Ts &&...ts) {base_type::augment(fw<Ts>(ts)...); refresh();}

    HittingPopulations propagate(HittingPopulations const &p, uint a, uint b, uint i, uint j, uint end) const {
        return p.propagate(R, Z, Z2, -1, a, b, i, j, end);
    }
};

/******************************************************************************************/

#pragma message("need to make sure about complex eigendecomposition!")
struct ContinuousReferenceSubnetwork : DiscreteReferenceSubnetwork {
    using base_type = DiscreteReferenceSubnetwork;

    using Cx = std::complex<real>;
    Mat<Cx> U, V;
    Col<Cx> e;

    EXTEND_REFLECT(ContinuousReferenceSubnetwork, base_type, U, V, e);

    void refresh() {
        ASSERT(la::eig_gen(e, U, R));
        V = la::inv(U).t();
        e = -e;
    }

    void shed() {base_type::shed(); refresh();}
    
    template <class ...Ts>
    void augment(Ts &&...ts) {base_type::augment(fw<Ts>(ts)...); refresh();}

    // O(N^2): Unconditional sample time given random float c in [0, 1]
    real sample_time(uint a, real c, ExponentialFitOptions const &ops={}) const {
        Col<Cx> alpha = (U.row(a) % la::sum(V, 0)).t();
        return exponential_fit(alpha, Col<Cx>(-e), c, ops);
    }

    // O(N): Conditional sample time given random float c in [0, 1]
    real conditional_sample_time(uint a, uint b, real c, ExponentialFitOptions const &ops={}) const {
        Col<Cx> alpha = (U.row(a).t() % V.row(b).t() / e) / Z(a, b);
        return exponential_fit(alpha, Col<Cx>(-e), c, ops);
    }

    // Escape probabilities given sample time
    Col<real> sample_escape_probabilities(uint a, real T) const {
        return la::real((V * la::diagmat(la::exp(-e * T)) * U.row(a).t()) % la::sum(R, 1) 
            / -la::as_scalar(U.row(a) * la::diagmat(e % la::exp(-e * T)) * la::sum(V, 0).t()));
    }

    // O(N^3): Occupation times given escape state and escape time
    Col<real> sample_occupation_times(uint a, uint b, real T) const {
        real const den = la::as_scalar(la::real(U.row(a) * la::diagmat(la::exp(-e * T)) * V.row(b).t()));
        Mat<Cx> Phi(length(), length());
        for (auto i : range(length())) for (auto j : range(length())) 
            Phi(i, j) = (i == j) ? (T * std::exp(-e(j) * T)) : (std::exp(-e(i) * T) - std::exp(-e(j) * T)) / (e(j) - e(i));
        return vmap<Col<real>>(range(length()), [&](auto c) {
            return la::as_scalar(la::real((U.row(a) % V.row(c)) * Phi * (U.row(c).t() % V.row(b).t()))) / den;
        });
    }
};

/******************************************************************************************/

struct DiscreteReversibleSubnetwork {
    Mat<real> R;
    Mat<real> Z;
    Col<real> p, ph;
    la::uword m_length = 0;

    REFLECT(DiscreteReversibleSubnetwork, R, Z, p, ph, m_length);

    DiscreteReversibleSubnetwork() = default;
    
    explicit DiscreteReversibleSubnetwork(la::uword capacity) :
        R(capacity, capacity), Z(capacity, capacity), p(capacity), ph(capacity) {}

    auto length() const {return m_length;}

    auto capacity() const {return R.n_rows;}

    auto all() const {return span(0, length());}

    template <class F>
    void augment_rates(real boltz, real total_rate, F &&rate) {
        QUICK_REQUIRE(total_rate, >, 0);
        QUICK_REQUIRE(boltz, >, 0);
        REQUIRE(length(), <, capacity());
        auto const b = length();

        R(b, b) = total_rate;
        p(b) = boltz;
        ph(b) = std::sqrt(boltz);
        real check = 0;
        for (auto i : range(b)) {
            check += rate(b, i);
            R(i, b)  = R(b, i) = rate(b, i) * -ph(b) / ph(i);
        }
        REQUIRE(total_rate - check, >=, -1e-5 * total_rate, total_rate, check);

        ++m_length;
    }

    void augment_inverse() {
        auto const b = length() - 1;
        if (b == 0) {
            Z(0, 0) = 1 / R(0, 0);
        } else {
            span A(0, b);
            real const r = R(b, b) - la::as_scalar(R(b, A) * Z(A, A) * R(A, b));
            Z(b, b) = 1 / r;
            Z(A, b) = Z(A, A) * R(A, b) * -Z(b, b);
            Z(b, A) = Z(A, b).t();
            Z(A, A) += r * Z(A, b) * Z(b, A);
        }
    }

    template <class F>
    void augment(real boltz, real total_rate, F &&rate) {
        augment_rates(boltz, total_rate, std::forward<F>(rate));
        augment_inverse();
    }

    void swap_states(uint i, uint j) {
        R.swap_cols(i, j); R.swap_rows(i, j);
        Z.swap_cols(i, j); Z.swap_rows(i, j);
        swap(p(i), p(j));
        swap(ph(i), ph(j));
    }

    void shed() {
        auto const b = length() - 1;
        span A(0, b);
        Z(A, A) += (-1 / Z(b, b)) * Z(A, b) * Z(b, A);
        --m_length;
    }

    Mat<real> rate_matrix() const {
        auto const A = all();
        return -la::diagmat(1/ph(A)) * R(A, A) * la::diagmat(ph(A));
    }

    Mat<real> fundamental_matrix() const {
        auto const A = all();
        return la::diagmat(1/ph(A)) * Z(A, A) * la::diagmat(ph(A));
    }

    // O(N)
    real escape_time(uint i) const {return la::dot(ph, Z.col(i)) / ph(i);}
    
    // O(N^2) Unconditional escape probabilities starting from state i
    Col<real> escape_probabilities(uint i) {
        auto const A = all();
        return Z(A, i) % (R(A, A) * ph(A)) / ph(i);
    }

    // O(N^2) 
    real conditional_occupation_variance(uint i, uint j) const {
        return 2 * la::as_scalar(Z.col(i).t() * Z * Z.col(j)) / Z(i, j);
    }

    // O(N) Occupation times conditioned on starting at state i and escaping at state j
    Col<real> conditional_occupation_times(uint i, uint j) {
        auto const A = all();
        return Z(A, i) % Z(A, j) / Z(i, j);
    }

    real conditional_time(uint i, uint j) const {
        auto const A = all();
        return la::dot(Z(A, i), Z(A, j)) / Z(i, j);
    }

    // O(N^2) How much occupancy time is lost, starting from i and removing state j
    Col<real> cut_time_losses(uint i) const {
        auto const A = all();
        return Z(i, A).t() / Z(A, A).diag() % (Z(A, A).t() * ph(A)) / std::sqrt(p(i));
    }

    Col<real> escape_probabilities_before(uint i, uint j) const {
        auto const A = all();
        return (Z(i, A) - Z(i, j) / Z(j, j) * Z(j, A)).t() % (R(A, A).t() * ph(A)) / ph(i);
    }
};

/******************************************************************************************/

struct DiscreteReversibleSquareSubnetwork : DiscreteReversibleSubnetwork {
    using base_type = DiscreteReversibleSubnetwork;
    Mat<real> Z2;

    DiscreteReversibleSquareSubnetwork() = default;

    explicit DiscreteReversibleSquareSubnetwork(std::size_t capacity) : base_type(capacity), Z2(capacity, capacity) {}

    Mat<real> squared_fundamental_matrix() const {
        auto const A = all();
        return la::diagmat(1/ph(A)) * Z2(A, A) * la::diagmat(ph(A));
    }

    void swap_states(uint i, uint j) {
        Z2.swap_cols(i, j); Z2.swap_rows(i, j);
        base_type::swap_states(i, j);
    }

    void shed() {
        auto const b = length() - 1;
        span A(0, b);
        real const s = std::sqrt(Z2(b, b)), q = s / Z(b, b), q2 = -1 / s;
        Z2(A, A) += la::outer(q * Z(A, b) + q2 * Z2(A, b), q * Z(b, A) + q2 * Z2(b, A)) 
                  + (-1 / Z2(b, b)) * la::outer(Z2(A, b), Z2(b, A));
        base_type::shed(); 
    }
    
    template <class F>
    void augment(real boltz, real total_rate, F &&rate) {
        base_type::augment_rates(boltz, total_rate, std::forward<F>(rate));
        auto const b = length() - 1;
        span A(0, b);
        if (b == 0) {
            Z2(b, b) = 1 / (R(b, b) * R(b, b));
        } else {
            Col<real> const zr = Z(A, A) * R(A, b);
            Col<real> const zr2 = Z2(A, A) * R(A, b);
            real const rm = -1 / (R(b, b) - la::dot(R(A, b), zr));
            real const kab = -rm * (1 + la::dot(R(A, b), zr2));
            Z2(A, A) += rm / kab * (zr2 * zr2.t()) - rm / kab * ((kab * zr + zr2) * (kab * zr.t() + zr2.t()));
            Z2(A, b) = rm * zr2 + kab * rm * zr;
            Z2(b, A) = Z2(A, b).t();
            Z2(b, b) = -kab * rm;
        }
        base_type::augment_inverse(); 
    }

    HittingPopulations propagate(HittingPopulations const &p, uint a, uint b, uint i, uint j, uint end) const {
        return p.propagate(R, Z, Z2, +1, a, b, i, j, end);
    }
};

/******************************************************************************************/

struct ContinuousReversibleSubnetwork : DiscreteReversibleSubnetwork {
    using base_type = DiscreteReversibleSubnetwork;
    using base_type::base_type;

    Mat<real> U, V;
    Col<real> e, me;

    EXTEND_REFLECT(ContinuousReversibleSubnetwork, base_type, U, V, e, me);

    void refresh() {
        auto const A = all();
        ASSERT(la::eig_sym(e, U, R(A, A)));
        V = U.each_col() % ph(A);
        U.each_col() /= ph(A);
        me = -e;
    }

    template <class ...Ts>
    void augment(Ts &&...ts) {
        base_type::augment(fw<Ts>(ts)...);
        refresh();
    }

    void augment() {
        base_type::shed();
        refresh();
    }

    void swap_states(uint i, uint j) {
        base_type::swap_states(i, j);
        U.swap_rows(i, j);
        V.swap_rows(i, j);
        swap(e(i), e(j));
        swap(me(i), me(j));
    }

    // O(N^2): Unconditional sample time given random float c in [0, 1]
    real sample_time(uint a, real c, ExponentialFitOptions const &ops={}) const {
        Col<real> const alpha = (U.row(a) % la::sum(V, 0)).t();
        return exponential_fit(alpha, me, c, ops);
    }

    // O(N): Conditional sample time given random float c in [0, 1]
    real conditional_sample_time(uint a, uint b, real c, ExponentialFitOptions const &ops={}) const {
        Col<real> const alpha = (U.row(a).t() % V.row(b).t() / e) * (ph(a) / ph(b) / Z(a, b));
        return exponential_fit(alpha, me, c, ops);
    }

    // Escape probabilities given sample time
    Col<real> sample_escape_probabilities(uint a, real T) const {
        auto const A = all();
        Col<real> const scale = la::exp(me * T);
        return (U.row(a) * la::diagmat(scale) * V.t()).t() % (la::diagmat(1/ ph(A)) * R(A, A) * ph(A)) 
           / la::as_scalar(U.row(a) * la::diagmat(e % scale) * la::sum(V, 0).t());
    }

    // O(N^3): Occupation times given escape state and escape time
    Col<real> sample_occupation_times(uint a, uint b, real T) const {
        Col<real> const scale = la::exp(me * T);
        Mat<real> Phi(length(), length());
        for (auto j : range(length())) for (auto i : range(length())) 
            Phi(i, j) = ((i == j) ? (T * scale(j)) : (scale(i) - scale(j)) / (e(j) - e(i))) * U(a, i) * V(b, j);
        return la::sum(U % (V * Phi.t()), 1) / la::as_scalar(U.row(a) * la::diagmat(scale) * V.row(b).t());
    }
};

/******************************************************************************************/

template <class System, class State>
struct SystemWithStates : System {
    vec<State> states;

    SystemWithStates() = default;
    explicit SystemWithStates(la::uword capacity) : System(capacity) {states.reserve(capacity);}

    void swap_states(uint i, uint j) {
        REQUIRE(i, !=, j);
        swap(states[i], states[j]);
        System::swap_states(i, j);
    }

    void shed() {
        states.pop_back();
        System::shed();
    }
};

/******************************************************************************************/

template <class System, class State, class Callback=NoOp>
real mean_time_else(System &sub, State &w, Callback &&callback={}, real const threshold=0) {
    QUICK_REQUIRE(sub.capacity(), >=, 2, "capacity must be at least 2 even for a subnetwork of 1 state");
    auto a = sub.length();
    QUICK_REQUIRE(a, ==, sub.index(w));
    sub.augment(w);

    // Pop a state if necessary
    if (threshold || sub.length() == sub.capacity()) while (true) {
        auto const losses = sub.cut_time_losses(a); 
        REQUIRE(la::accu(losses), >, 0, "Occupation times should be non-negative", losses);

        auto const least = losses.index_min();
        if (!Release) {
            // This is annoying, but for numerical precision issues sometimes cut points are chosen. Prevent them as follows:
            auto const cuts = math::cut_points(sub.R, sub.length());
            // for (auto const &c : cuts) losse.at(c) = inf<real>();
            ASSERT(!contains(cuts, least));
        }

        if (sub.length() == sub.capacity() || (threshold && losses(least) / sub.escape_time(a) < threshold)) {
            sub.swap_states(least, sub.length()-1);
            if (a == sub.length()-1) a = least;
            sub.shed();
        } else break;
    }

    auto const p_escape = sub.escape_probabilities(a); // uses I
    REQUIRE(p_escape.min(), >, -1e-8);
    REQUIRE(abs(la::accu(p_escape) - 1), <, 1e-7, "Escape probabilities should sum to 1", la::accu(p_escape), p_escape);

    // Choose the state to escape from
    uint const b = find_cumulative(p_escape, random_float()).first - p_escape.begin();

    // Calculate time until escape
    real total_time = sub.conditional_time(a, b);
    REQUIRE(total_time, >, 0, "Escape times should be positive", total_time);

    callback(sub, a, b);

    // Travel to that state
    sub.transport(w, b);
    sub.constrained_step(w, b);

    return total_time;
}

/******************************************************************************************/

}