#pragma once
#include <iostream>
#include <concepts>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <sstream>

// --- Concepts ---
template <typename T>
concept Field = requires(T a, T b) {
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a / b } -> std::convertible_to<T>;
    { -a }    -> std::convertible_to<T>;
    T(0); T(1);
};

// --- Forward Declarations (Required for Friend Templates) ---
template <Field T> class Polynomial;
template <typename T> class Monomial;

template <typename T>
std::ostream& operator<<(std::ostream& os, const Monomial<T>& m);

template <Field T>
std::ostream& operator<<(std::ostream& os, const Polynomial<T>& p);

template <Field T>
Polynomial<T> s_polynomial(const Polynomial<T>& f, const Polynomial<T>& g);

template <Field T>
Polynomial<T> multivariate_division(Polynomial<T> f, const std::vector<Polynomial<T>>& G);

template <Field T>
Polynomial<T> gcd(const Polynomial<T>& a, const Polynomial<T>& b);

template <Field T>
Polynomial<T> subresultant_gcd(const Polynomial<T>& a, const Polynomial<T>& b);

class Symbol;

// --- Helper Structs ---
template <typename T>
struct VarBinding {
    size_t id;
    T value;
};

template <typename Poly>
struct PolyBinding {
    size_t id;
    Poly value; 
};

// --- Monomial Class ---
template <typename T> 
class Monomial {
private:
    std::vector<int> exponents; 
public:
    int exponent_of(size_t var_id) const {
        if (var_id >= exponents.size()) {
            throw std::runtime_error("Monomial class exponent index out of bound.");
        }
        return exponents[var_id];
    }

    Monomial() = default;
    explicit Monomial(size_t var_id, int exponent = 1);

    Monomial operator*(const Monomial& other) const;
    Monomial operator/(const Monomial& other) const;

    std::pair<T, Monomial<T>> eval_partial(const std::vector<VarBinding<T>>& substitutions) const;

    bool is_divisible_by(const Monomial& other) const;
    int degree() const; 

    bool operator<(const Monomial& other) const; 
    bool operator==(const Monomial& other) const;
    
    // FIX: Friend declaration now specifies it is a template
    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const Monomial<U>& m);
};

// --- Polynomial Class ---
template <Field T>
class Polynomial {
private:
    struct Term {
        Monomial<T> monomial;
        T coefficient;
        bool operator<(const Term& other) const { return monomial < other.monomial; }
    };
    std::vector<Term> terms;

    void canonicalize();

public:
    Polynomial() = default;
    Polynomial(T scalar); 
    explicit Polynomial(size_t var_id);

    ~Polynomial() = default;
    Polynomial(const Polynomial& other) = default;
    Polynomial(Polynomial&& other) noexcept = default;
    Polynomial& operator=(const Polynomial& other) = default;
    Polynomial& operator=(Polynomial&& other) noexcept = default;

    Polynomial& operator+=(const Polynomial& other);
    Polynomial& operator-=(const Polynomial& other);
    Polynomial& operator*=(const Polynomial& other);
    Polynomial& operator/=(const Polynomial& other);

    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;
    Polynomial operator/(const Polynomial& other) const;
    
    Polynomial operator+(T scalar) const;
    Polynomial operator-(T scalar) const;
    Polynomial operator*(T scalar) const;
    Polynomial operator/(T scalar) const;
    Polynomial operator^(int exp) const;

    template <typename... Args>
    Polynomial eval(Args... args) const {
        std::vector<VarBinding<T>> substitutions = { args... };    
        Polynomial result; 
        return result;
    }

    template <typename... Args>
    void substitute(Args... args) {
        std::vector<VarBinding<T>> substitutions = { args... };
    }
    
    template <typename... Args>
    Polynomial compose(Args... args) const {
        std::vector<PolyBinding<Polynomial<T>>> substitutions = { args... };
        Polynomial result; 
        return result; 
    }
    
    template <Field U>
    friend Polynomial<U> s_polynomial(const Polynomial<U>& f, const Polynomial<U>& g);
    
    template <Field U>
    friend Polynomial<U> multivariate_division(Polynomial<U> f, const std::vector<Polynomial<U>>& G);
    
    template <Field U>
    friend Polynomial<U> gcd(const Polynomial<U>& a, const Polynomial<U>& b);
    
    template <Field U>
    friend Polynomial<U> subresultant_gcd(const Polynomial<U>& a, const Polynomial<U>& b);

    template <Field U>
    friend std::ostream& operator<<(std::ostream& os, const Polynomial<U>& p);

    Polynomial derivative(const Symbol& sd) const;
    Polynomial integrate(const Symbol& s) const;

    bool is_zero() const;
    T lead_coefficient() const;      
    Monomial<T> lead_monomial() const; 
};

// --- Symbol Class ---
class Symbol {
private:
    inline static std::vector<std::string> registry;
    //inline static size_t global_id_counter = 0;
    std::string name;
    size_t id;

public:
    Symbol(std::string n) : name(std::move(n)) {
        id = registry.size();
        registry.emplace_back(name);
    }

    static const std::string& get_name(size_t id) {
        if (id >= registry.size()) {
            throw std::runtime_error("Symbol::get_name id out of bound.");
        }
        return registry[id];
    }

    template <typename T>
    VarBinding<T> operator=(T val) const {
        return VarBinding<T>{id, val};
    }

    template <typename T>
    PolyBinding<Polynomial<T>> operator=(const Polynomial<T>& val) const {
        return PolyBinding<Polynomial<T>>{id, val};
    }

    template <Field T>
    operator Polynomial<T>() const;

    size_t get_id() const { return id; }
};

template <Field T>
Symbol::operator Polynomial<T>() const {
    return Polynomial<T>(this->id);
}

template <typename T>
Monomial<T>::Monomial(size_t var_id, int exponent) {
    if (exponent < 0) {
        throw std::runtime_error("Monomial exponents must be nonnegative.");
    }
    exponents.resize(var_id + 1, 0);
    exponents[var_id] = exponent;
}

template <typename T>
Monomial<T> Monomial<T>::operator*(const Monomial &other) const {
    Monomial result;
    size_t max_size = std::max(exponents.size(), other.exponents.size());
    result.exponents.resize(max_size, 0);

    for (size_t i = 0; i < max_size; i++) {
        int e = (i < exponents.size()) ? exponents[i] : 0;
        int f = (i < other.exponents.size()) ? other.exponents[i] : 0;
        result.exponents[i] = e + f;
    }
    return result;
}

template <typename T>
Monomial<T> Monomial<T>::operator/(const Monomial& other) const {
    Monomial result;
    size_t max_size = std::max(exponents.size(), other.exponents.size());
    result.exponents.resize(max_size, 0);

    for (size_t i = 0; i < max_size; i++) {
        int e = (i < exponents.size()) ? exponents[i] : 0;
        int f = (i < other.exponents.size()) ? other.exponents[i] : 0;
        if (e < f) {
            throw std::runtime_error("Monomial division failed: result not a polynomial.");
        }
        result.exponents[i] = e - f;
    }
    return result;
}

template <typename T>
std::pair<T, Monomial<T>> Monomial<T>::eval_partial(const std::vector<VarBinding<T>>& substitutions) const {
    Monomial result = *this;
    T coefficient = T(1);

    for (const auto& binding : substitutions) {
        if (binding.id < result.exponents.size()) {
            int exp = result.exponents[binding.id];
            if (exp > 0) {
                T base = binding.value;
                for(int k = 0; k < exp; k++) coefficient = coefficient * base;
                result.exponents[binding.id] = 0; 
            }
        }
    }
    return {coefficient, result};
}

template <typename T>
bool Monomial<T>::is_divisible_by(const Monomial& other) const {
    size_t size = other.exponents.size();
    if (exponents.size() < size) {
        for (size_t i = exponents.size(); i < size; i++) {
            if (other.exponents[i] > 0) return false;
        }
        size = exponents.size();
    }
    
    for (size_t i = 0; i < size; i++) {
        if (exponents[i] < other.exponents[i]) return false;
    }
    return true;
}

template <typename T>
int Monomial<T>::degree() const {
    int sum = 0;
    for (int e : exponents) sum += e;
    return sum;
}

// --- Comparison (Graded Reverse Lexicographic - GrevLex) ---
// Returns true if 'this' < 'other'
template <typename T>
bool Monomial<T>::operator<(const Monomial& other) const {
    int deg_this = this->degree();
    int deg_that = other.degree();
    
    if (deg_this != deg_that) {
        return deg_this < deg_that;
    }
    size_t max_size = std::max(exponents.size(), other.exponents.size());

    for (size_t i = max_size; i-- > 0; ) {
        int e = (i < exponents.size()) ? exponents[i] : 0;
        int f = (i < other.exponents.size()) ? other.exponents[i] : 0;
        if (e != f) {
            return e > f; 
        }
    }

    return false; // They are strictly equal
}

template <typename T>
bool Monomial<T>::operator==(const Monomial& other) const {
    size_t max_size = std::max(exponents.size(), other.exponents.size());
    for (size_t i = 0; i < max_size; i++) {
        int e = (i < exponents.size()) ? exponents[i] : 0;
        int f = (i < other.exponents.size()) ? other.exponents[i] : 0;
        if (e != f) return false;
    }
    return true;
}

template <typename U>
std::ostream& operator<<(std::ostream& os, const Monomial<U>& m) {
    bool first = true;
    bool is_one = true;

    for (size_t i = 0; i < m.exponents.size(); ++i) {
        if (m.exponents[i] > 0) {
            if (!first) os << "*";
            os << Symbol::get_name(i); 
            if (m.exponents[i] > 1) {
                os << "^" << m.exponents[i];
            }
            first = false;
            is_one = false;
        }
    }
    if (is_one) os << "1";
    return os;
}