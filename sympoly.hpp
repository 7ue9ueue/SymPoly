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

// --- Monomial Implementation ---
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

// --- Polynomial Implementation ---

template <Field T>
void Polynomial<T>::canonicalize() {
    if (terms.empty()) return;

    // Sort terms based on Monomial ordering (Strict Weak Ordering)
    // Note: Since Monomial::operator< defines the order, we rely on it.
    std::sort(terms.begin(), terms.end());

    // Merge duplicate monomials and sum coefficients
    size_t write_index = 0;
    for (size_t read_index = 1; read_index < terms.size(); ++read_index) {
        if (terms[write_index].monomial == terms[read_index].monomial) {
            terms[write_index].coefficient = terms[write_index].coefficient + terms[read_index].coefficient;
        } else {
            ++write_index;
            if (write_index != read_index) {
                terms[write_index] = terms[read_index];
            }
        }
    }
    
    // Resize to remove merged elements
    terms.resize(write_index + 1);

    // 3. Remove terms with zero coefficients (Sparsity)
    // Using erase-remove idiom
    terms.erase(std::remove_if(terms.begin(), terms.end(), 
        [](const Term& term) {
            // Assuming T(0) constructs the additive identity
            // For floating point, strict equality might be dangerous, 
            // but for exact fields (PolySolve++ goal), this is correct.
            return term.coefficient == T(0); 
        }), 
    terms.end());
}

template <Field T>
Polynomial<T>::Polynomial(T scalar) {
    if (scalar != T(0)) {
        terms.push_back({Monomial<T>(), scalar});
    }
}

template <Field T>
Polynomial<T>::Polynomial(size_t var_id) {
    terms.push_back({Monomial<T>(var_id, 1), T(1)});
}

template <Field T>
bool Polynomial<T>::is_zero() const {
    return terms.empty();
}

template <Field T>
T Polynomial<T>::lead_coefficient() const {
    if (terms.empty()) return T(0);
    // Terms are sorted ascending by default in Monomial < operator (usually graded reverse lex).
    // The "Largest" monomial is at the back.
    return terms.back().coefficient;
}

template <Field T>
Monomial<T> Polynomial<T>::lead_monomial() const {
    if (terms.empty()) return Monomial<T>();
    return terms.back().monomial;
}

// --- Arithmetic Operators ---

template <Field T>
Polynomial<T> Polynomial<T>::operator+(const Polynomial& other) const {
    Polynomial result;
    result.terms.reserve(terms.size() + other.terms.size());

    // Merge sorted sequences
    auto this_it = terms.begin();
    auto that_it = other.terms.begin();

    while (this_it != terms.end() && that_it != other.terms.end()) {
        if (this_it->monomial < that_it->monomial) {
            result.terms.push_back(*this_it);
            ++this_it;
        } else if (that_it->monomial < this_it->monomial) {
            result.terms.push_back(*that_it);
            ++that_it;
        } else {
            // Monomials are equal, add coefficients
            T sum_coefficient = this_it->coefficient + that_it->coefficient;
            if (sum_coefficient != T(0)) {
                result.terms.push_back({this_it->monomial, sum_coefficient});
            }
            ++this_it;
            ++that_it;
        }
    }

    // Append remaining elements
    while (this_it != terms.end()) {
        result.terms.push_back(*this_it++);
    }
    while (that_it != other.terms.end()) {
        result.terms.push_back(*that_it++);
    }

    // No need to canonicalize as we merged strictly sorted lists
    return result;
}

template <Field T>
Polynomial<T> Polynomial<T>::operator-(const Polynomial& other) const {
    // Optimization: P - Q is P + (-1 * Q)
    // But direct implementation avoids copying 'other' just to negate it.
    Polynomial result;
    result.terms.reserve(terms.size() + other.terms.size());

    auto this_it = terms.begin();
    auto that_it = other.terms.begin();

    while (this_it != terms.end() && that_it != other.terms.end()) {
        if (this_it->monomial < that_it->monomial) {
            result.terms.push_back(*this_it);
            ++this_it;
        } else if (that_it->monomial < this_it->monomial) {
            // Subtracting the other term means negating coefficient
            result.terms.push_back({that_it->monomial, -that_it->coefficient});
            ++that_it;
        } else {
            T diff_coefficient = this_it->coefficient - that_it->coefficient;
            if (diff_coefficient != T(0)) {
                result.terms.push_back({this_it->monomial, diff_coefficient});
            }
            ++this_it;
            ++that_it;
        }
    }

    while (this_it != terms.end()) {
        result.terms.push_back(*this_it++);
    }
    while (that_it != other.terms.end()) {
        result.terms.push_back({that_it->monomial, -that_it->coefficient});
        ++that_it;
    }
    return result;
}

template <Field T>
Polynomial<T> Polynomial<T>::operator*(const Polynomial& other) const {
    if (this->is_zero() || other.is_zero()) return Polynomial<T>();

    Polynomial result;
    // Heuristic reservation
    result.terms.reserve(terms.size() * other.terms.size());

    for (const auto& term_lhs : terms) {
        for (const auto& term_rhs : other.terms) {
            Monomial<T> product_monomial = term_lhs.monomial * term_rhs.monomial;
            T product_coefficient = term_lhs.coefficient * term_rhs.coefficient;
            result.terms.push_back({product_monomial, product_coefficient});
        }
    }

    result.canonicalize();
    return result;
}

template <Field T>
Polynomial<T> Polynomial<T>::operator/(const Polynomial& other) const {
    // Currently strictly supporting scalar division via this operator
    // Full polynomial division (Euclidean/Groebner) is handled by multivariate_division friend.
    if (other.terms.size() == 1 && other.terms[0].monomial == Monomial<T>()) {
        return *this / other.terms[0].coefficient;
    }
    throw std::runtime_error("Direct Polynomial/Polynomial division not supported via operator/. Use multivariate_division or pseudo-division.");
}

// --- Scalar Operators ---

template <Field T>
Polynomial<T> Polynomial<T>::operator+(T scalar) const {
    Polynomial result = *this;
    // We can optimize this by inserting into sorted position, but this is safe:
    result.terms.push_back({Monomial<T>(), scalar}); 
    result.canonicalize();
    return result;
}

template <Field T>
Polynomial<T> Polynomial<T>::operator-(T scalar) const {
    return *this + (-scalar);
}

template <Field T>
Polynomial<T> Polynomial<T>::operator*(T scalar) const {
    if (scalar == T(0)) return Polynomial<T>();
    
    Polynomial result = *this;
    for (auto& term : result.terms) {
        term.coefficient = term.coefficient * scalar;
    }
    return result;
}

template <Field T>
Polynomial<T> Polynomial<T>::operator/(T scalar) const {
    if (scalar == T(0)) throw std::runtime_error("Division by zero.");
    
    // In a Field, division is multiplication by inverse.
    // However, we use the / operator provided by the type T.
    Polynomial result = *this;
    for (auto& term : result.terms) {
        term.coefficient = term.coefficient / scalar;
    }
    return result;
}

// --- Assignment Operators ---

template <Field T>
Polynomial<T>& Polynomial<T>::operator+=(const Polynomial& other) {
    *this = *this + other;
    return *this;
}

template <Field T>
Polynomial<T>& Polynomial<T>::operator-=(const Polynomial& other) {
    *this = *this - other;
    return *this;
}

template <Field T>
Polynomial<T>& Polynomial<T>::operator*=(const Polynomial& other) {
    *this = *this * other;
    return *this;
}

template <Field T>
Polynomial<T>& Polynomial<T>::operator/=(const Polynomial& other) {
    *this = *this / other;
    return *this;
}

// --- Exponentiation (Binary Exponentiation) ---

template <Field T>
Polynomial<T> Polynomial<T>::operator^(int exp) const {
    if (exp < 0) throw std::runtime_error("Negative polynomial exponentiation not supported.");
    if (exp == 0) return Polynomial<T>(T(1));
    if (exp == 1) return *this;

    Polynomial<T> base = *this;
    Polynomial<T> result(T(1));
    
    while (exp > 0) {
        if (exp % 2 == 1) {
            result *= base;
        }
        base *= base;
        exp /= 2;
    }
    return result;
}

// --- Calculus ---

template <Field T>
Polynomial<T> Polynomial<T>::derivative(const Symbol& symbol) const {
    Polynomial result;
    result.terms.reserve(terms.size());
    size_t var_id = symbol.get_id();
    
    // Helper unit monomial for division to reduce degree
    // d/dx (x^n) -> n * x^(n-1)
    // x^(n-1) is achieved by x^n / x^1
    Monomial<T> var_unit(var_id, 1);

    for (const auto& term : terms) {
        int exponent = term.monomial.exponent_of(var_id);
        
        if (exponent > 0) {
            T new_coefficient = term.coefficient * T(exponent);
            
            // Technically unsafe if Monomial division fails, but here we guarantee strict divisibility
            // because exponent > 0 implies divisibility by x^1.
            Monomial<T> new_monomial = term.monomial / var_unit;
            
            result.terms.push_back({new_monomial, new_coefficient});
        }
        // If exponent is 0, the term is constant wrt this variable, derivative is 0 (omitted)
    }
    
    // Result is naturally sorted because strictly decreasing degree in one var 
    // preserves order of other vars in GrevLex usually, but to be safe:
    result.canonicalize(); 
    return result;
}

template <Field T>
Polynomial<T> Polynomial<T>::integrate(const Symbol& symbol) const {
    Polynomial result;
    result.terms.reserve(terms.size());
    size_t var_id = symbol.get_id();

    // Helper unit monomial to increase degree
    Monomial<T> var_unit(var_id, 1);

    for (const auto& term : terms) {
        int old_exponent = term.monomial.exponent_of(var_id);
        int new_exponent = old_exponent + 1;

        // Integration Rule: x^n -> x^(n+1) / (n+1)
        Monomial<T> new_monomial = term.monomial * var_unit;
        T div_scalar = T(new_exponent);
        T new_coefficient = term.coefficient / div_scalar;

        result.terms.push_back({new_monomial, new_coefficient});
    }

    result.canonicalize();
    return result;
}

// --- IO Stream ---

template <Field T>
std::ostream& operator<<(std::ostream& os, const Polynomial<T>& p) {
    if (p.is_zero()) {
        os << "0";
        return os;
    }

    bool first = true;
    // Iterate in reverse because terms are sorted ascending (smallest first), 
    // but humans read polynomials largest degree first.
    for (auto it = p.terms.rbegin(); it != p.terms.rend(); ++it) {
        T coeff = it->coefficient;
        
        if (!first) {
            // Handle sign for formatting
            // Assuming T supports comparison with 0
            // If T is generic (like complex), this check might be simplistic, 
            // but standard for Real fields.
            // Note: We use + for everything unless we inspect the type, 
            // but checking < 0 makes output cleaner.
            os << " + "; 
        }
        
        // Print Coefficient
        // Logic: Print coeff if it's not 1, OR if the monomial is degree 0 (constant term)
        bool is_constant_term = (it->monomial.degree() == 0);
        
        // Special printing logic for 1/-1 to avoid "1*x" or "-1*x"
        // This requires T(1) comparison which might be costly for BigInts, but necessary for clean output.
        if (is_constant_term) {
             os << coeff;
        } else {
             if (coeff == T(-1)) os << "-";
             else if (coeff != T(1)) os << coeff << "*";
        }

        // Print Monomial
        if (!is_constant_term) {
            os << it->monomial;
        }
        
        first = false;
    }
    return os;
}