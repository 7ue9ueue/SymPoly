#include <iostream>
#include <concepts>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

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
    Polynomial operator^(unsigned int exp) const;

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
    
    // FIX: Use different template parameters (U) to resolve warnings
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
    inline static size_t global_id_counter = 0;
    std::string name;
    size_t id;

public:
    Symbol(std::string n) : name(std::move(n)), id(global_id_counter++) {}

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

int main() {
    return 0;
}