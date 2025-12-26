#include <iostream>
#include <cassert>
#include <vector>
#include "sympoly.hpp"

struct TestSuiteResults {
    int passed = 0;
    int total = 0;
    void record(bool success) {
        total++;
        if (success) passed++;
    }
};
TestSuiteResults global_results;

#define ASSERT_TEST(cond) \
    do { \
        if (!(cond)) { \
            std::cerr << "Assertion failed: " << #cond << " at line " << __LINE__ << "\n"; \
            global_results.record(false); \
            return; \
        } \
    } while (0)

// --------------------------------------------------------------------------
// Test Modules
// --------------------------------------------------------------------------

void test_polynomial_eval_substitute() {
    std::cout << "[Test] Eval & Substitute (10 cases)... ";
    
    Symbol<double> x("x"), y("y");
    Polynomial<double> p_xy = x*x + y + 5.0; // x^2 + y + 5

    // 1. Full Evaluation to scalar
    auto res1 = p_xy.eval(x = 2.0, y = 3.0); // 4 + 3 + 5 = 12
    ASSERT_TEST(res1.is_constant() && res1.to_scalar() == 12.0);

    // 2. Partial Evaluation (x only)
    auto res2 = p_xy.eval(x = 3.0); // 9 + y + 5 = y + 14
    ASSERT_TEST(!res2.is_constant());
    ASSERT_TEST(res2.eval(y = 1.0).to_scalar() == 15.0);

    // 3. Partial Evaluation (y only)
    auto res3 = p_xy.eval(y = -5.0); // x^2 - 5 + 5 = x^2
    ASSERT_TEST(res3.lead_monomial().degree() == 2);
    ASSERT_TEST(res3.lead_coefficient() == 1.0);

    // 4. Eval on zero polynomial
    Polynomial<double> zero;
    ASSERT_TEST(zero.eval(x = 100.0).is_zero());

    // 5. Eval on constant polynomial
    Polynomial<double> con(42.0);
    ASSERT_TEST(con.eval(x = 0.0).to_scalar() == 42.0);

    // 6. Substitution (In-place)
    Polynomial<double> p_sub = x + 10.0;
    p_sub.substitute(x = 5.0);
    ASSERT_TEST(p_sub.is_constant() && p_sub.to_scalar() == 15.0);

    // 7. Evaluating variable not in polynomial
    Polynomial<double> p_x = x * 2.0;
    auto res7 = p_x.eval(y = 10.0); // Should remain 2x
    ASSERT_TEST(res7.lead_coefficient() == 2.0);

    // 8. Evaluation with zero
    ASSERT_TEST(p_xy.eval(x = 0.0, y = 0.0).to_scalar() == 5.0);

    // 9. Large power evaluation
    Polynomial<double> p_pwr = (x ^ 4); 
    ASSERT_TEST(p_pwr.eval(x = 2.0).to_scalar() == 16.0);

    // 10. Multiple Substitutions sequence
    Polynomial<double> p_seq = x + y;
    p_seq.substitute(x = 1.0);
    p_seq.substitute(y = 2.0);
    ASSERT_TEST(p_seq.to_scalar() == 3.0);

    std::cout << "Passed.\n";
    global_results.record(true);
}

void test_polynomial_compose() {
    std::cout << "[Test] Composition (10 cases)... ";

    Symbol<double> x("x"), y("y"), u("u"), v("v");
    
    // 1. Basic Linear Composition: f(x)=x+1, g(u)=u+2 => f(g(u))=u+3
    Polynomial<double> f1 = x + 1.0;
    Polynomial<double> g1 = u + 2.0;
    auto res1 = f1.compose(x = g1);
    ASSERT_TEST(res1.eval(u = 0.0).to_scalar() == 3.0);

    // 2. Quadratic Composition: f(x)=x^2, g(u)=u+1 => (u+1)^2 = u^2+2u+1
    Polynomial<double> f2 = x ^ 2;
    Polynomial<double> g2 = u + 1.0;
    auto res2 = f2.compose(x = g2);
    ASSERT_TEST(res2.lead_monomial().degree() == 2);
    ASSERT_TEST(res2.eval(u = 1.0).to_scalar() == 4.0);

    // 3. Multivariate Composition: f(x,y)=x*y, x=u+v, y=u-v => u^2 - v^2
    Polynomial<double> f3 = x * y;
    auto res3 = f3.compose(x = (u + v), y = (u - v));
    ASSERT_TEST(res3.eval(u = 3.0, v = 2.0).to_scalar() == 5.0); // 9 - 4

    // 4. Compose with Constant
    Polynomial<double> f4 = x * 2.0;
    auto res4 = f4.compose(x = Polynomial<double>(5.0));
    ASSERT_TEST(res4.is_constant() && res4.to_scalar() == 10.0);

    // 5. Identity Composition
    Polynomial<double> f5 = x*x + x;
    auto res5 = f5.compose(x = Polynomial<double>(x));
    ASSERT_TEST((res5 - f5).is_zero());

    // 6. Compose into Zero
    Polynomial<double> zero;
    ASSERT_TEST(zero.compose(x = f5).is_zero());

    // 7. Compose with Zero
    ASSERT_TEST(f5.compose(x = zero).is_zero());

    // 8. Partial Composition (Compose x, leave y)
    Polynomial<double> f8 = x + y;
    auto res8 = f8.compose(x = u*u); // u^2 + y
    ASSERT_TEST(res8.eval(u = 2.0, y = 3.0).to_scalar() == 7.0);

    // 9. Deep nested logic: f(g(h(x)))
    Polynomial<double> h = x + 1.0;
    Polynomial<double> g = x ^ 2;
    Polynomial<double> f = x - 1.0;
    auto res9 = f.compose(x = g.compose(x = h)); // ((x+1)^2) - 1
    ASSERT_TEST(res9.eval(x = 1.0).to_scalar() == 3.0);

    // 10. Self-composition: f(x) = x^2, f(f(x)) = x^4
    Polynomial<double> f10 = x ^ 2;
    auto res10 = f10.compose(x = f10);
    ASSERT_TEST(res10.lead_monomial().degree() == 4);

    std::cout << "Passed.\n";
    global_results.record(true);
}

void test_polynomial_div_mod() {
    std::cout << "[Test] Div & Mod (10 cases)... ";

    Symbol<double> x("x"), y("y");

    // 1. Exact Monomial Division: x^2 / x = x
    Polynomial<double> p1 = x^2;
    auto [q1, r1] = p1.div_mod(x);
    ASSERT_TEST(q1.lead_monomial().degree() == 1 && r1.is_zero());

    // 2. Division with Remainder: (x^2 + 1) / x = x, rem 1
    Polynomial<double> p2 = (x^2) + 1.0;
    auto [q2, r2] = p2.div_mod(x);
    ASSERT_TEST(q2.lead_monomial().degree() == 1 && r2.to_scalar() == 1.0);

    // 3. Difference of Squares: (x^2 - 1) / (x - 1) = x + 1
    Polynomial<double> p3 = (x^2) - 1.0;
    auto [q3, r3] = p3.div_mod(x - 1.0);
    ASSERT_TEST(r3.is_zero());
    ASSERT_TEST(q3.eval(x = 0.0).to_scalar() == 1.0);

    // 4. Division by Scalar (Field property)
    Polynomial<double> p4 = x * 4.0;
    auto [q4, r4] = p4.div_mod(Polynomial<double>(2.0));
    ASSERT_TEST(q4.lead_coefficient() == 2.0 && r4.is_zero());

    // 5. Divisor Degree > Dividend Degree
    Polynomial<double> p5 = x;
    auto [q5, r5] = p5.div_mod(x^2);
    ASSERT_TEST(q5.is_zero());
    ASSERT_TEST((r5 - x).is_zero());

    // 6. Property Check: Dividend = Q*D + R
    Polynomial<double> f = (x^3) + (x^2)*2.0 + x + 5.0;
    Polynomial<double> g = (x^2) + 1.0;
    auto [q6, r6] = f.div_mod(g);
    Polynomial<double> check6 = q6 * g + r6;
    ASSERT_TEST((f - check6).is_zero());

    // 7. Division by Zero (Exception)
    bool caught = false;
    try { f.div_mod(Polynomial<double>(0.0)); } catch(...) { caught = true; }
    ASSERT_TEST(caught);

    // 8. Multivariate Reduction (Simple): (x + y) / x
    // Based on GrevLex (x > y), LT(f) is x, LT(g) is x. 
    // Quotient 1, Remainder y.
    auto [q8, r8] = (x + y).div_mod(x);
    ASSERT_TEST(q8.to_scalar() == 1.0 && !r8.is_zero());

    // 9. Large Remainder
    Polynomial<double> p9 = x^5;
    auto [q9, r9] = p9.div_mod(x - 2.0);
    // x^5 / (x-2) -> Remainder is 2^5 = 32 (Remainder Theorem)
    ASSERT_TEST(r9.to_scalar() == 32.0);

    // 10. Non-trivial reduction: (x^2*y^2) / (x*y)
    auto [q10, r10] = ((x^2)*(y^2)).div_mod(x*y);
    ASSERT_TEST(r10.is_zero());
    ASSERT_TEST(q10.lead_monomial().degree() == 2);

    std::cout << "Passed.\n";
    global_results.record(true);
}

int main() {
    std::cout << "=============================================\n";
    std::cout << "PolyNum Algebraic Engine: Extended Test Suite\n";
    std::cout << "=============================================\n";
    
    test_polynomial_eval_substitute();
    test_polynomial_compose();
    test_polynomial_div_mod();

    std::cout << "=============================================\n";
    std::cout << "FINAL VERDICT:\n";
    std::cout << "Tests Modules Passed: " << global_results.passed << " / " << global_results.total << "\n";
    
    if (global_results.passed == global_results.total) {
        std::cout << "STATUS: SUCCESS - All components functioning correctly.\n";
    } else {
        std::cout << "STATUS: FAILURE - Check assertion logs above.\n";
        return 1;
    }
    std::cout << "=============================================\n";
    
    return 0;
}