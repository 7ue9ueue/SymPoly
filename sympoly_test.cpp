#include <iostream>
#include <cassert>
#include <sstream>
#include <vector>
#include "sympoly.hpp"

// --- Verdict Tracking ---
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
// Expanded Test Functions
// --------------------------------------------------------------------------

void test_monomial_construction() {
    std::cout << "[Test] Construction (10 cases)... ";
    try {
        Monomial<int> m1; ASSERT_TEST(m1.degree() == 0);                         // Case 1: Default
        Monomial<int> m2(0, 5); ASSERT_TEST(m2.degree() == 5);                   // Case 2: x0^5
        Monomial<int> m3(100, 1); ASSERT_TEST(m3.degree() == 1);                 // Case 3: High index
        Monomial<int> m4(5, 0); ASSERT_TEST(m4.degree() == 0);                   // Case 4: Zero exp
        
        Monomial<int> m5(1, 1);
        Monomial<int> m6 = m5; ASSERT_TEST(m6.degree() == 1);                    // Case 5: Copy
        
        // Exceptional cases
        bool caught = false;
        try { Monomial<int> m_neg(1, -1); } catch(...) { caught = true; }
        ASSERT_TEST(caught);                                                     // Case 6: Negative exp
        
        // Large values
        Monomial<int> m7(0, 1000); ASSERT_TEST(m7.degree() == 1000);             // Case 7: Large degree
        
        Monomial<int> m8;
        Monomial<int> m9 = std::move(m1); ASSERT_TEST(m9.degree() == 0);         // Case 8: Move
        
        Monomial<int> m10(2, 1); ASSERT_TEST(!(m10 == m1));                      // Case 9: Inequality
        Monomial<int> m11(0, 5); ASSERT_TEST(m11 == m2);                         // Case 10: Equality
        
        std::cout << "Passed.\n";
        global_results.record(true);
    } catch (...) { global_results.record(false); }
}

void test_monomial_multiplication() {
    std::cout << "[Test] Multiplication (10 cases)... ";
    Monomial<int> x0(0, 1), x1(1, 1), x2(2, 1), one;

    ASSERT_TEST((x0 * one) == x0);                             // 1. Identity
    ASSERT_TEST((x0 * x1) == (x1 * x0));                       // 2. Commutativity
    ASSERT_TEST((x0 * x0).degree() == 2);                      // 3. Same variable
    ASSERT_TEST(((x0 * x1) * x2) == (x0 * (x1 * x2)));         // 4. Associativity
    ASSERT_TEST((x0 * Monomial<int>(0, 0)) == x0);             // 5. Mult by deg 0
    ASSERT_TEST((Monomial<int>(10, 2) * Monomial<int>(10, 3)).degree() == 5); // 6. High index add
    ASSERT_TEST((x0 * x1 * x2).degree() == 3);                 // 7. Chain multiplication
    
    Monomial<int> large(0, 50);
    ASSERT_TEST((large * large).degree() == 100);              // 8. Large exponents
    
    Monomial<int> mixed1(0, 2); mixed1 = mixed1 * Monomial<int>(5, 3);
    Monomial<int> mixed2(5, 1); mixed2 = mixed2 * Monomial<int>(0, 1);
    ASSERT_TEST((mixed1 * mixed2).degree() == 7);              // 9. Mixed indices
    ASSERT_TEST((one * one) == one);                           // 10. 1 * 1 = 1

    std::cout << "Passed.\n";
    global_results.record(true);
}

void test_monomial_division() {
    std::cout << "[Test] Division (10 cases)... ";
    Monomial<int> x0_5(0, 5), x0_2(0, 2), x1_2(1, 2), one;

    ASSERT_TEST((x0_5 / x0_2).degree() == 3);                  // 1. Simple division
    ASSERT_TEST((x0_5 / one) == x0_5);                         // 2. Divide by 1
    ASSERT_TEST((x0_5 / x0_5) == one);                         // 3. Self division
    
    auto mixed = x0_5 * x1_2;
    ASSERT_TEST((mixed / x1_2) == x0_5);                       // 4. Multi-var division
    ASSERT_TEST((mixed / x0_5) == x1_2);                       // 5. Multi-var division reverse
    
    // Failures
    bool c1 = false; try { x0_2 / x0_5; } catch(...) { c1 = true; }
    ASSERT_TEST(c1);                                           // 6. Underflow exponent
    
    bool c2 = false; try { x0_5 / x1_2; } catch(...) { c2 = true; }
    ASSERT_TEST(c2);                                           // 7. Missing variable in numerator
    
    Monomial<int> high(100, 10);
    ASSERT_TEST((high / Monomial<int>(100, 4)).degree() == 6); // 8. High index division
    
    ASSERT_TEST((one / one) == one);                           // 9. 1/1
    
    Monomial<int> x0x1x2 = (Monomial<int>(0,1)*Monomial<int>(1,1)*Monomial<int>(2,1));
    ASSERT_TEST((x0x1x2 / Monomial<int>(1,1)).degree() == 2);  // 10. Middle var removal

    std::cout << "Passed.\n";
    global_results.record(true);
}

void test_divisibility_check() {
    std::cout << "[Test] Divisibility (10 cases)... ";
    Monomial<int> x0_2(0, 2), x1_2(1, 2), one;
    auto x0x1 = x0_2 * x1_2;

    ASSERT_TEST(x0x1.is_divisible_by(x0_2));                   // 1. Basic factor
    ASSERT_TEST(x0x1.is_divisible_by(one));                    // 2. All divisible by 1
    ASSERT_TEST(one.is_divisible_by(one));                     // 3. 1 div by 1
    ASSERT_TEST(!one.is_divisible_by(x0_2));                   // 4. 1 not div by x
    ASSERT_TEST(!x0_2.is_divisible_by(x1_2));                  // 5. Different vars
    ASSERT_TEST(x0_2.is_divisible_by(Monomial<int>(0, 1)));    // 6. Higher power
    ASSERT_TEST(!Monomial<int>(0, 1).is_divisible_by(x0_2));   // 7. Lower power
    ASSERT_TEST(x0x1.is_divisible_by(x0x1));                   // 8. Self
    
    Monomial<int> high(50, 1);
    ASSERT_TEST(high.is_divisible_by(one));                    // 9. High index vs 1
    ASSERT_TEST(!one.is_divisible_by(high));                   // 10. 1 vs high index

    std::cout << "Passed.\n";
    global_results.record(true);
}

void test_grevlex_ordering() {
    std::cout << "[Test] GrevLex (10 cases)... ";
    Monomial<int> x(0, 1), y(1, 1), z(2, 1), one;

    // Standard properties
    ASSERT_TEST(one < x);                                      // 1. Constant is smallest
    ASSERT_TEST(x < (x * x));                                  // 2. Degree weight
    ASSERT_TEST(y < x);                                        // 3. GrevLex: x > y
    ASSERT_TEST(z < y);                                        // 4. GrevLex: y > z
    
    // Degree ties
    // x^2 vs x*y -> x^2 is larger
    ASSERT_TEST((x * y) < (x * x));                            // 5. Degree tie
    
    // The GrevLex special: x*z vs y^2
    // Both degree 2. Last var is z. x*z has z^1, y^2 has z^0.
    // Larger trailing exponent makes it SMALLER.
    ASSERT_TEST((x * z) < (y * y));                            // 6. GrevLex specific tie-break
    
    // x*y^2*z (1,2,1) vs y^3*z (0,3,1)
    auto A = x * (y*y) * z;
    auto B = (y*y*y) * z;
    ASSERT_TEST(B < A);                                        // 7. Classic Buchberger test
    
    ASSERT_TEST(!(x < x));                                     // 8. Strictness
    ASSERT_TEST(one < Monomial<int>(100, 1));                  // 9. High index
    ASSERT_TEST((x) < (z*z*z));                                // 10. Higher degree always wins

    std::cout << "Passed.\n";
    global_results.record(true);
}

void test_evaluation() {
    std::cout << "[Test] Evaluation (10 cases)... ";
    Monomial<double> m(0, 2); m = m * Monomial<double>(1, 2); // x0^2 * x1^2

    auto res1 = m.eval_partial({{0, 2.0}});
    ASSERT_TEST(res1.first == 4.0 && res1.second.degree() == 2); // 1. Partial x0
    
    auto res2 = m.eval_partial({{1, 3.0}});
    ASSERT_TEST(res2.first == 9.0 && res2.second.degree() == 2); // 2. Partial x1
    
    auto res3 = m.eval_partial({{0, 2.0}, {1, 2.0}});
    ASSERT_TEST(res3.first == 16.0 && res3.second.degree() == 0);// 3. Full evaluation
    
    Monomial<double> one;
    ASSERT_TEST(one.eval_partial({{0, 5.0}}).first == 1.0);      // 4. Eval constant
    
    auto res5 = m.eval_partial({{5, 10.0}});
    ASSERT_TEST(res5.first == 1.0 && res5.second == m);          // 5. Eval non-existent var
    
    auto res6 = m.eval_partial({{0, 0.0}});
    ASSERT_TEST(res6.first == 0.0);                              // 6. Zero multiplication
    
    Monomial<double> high(10, 3);
    ASSERT_TEST(high.eval_partial({{10, 2.0}}).first == 8.0);    // 7. High index eval
    
    Monomial<double> linear(0, 1);
    ASSERT_TEST(linear.eval_partial({{0, -1.0}}).first == -1.0); // 8. Negative field values
    
    auto res9 = m.eval_partial({{0, 1.0}, {1, 1.0}});
    ASSERT_TEST(res9.first == 1.0);                              // 9. Identity eval
    
    auto res10 = (Monomial<double>(0,1)*Monomial<double>(0,1)).eval_partial({{0, 3.0}});
    ASSERT_TEST(res10.first == 9.0);                             // 10. Combined power eval

    std::cout << "Passed.\n";
    global_results.record(true);
}
void test_formatting() {
    std::cout << "[Test] Formatting (10 cases)...\n";

    // Setup: Initialize the registry with specific names
    // ID 0 = x, ID 1 = y, ID 2 = z, ... ID 10 = alpha
    Symbol x("x"), y("y"), z("z"), w("w"), p("p"), q("q"), 
           r("r"), s("s"), t("t"), u("u"), alpha("alpha");

    auto check = [](Monomial<int> m, std::string expected, int caseNum) {
        std::stringstream ss; 
        ss << m;
        std::string result = ss.str();
        
        // Print to cout for visual verification as requested
        std::cout << "  Case " << caseNum << ": Expected \"" << expected 
                  << "\", Got \"" << result << "\"";
        
        if (result == expected) {
            std::cout << " [PASS]\n";
            return true;
        } else {
            std::cout << " [FAIL]\n";
            return false;
        }
    };

    // 1. Identity / Constant
    ASSERT_TEST(check(Monomial<int>(), "1", 1));

    // 2. Single variable (using ID 0 -> "x")
    ASSERT_TEST(check(Monomial<int>(0, 1), "x", 2));

    // 3. Single variable with power
    ASSERT_TEST(check(Monomial<int>(0, 2), "x^2", 3));

    // 4. Multiple variables (x * y)
    ASSERT_TEST(check(Monomial<int>(0, 1) * Monomial<int>(1, 1), "x*y", 4));

    // 5. Gap in indices (ID 0 and ID 2 -> "x*z")
    ASSERT_TEST(check(Monomial<int>(0, 1) * Monomial<int>(2, 1), "x*z", 5));

    // 6. Natural Sorting (IDs are 0, 1; should print x then y regardless of mult order)
    ASSERT_TEST(check(Monomial<int>(1, 2) * Monomial<int>(0, 1), "x*y^2", 6));

    // 7. High power formatting
    ASSERT_TEST(check(Monomial<int>(0, 10), "x^10", 7));

    // 8. Named multi-character symbols (ID 10 -> "alpha")
    ASSERT_TEST(check(Monomial<int>(10, 1), "alpha", 8));

    // 9. Complex mixed expression (x^2 * z^3)
    Monomial<int> complex = Monomial<int>(0, 2) * Monomial<int>(2, 3);
    ASSERT_TEST(check(complex, "x^2*z^3", 9));

    // 10. Large gap with multiple powers (y^2 * alpha^5)
    Monomial<int> large_gap = Monomial<int>(1, 2) * Monomial<int>(10, 5);
    ASSERT_TEST(check(large_gap, "y^2*alpha^5", 10));

    std::cout << "Formatting Tests Passed.\n";
    global_results.record(true);
}

// --------------------------------------------------------------------------
// Main Execution
// --------------------------------------------------------------------------

int main() {
    std::cout << "=============================================\n";
    std::cout << "PolyNum Algebraic Engine: Monomial Test Suite\n";
    std::cout << "=============================================\n";
    
    test_monomial_construction();
    test_monomial_multiplication();
    test_monomial_division();
    test_divisibility_check();
    test_grevlex_ordering();
    test_evaluation();
    test_formatting();

    std::cout << "=============================================\n";
    std::cout << "FINAL VERDICT:\n";
    std::cout << "Tests Modules Passed: " << global_results.passed << " / " << global_results.total << "\n";
    
    if (global_results.passed == global_results.total) {
        std::cout << "STATUS: SUCCESS - Monomial class is mathematically sound.\n";
    } else {
        std::cout << "STATUS: FAILURE - Check assertion logs above.\n";
        return 1;
    }
    std::cout << "=============================================\n";
    
    return 0;
}