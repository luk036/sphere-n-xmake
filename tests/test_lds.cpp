#include <doctest/doctest.h>  // for Approx, ResultBuilder, TestCase, CHECK

#include <gsl/span>     // for span
#include <sphere_n/sphere_n.hpp>  // for Circle, Halton, Sphere, Sphere3Hopf

TEST_CASE("Circle") {
    auto cgen = lds2::Circle(2);
    const auto res = cgen.pop();
    CHECK_EQ(res[0], doctest::Approx(0.0));
}

TEST_CASE("Halton") {
    const size_t base[] = {2, 3};
    auto hgen = lds2::Halton(base);
    const auto res = hgen.pop();
    CHECK_EQ(res[0], doctest::Approx(0.5));
}

TEST_CASE("Sphere") {
    const size_t base[] = {2, 3};
    auto sgen = lds2::Sphere(base);
    const auto res = sgen.pop();
    CHECK_EQ(res[0], doctest::Approx(0.8660254038));
}

TEST_CASE("Sphere3Hopf") {
    const size_t base[] = {2, 3, 5};
    auto shfgen = lds2::Sphere3Hopf(base);
    const auto res = shfgen.pop();
    CHECK_EQ(res[0], doctest::Approx(-0.2236067977));
}
