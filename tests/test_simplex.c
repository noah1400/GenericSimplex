#include "../src/simplex.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures;

#define CHECK(cond) do { \
    if (!(cond)) { \
        printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        failures++; \
    } \
} while (0)

#define NEAR(a, b) (fabs((a) - (b)) < 1e-6)

static void test_max_le(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {3, 2};
    double a1[] = {1, 1}, a2[] = {1, 3};

    sx_problem_init(&p, 2, 2);
    sx_set_objective(&p, SX_MAXIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_LE, 4);
    sx_set_constraint(&p, 1, a2, SX_LE, 6);
    CHECK(sx_solve(&p, NULL, &s) == SX_OPTIMAL);
    CHECK(NEAR(s.objective, 12));
    CHECK(NEAR(s.x[0], 4) && NEAR(s.x[1], 0));
    sx_solution_free(&s);
    sx_problem_free(&p);
}

static void test_min_ge(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {2, 3};
    double a1[] = {1, 1}, a2[] = {1, 3};

    sx_problem_init(&p, 2, 2);
    sx_set_objective(&p, SX_MINIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_GE, 4);
    sx_set_constraint(&p, 1, a2, SX_GE, 6);
    CHECK(sx_solve(&p, NULL, &s) == SX_OPTIMAL);
    CHECK(NEAR(s.objective, 9));
    CHECK(NEAR(s.x[0], 3) && NEAR(s.x[1], 1));
    sx_solution_free(&s);
    sx_problem_free(&p);
}

static void test_mixed_eq(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {1, 2, 3};
    double a1[] = {1, 1, 1}, a2[] = {1, 0, 0}, a3[] = {0, 1, 0};

    sx_problem_init(&p, 3, 3);
    sx_set_objective(&p, SX_MAXIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_EQ, 10);
    sx_set_constraint(&p, 1, a2, SX_GE, 2);
    sx_set_constraint(&p, 2, a3, SX_LE, 3);
    CHECK(sx_solve(&p, NULL, &s) == SX_OPTIMAL);
    CHECK(NEAR(s.objective, 26));
    CHECK(NEAR(s.x[0], 2) && NEAR(s.x[1], 0) && NEAR(s.x[2], 8));
    sx_solution_free(&s);
    sx_problem_free(&p);
}

static void test_infeasible(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {1};
    double a1[] = {1}, a2[] = {1};

    sx_problem_init(&p, 1, 2);
    sx_set_objective(&p, SX_MAXIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_LE, 1);
    sx_set_constraint(&p, 1, a2, SX_GE, 2);
    CHECK(sx_solve(&p, NULL, &s) == SX_INFEASIBLE);
    CHECK(s.x == NULL);
    sx_problem_free(&p);
}

static void test_unbounded(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {1, 1};
    double a1[] = {1, -1};

    sx_problem_init(&p, 2, 1);
    sx_set_objective(&p, SX_MAXIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_LE, 1);
    CHECK(sx_solve(&p, NULL, &s) == SX_UNBOUNDED);
    sx_problem_free(&p);

    sx_problem_init(&p, 1, 0);
    sx_set_objective(&p, SX_MAXIMIZE, c);
    CHECK(sx_solve(&p, NULL, &s) == SX_UNBOUNDED);
    sx_problem_free(&p);
}

static void test_bounds(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {1, 1};
    double a1[] = {1, 1};

    sx_problem_init(&p, 2, 1);
    sx_set_var(&p, 0, "x", 1, 2);
    sx_set_var(&p, 1, "y", -1, 1);
    sx_set_objective(&p, SX_MAXIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_LE, 2.5);
    CHECK(sx_solve(&p, NULL, &s) == SX_OPTIMAL);
    CHECK(NEAR(s.objective, 2.5));
    CHECK(s.x[0] <= 2 + 1e-9 && s.x[0] >= 1 - 1e-9);
    CHECK(s.x[1] <= 1 + 1e-9 && s.x[1] >= -1 - 1e-9);
    sx_solution_free(&s);

    sx_set_var(&p, 0, "x", -INFINITY, 7);
    sx_set_var(&p, 1, "y", -3, -3);
    sx_set_constraint(&p, 0, a1, SX_LE, 100);
    CHECK(sx_solve(&p, NULL, &s) == SX_OPTIMAL);
    CHECK(NEAR(s.x[0], 7) && NEAR(s.x[1], -3));
    CHECK(NEAR(s.objective, 4));
    sx_solution_free(&s);
    sx_problem_free(&p);
}

static void test_shifted_rhs(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {1, 1};
    double a1[] = {-1, -1};

    sx_problem_init(&p, 2, 1);
    sx_set_var(&p, 0, "x", 2, INFINITY);
    sx_set_var(&p, 1, "y", 2, INFINITY);
    sx_set_objective(&p, SX_MINIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_LE, -5);
    CHECK(sx_solve(&p, NULL, &s) == SX_OPTIMAL);
    CHECK(NEAR(s.objective, 5));
    sx_solution_free(&s);
    sx_problem_free(&p);
}

static void test_free_var(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {1, 0};
    double a1[] = {1, 0}, a2[] = {1, 1};

    sx_problem_init(&p, 2, 2);
    sx_set_var(&p, 0, "x", -INFINITY, INFINITY);
    sx_set_objective(&p, SX_MINIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_GE, -3);
    sx_set_constraint(&p, 1, a2, SX_EQ, 0);
    CHECK(sx_solve(&p, NULL, &s) == SX_OPTIMAL);
    CHECK(NEAR(s.x[0], -3) && NEAR(s.x[1], 3));
    CHECK(NEAR(s.objective, -3));
    sx_solution_free(&s);
    sx_problem_free(&p);
}

static void test_redundant(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {1, 0};
    double a1[] = {1, 1}, a2[] = {2, 2};

    sx_problem_init(&p, 2, 2);
    sx_set_objective(&p, SX_MAXIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_EQ, 4);
    sx_set_constraint(&p, 1, a2, SX_EQ, 8);
    CHECK(sx_solve(&p, NULL, &s) == SX_OPTIMAL);
    CHECK(NEAR(s.x[0], 4) && NEAR(s.x[1], 0));
    sx_solution_free(&s);
    sx_problem_free(&p);
}

/* Beale's cycling example; terminates with Bland's rule. */
static void test_beale(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {0, 0, 0, -0.75, 20, -0.5, 6};
    double a1[] = {1, 0, 0, 0.25, -8, -1, 9};
    double a2[] = {0, 1, 0, 0.5, -12, -0.5, 3};
    double a3[] = {0, 0, 1, 0, 0, 1, 0};

    sx_problem_init(&p, 7, 3);
    sx_set_objective(&p, SX_MINIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_EQ, 0);
    sx_set_constraint(&p, 1, a2, SX_EQ, 0);
    sx_set_constraint(&p, 2, a3, SX_EQ, 1);
    CHECK(sx_solve(&p, NULL, &s) == SX_OPTIMAL);
    CHECK(NEAR(s.objective, -1.25));
    sx_solution_free(&s);
    sx_problem_free(&p);
}

static void test_invalid(void)
{
    sx_problem p;
    sx_solution s;
    double c[] = {1};

    sx_problem_init(&p, 1, 0);
    sx_set_var(&p, 0, "x", 2, 1);
    sx_set_objective(&p, SX_MINIMIZE, c);
    CHECK(sx_solve(&p, NULL, &s) == SX_ERROR);
    sx_problem_free(&p);
}

static void test_iter_limit(void)
{
    sx_problem p;
    sx_solution s;
    sx_options o = {0, 0, 1};
    double c[] = {2, 3};
    double a1[] = {1, 1}, a2[] = {1, 3};

    sx_problem_init(&p, 2, 2);
    sx_set_objective(&p, SX_MINIMIZE, c);
    sx_set_constraint(&p, 0, a1, SX_GE, 4);
    sx_set_constraint(&p, 1, a2, SX_GE, 6);
    CHECK(sx_solve(&p, &o, &s) == SX_ITER_LIMIT);
    sx_problem_free(&p);
}

int main(void)
{
    test_max_le();
    test_min_ge();
    test_mixed_eq();
    test_infeasible();
    test_unbounded();
    test_bounds();
    test_shifted_rhs();
    test_free_var();
    test_redundant();
    test_beale();
    test_invalid();
    test_iter_limit();
    if (failures) {
        printf("%d failure(s)\n", failures);
        return 1;
    }
    printf("all tests passed\n");
    return 0;
}
