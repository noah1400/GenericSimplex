#include "src/simplex.h"

#include <math.h>
#include <stdio.h>

/*
 * maximize   3x + 5y
 * subject to  x       <= 4
 *                 2y  <= 12
 *             3x + 2y <= 18
 *             x, y >= 0
 */
int main(void)
{
    sx_problem p;
    sx_solution s;
    size_t i;
    double obj[] = {3, 5};
    double r0[] = {1, 0};
    double r1[] = {0, 2};
    double r2[] = {3, 2};

    if (sx_problem_init(&p, 2, 3) != 0)
        return 1;
    sx_set_var(&p, 0, "x", 0, INFINITY);
    sx_set_var(&p, 1, "y", 0, INFINITY);
    sx_set_objective(&p, SX_MAXIMIZE, obj);
    sx_set_constraint(&p, 0, r0, SX_LE, 4);
    sx_set_constraint(&p, 1, r1, SX_LE, 12);
    sx_set_constraint(&p, 2, r2, SX_LE, 18);

    sx_solve(&p, NULL, &s);
    printf("status: %s\n", sx_status_str(s.status));
    if (s.status == SX_OPTIMAL) {
        printf("objective: %g\n", s.objective);
        for (i = 0; i < s.nvars; i++)
            printf("  %s = %g\n", p.vars[i].name, s.x[i]);
    }
    printf("iterations: %zu\n", s.iterations);

    sx_solution_free(&s);
    sx_problem_free(&p);
    return 0;
}
