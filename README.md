# GenericSimplex

Two-phase simplex solver for linear programs in C99. No dependencies beyond libm.

## Problem form

```
min / max   c^T x
subject to  a_i^T x  (<=, >=, =)  b_i     for each constraint i
            l_k <= x_k <= u_k              for each variable k
```

Bounds may be `-INFINITY` / `INFINITY`; free variables are supported.

## Usage

```c
#include <math.h>
#include "simplex.h"

sx_problem p;
sx_solution s;
double obj[] = {3, 5};
double r0[] = {1, 0}, r1[] = {0, 2}, r2[] = {3, 2};

sx_problem_init(&p, 2, 3);                 /* 2 variables, 3 constraints */
sx_set_var(&p, 0, "x", 0, INFINITY);
sx_set_var(&p, 1, "y", 0, INFINITY);
sx_set_objective(&p, SX_MAXIMIZE, obj);
sx_set_constraint(&p, 0, r0, SX_LE, 4);
sx_set_constraint(&p, 1, r1, SX_LE, 12);
sx_set_constraint(&p, 2, r2, SX_LE, 18);

if (sx_solve(&p, NULL, &s) == SX_OPTIMAL)
    printf("%g at x=%g y=%g\n", s.objective, s.x[0], s.x[1]);

sx_solution_free(&s);
sx_problem_free(&p);
```

`sx_solve` returns `SX_OPTIMAL`, `SX_INFEASIBLE`, `SX_UNBOUNDED`, `SX_ITER_LIMIT`
or `SX_ERROR` (invalid input or allocation failure). `s.x` is only set on
`SX_OPTIMAL`. Pass an `sx_options` to override tolerances or the iteration limit.

Pivoting uses Bland's rule, so the solver cannot cycle on degenerate problems.

## Build

```
make          # libsimplex.a and example
make test
```
