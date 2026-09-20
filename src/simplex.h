#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <stddef.h>

typedef enum { SX_MAXIMIZE, SX_MINIMIZE } sx_sense;
typedef enum { SX_LE, SX_GE, SX_EQ } sx_rel;

typedef enum {
    SX_OPTIMAL,
    SX_INFEASIBLE,
    SX_UNBOUNDED,
    SX_ITER_LIMIT,
    SX_ERROR
} sx_status;

typedef struct {
    const char *name;
    double lower;   /* -INFINITY if unbounded below */
    double upper;   /*  INFINITY if unbounded above */
} sx_var;

typedef struct {
    double *coef;   /* nvars entries */
    sx_rel rel;
    double rhs;
} sx_constraint;

typedef struct {
    size_t nvars;
    size_t ncons;
    sx_sense sense;
    sx_var *vars;
    double *objective;      /* nvars entries */
    sx_constraint *cons;    /* ncons entries */
} sx_problem;

typedef struct {
    double eps;         /* pivot / zero tolerance, default 1e-9 */
    double tol;         /* optimality / feasibility tolerance, default 1e-7 */
    size_t max_iter;    /* 0 = automatic */
} sx_options;

typedef struct {
    sx_status status;
    double objective;
    double *x;          /* nvars entries, NULL unless SX_OPTIMAL */
    size_t nvars;
    size_t iterations;
} sx_solution;

int sx_problem_init(sx_problem *p, size_t nvars, size_t ncons);
void sx_problem_free(sx_problem *p);

void sx_set_var(sx_problem *p, size_t i, const char *name, double lower, double upper);
void sx_set_objective(sx_problem *p, sx_sense sense, const double *coef);
void sx_set_constraint(sx_problem *p, size_t i, const double *coef, sx_rel rel, double rhs);

sx_status sx_solve(const sx_problem *p, const sx_options *opts, sx_solution *s);
void sx_solution_free(sx_solution *s);

const char *sx_status_str(sx_status st);

#endif
