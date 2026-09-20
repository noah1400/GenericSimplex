#include "simplex.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define NONE ((size_t)-1)

typedef struct {
    double *t;      /* (m + 1) rows, w columns, last row is objective, last column rhs */
    size_t m, w;
    size_t *basis;
} tableau;

#define T(tb, i, j) ((tb)->t[(i) * (tb)->w + (j)])

/* original variable k = offset + sign * y[col] - y[col2] */
typedef struct {
    size_t col, col2;
    double offset, sign;
} vmap;

int sx_problem_init(sx_problem *p, size_t nvars, size_t ncons)
{
    size_t i;

    memset(p, 0, sizeof *p);
    p->nvars = nvars;
    p->ncons = ncons;
    p->sense = SX_MINIMIZE;
    p->vars = calloc(nvars ? nvars : 1, sizeof *p->vars);
    p->objective = calloc(nvars ? nvars : 1, sizeof *p->objective);
    p->cons = calloc(ncons ? ncons : 1, sizeof *p->cons);
    if (!p->vars || !p->objective || !p->cons)
        goto fail;
    for (i = 0; i < nvars; i++) {
        p->vars[i].lower = 0.0;
        p->vars[i].upper = INFINITY;
    }
    for (i = 0; i < ncons; i++) {
        p->cons[i].coef = calloc(nvars ? nvars : 1, sizeof *p->cons[i].coef);
        if (!p->cons[i].coef)
            goto fail;
    }
    return 0;
fail:
    sx_problem_free(p);
    return -1;
}

void sx_problem_free(sx_problem *p)
{
    size_t i;

    if (p->cons)
        for (i = 0; i < p->ncons; i++)
            free(p->cons[i].coef);
    free(p->cons);
    free(p->vars);
    free(p->objective);
    memset(p, 0, sizeof *p);
}

void sx_set_var(sx_problem *p, size_t i, const char *name, double lower, double upper)
{
    p->vars[i].name = name;
    p->vars[i].lower = lower;
    p->vars[i].upper = upper;
}

void sx_set_objective(sx_problem *p, sx_sense sense, const double *coef)
{
    p->sense = sense;
    memcpy(p->objective, coef, p->nvars * sizeof *coef);
}

void sx_set_constraint(sx_problem *p, size_t i, const double *coef, sx_rel rel, double rhs)
{
    memcpy(p->cons[i].coef, coef, p->nvars * sizeof *coef);
    p->cons[i].rel = rel;
    p->cons[i].rhs = rhs;
}

void sx_solution_free(sx_solution *s)
{
    free(s->x);
    s->x = NULL;
    s->nvars = 0;
}

const char *sx_status_str(sx_status st)
{
    switch (st) {
    case SX_OPTIMAL:    return "optimal";
    case SX_INFEASIBLE: return "infeasible";
    case SX_UNBOUNDED:  return "unbounded";
    case SX_ITER_LIMIT: return "iteration limit";
    case SX_ERROR:      return "error";
    }
    return "unknown";
}

static void pivot(tableau *tb, size_t r, size_t c)
{
    size_t i, j;
    double p = T(tb, r, c);

    for (j = 0; j < tb->w; j++)
        T(tb, r, j) /= p;
    for (i = 0; i <= tb->m; i++) {
        double f = T(tb, i, c);
        if (i == r || f == 0.0)
            continue;
        for (j = 0; j < tb->w; j++)
            T(tb, i, j) -= f * T(tb, r, j);
    }
    tb->basis[r] = c;
}

/* Minimises the objective row over the first ncols columns using Bland's rule. */
static sx_status run(tableau *tb, size_t ncols, const sx_options *o, size_t *iters)
{
    size_t rhs = tb->w - 1;

    for (;;) {
        size_t i, j, enter = ncols, leave = tb->m;
        double best = 0.0;

        for (j = 0; j < ncols; j++) {
            if (T(tb, tb->m, j) < -o->tol) {
                enter = j;
                break;
            }
        }
        if (enter == ncols)
            return SX_OPTIMAL;

        for (i = 0; i < tb->m; i++) {
            double a = T(tb, i, enter), ratio;
            if (a <= o->eps)
                continue;
            ratio = T(tb, i, rhs) / a;
            if (leave == tb->m || ratio < best - o->eps ||
                (fabs(ratio - best) <= o->eps && tb->basis[i] < tb->basis[leave])) {
                leave = i;
                best = ratio;
            }
        }
        if (leave == tb->m)
            return SX_UNBOUNDED;
        if (*iters >= o->max_iter)
            return SX_ITER_LIMIT;

        pivot(tb, leave, enter);
        (*iters)++;
    }
}

static void set_objective_row(tableau *tb, const double *cost)
{
    size_t i, j, rhs = tb->w - 1;

    for (j = 0; j < rhs; j++)
        T(tb, tb->m, j) = cost[j];
    T(tb, tb->m, rhs) = 0.0;
    for (i = 0; i < tb->m; i++) {
        double cb = cost[tb->basis[i]];
        if (cb == 0.0)
            continue;
        for (j = 0; j < tb->w; j++)
            T(tb, tb->m, j) -= cb * T(tb, i, j);
    }
}

/* Pivots artificial variables out of the basis, dropping redundant rows. */
static void remove_artificials(tableau *tb, size_t art_start, double eps)
{
    size_t i = 0, j;

    while (i < tb->m) {
        if (tb->basis[i] < art_start) {
            i++;
            continue;
        }
        for (j = 0; j < art_start; j++)
            if (fabs(T(tb, i, j)) > eps)
                break;
        if (j < art_start) {
            pivot(tb, i, j);
            i++;
            continue;
        }
        tb->m--;
        memmove(&T(tb, i, 0), &T(tb, tb->m, 0), tb->w * sizeof(double));
        memmove(&T(tb, tb->m, 0), &T(tb, tb->m + 1, 0), tb->w * sizeof(double));
        tb->basis[i] = tb->basis[tb->m];
    }
}

static int problem_valid(const sx_problem *p)
{
    size_t i, k;

    if (!p->vars || !p->objective || (p->ncons && !p->cons))
        return 0;
    for (k = 0; k < p->nvars; k++) {
        const sx_var *v = &p->vars[k];
        if (isnan(v->lower) || isnan(v->upper) || v->lower > v->upper)
            return 0;
        if (v->lower == INFINITY || v->upper == -INFINITY)
            return 0;
        if (!isfinite(p->objective[k]))
            return 0;
    }
    for (i = 0; i < p->ncons; i++) {
        if (!p->cons[i].coef || !isfinite(p->cons[i].rhs))
            return 0;
        for (k = 0; k < p->nvars; k++)
            if (!isfinite(p->cons[i].coef[k]))
                return 0;
    }
    return 1;
}

static size_t map_vars(const sx_problem *p, vmap *vm, size_t *nbounded)
{
    size_t k, col = 0;

    *nbounded = 0;
    for (k = 0; k < p->nvars; k++) {
        const sx_var *v = &p->vars[k];
        vm[k].col = col++;
        vm[k].col2 = NONE;
        vm[k].sign = 1.0;
        vm[k].offset = 0.0;
        if (v->lower == -INFINITY && v->upper == INFINITY) {
            vm[k].col2 = col++;
        } else if (v->lower == -INFINITY) {
            vm[k].offset = v->upper;
            vm[k].sign = -1.0;
        } else {
            vm[k].offset = v->lower;
            if (v->upper != INFINITY)
                (*nbounded)++;
        }
    }
    return col;
}

static void fill_row(double *row, double *rhs, const vmap *vm, size_t nvars, const double *coef)
{
    size_t k;

    for (k = 0; k < nvars; k++) {
        double a = coef[k];
        if (a == 0.0)
            continue;
        if (row) {
            row[vm[k].col] += a * vm[k].sign;
            if (vm[k].col2 != NONE)
                row[vm[k].col2] -= a;
        }
        *rhs -= a * vm[k].offset;
    }
}

static sx_rel flip(sx_rel rel)
{
    return rel == SX_LE ? SX_GE : rel == SX_GE ? SX_LE : SX_EQ;
}

sx_status sx_solve(const sx_problem *p, const sx_options *opts, sx_solution *s)
{
    sx_options o;
    tableau tb;
    vmap *vm = NULL;
    sx_rel *rels = NULL;
    double *cost = NULL, *y = NULL;
    size_t nc, ns, na, nb, m, i, k, rhs, slack, art;
    sx_status st = SX_ERROR;

    memset(s, 0, sizeof *s);
    memset(&tb, 0, sizeof tb);
    s->status = SX_ERROR;
    s->nvars = p->nvars;

    o.eps = opts && opts->eps > 0 ? opts->eps : 1e-9;
    o.tol = opts && opts->tol > 0 ? opts->tol : 1e-7;
    o.max_iter = opts ? opts->max_iter : 0;

    if (!problem_valid(p))
        return SX_ERROR;

    vm = malloc((p->nvars ? p->nvars : 1) * sizeof *vm);
    if (!vm)
        return SX_ERROR;
    nc = map_vars(p, vm, &nb);
    m = p->ncons + nb;

    rels = malloc((m ? m : 1) * sizeof *rels);
    if (!rels)
        goto done;
    ns = na = 0;
    for (i = 0; i < m; i++) {
        if (i < p->ncons) {
            double b = p->cons[i].rhs;
            fill_row(NULL, &b, vm, p->nvars, p->cons[i].coef);
            rels[i] = b < 0.0 ? flip(p->cons[i].rel) : p->cons[i].rel;
        } else {
            rels[i] = SX_LE;
        }
        if (rels[i] != SX_EQ)
            ns++;
        if (rels[i] != SX_LE)
            na++;
    }

    tb.m = m;
    tb.w = nc + ns + na + 1;
    tb.t = calloc((m + 1) * tb.w, sizeof *tb.t);
    tb.basis = malloc((m ? m : 1) * sizeof *tb.basis);
    cost = calloc(tb.w, sizeof *cost);
    y = calloc(nc ? nc : 1, sizeof *y);
    if (!tb.t || !tb.basis || !cost || !y)
        goto done;
    if (o.max_iter == 0)
        o.max_iter = 1000 * (m + tb.w) + 1000;

    rhs = tb.w - 1;
    slack = nc;
    art = nc + ns;

    for (i = 0; i < m; i++) {
        double *row = &T(&tb, i, 0);

        if (i < p->ncons) {
            row[rhs] = p->cons[i].rhs;
            fill_row(row, &row[rhs], vm, p->nvars, p->cons[i].coef);
        } else {
            size_t b = i - p->ncons;
            for (k = 0; k < p->nvars; k++) {
                if (p->vars[k].lower != -INFINITY && p->vars[k].upper != INFINITY) {
                    if (b == 0)
                        break;
                    b--;
                }
            }
            row[vm[k].col] = 1.0;
            row[rhs] = p->vars[k].upper - p->vars[k].lower;
        }

        if (row[rhs] < 0.0)
            for (k = 0; k < tb.w; k++)
                row[k] = -row[k];

        if (rels[i] == SX_LE) {
            row[slack] = 1.0;
            tb.basis[i] = slack++;
        } else {
            if (rels[i] == SX_GE)
                row[slack++] = -1.0;
            row[art] = 1.0;
            cost[art] = 1.0;
            tb.basis[i] = art++;
        }
    }

    if (na > 0) {
        set_objective_row(&tb, cost);
        st = run(&tb, tb.w - 1, &o, &s->iterations);
        if (st != SX_OPTIMAL)
            goto done;
        if (-T(&tb, tb.m, rhs) > o.tol) {
            st = SX_INFEASIBLE;
            goto done;
        }
        remove_artificials(&tb, nc + ns, o.eps);
    }

    memset(cost, 0, tb.w * sizeof *cost);
    for (k = 0; k < p->nvars; k++) {
        double c = p->sense == SX_MAXIMIZE ? -p->objective[k] : p->objective[k];
        cost[vm[k].col] += c * vm[k].sign;
        if (vm[k].col2 != NONE)
            cost[vm[k].col2] -= c;
    }
    set_objective_row(&tb, cost);
    st = run(&tb, nc + ns, &o, &s->iterations);
    if (st != SX_OPTIMAL)
        goto done;

    s->x = malloc((p->nvars ? p->nvars : 1) * sizeof *s->x);
    if (!s->x) {
        st = SX_ERROR;
        goto done;
    }
    for (i = 0; i < tb.m; i++)
        if (tb.basis[i] < nc)
            y[tb.basis[i]] = T(&tb, i, rhs);
    s->objective = 0.0;
    for (k = 0; k < p->nvars; k++) {
        double v = vm[k].offset + vm[k].sign * y[vm[k].col];
        if (vm[k].col2 != NONE)
            v -= y[vm[k].col2];
        s->x[k] = v;
        s->objective += p->objective[k] * v;
    }

done:
    free(tb.t);
    free(tb.basis);
    free(cost);
    free(y);
    free(vm);
    free(rels);
    s->status = st;
    return st;
}
