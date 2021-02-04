# include "gradient.h"


void __gsl_vector_to_phenotype(
    const gsl_vector *v,
    double * c1,
    Parameters * c2
) {
    c1[1] = gsl_vector_get(v, 0);  // crom2IC(genome -> c1[0]);
    c1[2] = gsl_vector_get(v, 1);  // crom2IC(genome -> c1[1]);
    c1[3] = gsl_vector_get(v, 2);  // crom2IC(genome -> c1[2]);
    c1[4] = DATA[0][0];
    c1[5] = DATA[0][1];
    c1[6] = DATA[0][2];
    c1[7] = DATA[0][3];
    c1[0] = POP_SIZE - (c1[1] + c1[2] + c1[3] + c1[4] + c1[5] + c1[6]);

    c2 -> beta = gsl_vector_get(v, 3);  // crom2HSPar(genome -> c2[0]);
    c2 -> phi = gsl_vector_get(v, 4);  // crom2Par(genome -> c2[1]);
    c2 -> e1 = gsl_vector_get(v, 5);  // crom2LSPar(genome -> c2[2]);
    c2 -> eY = gsl_vector_get(v, 6);  // crom2LSPar(genome -> c2[3]);
    c2 -> sigma = gsl_vector_get(v, 7);  // crom2Par(genome -> c2[4]);
    c2 -> gamma1 = gsl_vector_get(v, 8);  // crom2Par(genome -> c2[5]);
    c2 -> gamma2 = gsl_vector_get(v, 9);  // crom2Par(genome -> c2[6]);
    c2 -> kappa = gsl_vector_get(v, 10);  // crom2LSPar(genome -> c2[7]);
    c2 -> p = gsl_vector_get(v, 11);  // crom2Par(genome -> c2[8]);
    c2 -> alpha = gsl_vector_get(v, 12);  // crom2HSPar(genome -> c2[9]);
    c2 -> delta = gsl_vector_get(v, 13);  // crom2HSPar(genome -> c2[10]);
}

void __genotype_to_gsl_vector(
    Genome * g,
    gsl_vector * v
) {
    gsl_vector_set(v, 0, crom2IC(g -> c1[0]));
    gsl_vector_set(v, 1, crom2IC(g -> c1[1]));
    gsl_vector_set(v, 2, crom2IC(g -> c1[2]));

    gsl_vector_set(v, 3, crom2HSPar(g -> c2[0]));
    gsl_vector_set(v, 4, crom2Par(g -> c2[1]));
    gsl_vector_set(v, 5, crom2LSPar(g -> c2[2]));
    gsl_vector_set(v, 6, crom2LSPar(g -> c2[3]));
    gsl_vector_set(v, 7, crom2Par(g -> c2[4]));
    gsl_vector_set(v, 8, crom2Par(g -> c2[5]));
    gsl_vector_set(v, 9, crom2Par(g -> c2[6]));
    gsl_vector_set(v, 10, crom2LSPar(g -> c2[7]));
    gsl_vector_set(v, 11, crom2Par(g -> c2[8]));
    gsl_vector_set(v, 12, crom2HSPar(g -> c2[9]));
    gsl_vector_set(v, 13, crom2HSPar(g -> c2[10]));
}

void __gsl_vector_to_genotype(
    Genome * g,
    gsl_vector * v
) {
    g->c1[0] = (unsigned long) (gsl_vector_get(v, 0) * 1000.0);
    g->c1[1] = (unsigned long) (gsl_vector_get(v, 1) * 1000.0);
    g->c1[2] = (unsigned long) (gsl_vector_get(v, 2) * 1000.0);

    g->c2[0] = (unsigned long) (gsl_vector_get(v, 3) * 1099511627776.0);
    g->c2[1] = (unsigned long) (gsl_vector_get(v, 4) * 1048576.0);
    g->c2[2] = (unsigned long) (gsl_vector_get(v, 5) * 1024.0);
    g->c2[3] = (unsigned long) (gsl_vector_get(v, 6) * 1024.0);
    g->c2[4] = (unsigned long) (gsl_vector_get(v, 7) * 1048576.0);
    g->c2[5] = (unsigned long) (gsl_vector_get(v, 8) * 1048576.0);
    g->c2[6] = (unsigned long) (gsl_vector_get(v, 9) * 1048576.0);
    g->c2[7] = (unsigned long) (gsl_vector_get(v, 10) * 1024.0);
    g->c2[8] = (unsigned long) (gsl_vector_get(v, 11) * 1048576.0);
    g->c2[9] = (unsigned long) (gsl_vector_get(v, 12) * 1099511627776.0);
    g->c2[10] = (unsigned long) (gsl_vector_get(v, 13) * 1099511627776.0);
}


double __fitness(
    const gsl_vector *v,
    void * p
) {
    fitness_func * func = (fitness_func *) p;

    double * ic;
	ic = (double *) malloc(CoreModelDIM * sizeof(double));
	Parameters * params;
	params = (Parameters *) malloc(sizeof(Parameters));
	__gsl_vector_to_phenotype(v, ic, params);

	double fitness = 0.0;
	if (run_runge_putta(ic, params, *func, &fitness))
		return DBL_MAX;

	free(ic);
	free(params);
	return fitness;
}

double __central_difference(
    double x,
    void * p
) {
    GradientParams * params = (GradientParams *) p;
    double fitness = 0.0;
    int param_to_change = params->i;

    double val = gsl_vector_get(params->v, param_to_change);
    gsl_vector_set(params->v, param_to_change, x);
    __gsl_vector_to_phenotype(params->v, params->ic, params->params);

    if (run_runge_putta(params->ic, params->params, *(params->ff), &fitness))
        return DBL_MAX;

    gsl_vector_set(params->v, param_to_change, val); // take it back so it can be reused
    return fitness;
}


/* The gradient of f, df = (df/dx, df/dy). */
void __fitness_gradient(
    const gsl_vector *v,
    void * p,
    gsl_vector *df
) {
  fitness_func * func = (fitness_func *) p;

  double * ic;
  ic = (double *) malloc(CoreModelDIM * sizeof(double));
  Parameters * params;
  params = (Parameters *) malloc(sizeof(Parameters));
  gsl_vector * w = gsl_vector_alloc(GRADIENT_DIM);
  gsl_vector_memcpy(w, v);
  GradientParams g_params = {0, w, ic, params, func};

  gsl_function F = {__central_difference, (void *) &g_params};

  double df_dxi, abserr;
  for (int i = 0; i < GRADIENT_DIM; i++) {
      g_params.i = i;
      gsl_deriv_central(&F, gsl_vector_get(v, i), 1e-8, &df_dxi, &abserr);
      gsl_vector_set(df, i, df_dxi);
  }

  gsl_vector_free(w);
  free(ic);
  free(params);
}

/* Compute both f and df together. */
void __fitness_f_gradient(
    const gsl_vector * v,
    void * p,
    double * f,
    gsl_vector * df
) {
    fitness_func * func = (fitness_func *) p;

    double * ic;
	ic = (double *) malloc(CoreModelDIM * sizeof(double));
	Parameters * params;
	params = (Parameters *) malloc(sizeof(Parameters));
	__gsl_vector_to_phenotype(v, ic, params);

	if (run_runge_putta(ic, params, *func, f))
		*f = DBL_MAX;

    gsl_vector * w = gsl_vector_alloc(GRADIENT_DIM);
    gsl_vector_memcpy(w, v);
    GradientParams g_params = {0, w, ic, params, func};

    gsl_function F = {__central_difference, (void *) &g_params};

    double df_dxi, abserr;
    for (int i = 0; i < GRADIENT_DIM; i++) {
        g_params.i = i;
        gsl_deriv_central(&F, gsl_vector_get(v, i), 1e-8, &df_dxi, &abserr);
        gsl_vector_set(df, i, df_dxi);
    }

    gsl_vector_free(w);
	free(ic);
	free(params);
}


/*
 * Returns 1 if the optimisation could decrease the fidelity.
 */
int optimise_parameters(Genome * genome, fitness_func func) {
    gsl_vector * v = gsl_vector_alloc(GRADIENT_DIM); // this vector will allocated the optimised parameters
    __genotype_to_gsl_vector(genome, v); // it is initialised with the parameters from the genome

    int status, iter = 0;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    gsl_multimin_function_fdf my_func;
    my_func.n = GRADIENT_DIM;
    my_func.f = __fitness;
    my_func.df = __fitness_gradient;
    my_func.fdf = __fitness_f_gradient;
    my_func.params = (void *) &func;

    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc(T, GRADIENT_DIM);
    gsl_multimin_fdfminimizer_set(s, &my_func, v, 0.1, 1e-4);

    do {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);
        if (status) break;
        status = gsl_multimin_test_gradient(s->gradient, 1e-3);
    } while (status == GSL_CONTINUE && iter < GRADIENT_MAX_ITERS);

    if (s->f < genome->fitness) {
        printf("Optimising the genome decreased the fidelity from %.2f to %.2f\n", genome->fitness, s->f);
        genome->fitness = s->f;
        __gsl_vector_to_genotype(genome, v);

        gsl_multimin_fdfminimizer_free(s);
        gsl_vector_free(v);
        return 1;
    }

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(v);

    return 0;
}
