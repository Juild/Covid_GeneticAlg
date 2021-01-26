# include "sim.h"

void phenotype_to_genotype(Genome * genome, IC * c1, Parameters * c2) {
	c1 -> E = crom2IC(genome -> c1[0]);
	c1 -> I1 = crom2IC(genome -> c1[1]);
	c1 -> A = crom2IC(genome -> c1[2]);

	c2 -> beta = crom2HSPar(genome -> c2[0]);
	c2 -> phi = crom2Par(genome -> c2[1]);
	c2 -> e1 = crom2LSPar(genome -> c2[2]);
	c2 -> eY = crom2LSPar(genome -> c2[3]);
	c2 -> sigma = crom2Par(genome -> c2[4]);
	c2 -> gamma1 = crom2Par(genome -> c2[5]);
	c2 -> gamma2 = crom2Par(genome -> c2[6]);
	c2 -> kappa = crom2LSPar(genome -> c2[7]);
	c2 -> p = crom2Par(genome -> c2[8]);
	c2 -> alpha = crom2HSPar(genome -> c2[9]);
	c2 -> delta = crom2HSPar(genome -> c2[10]);
}

void fitness_uniform(int day, double * rk_data, double * f) {
	for (int j = 0; j < N_PARAMS; j++) *f += gsl_pow_uint(rk_data[j] - DATA[day][j], 2);
}

void fitness_linear(int day, double * rk_data, double * f) {
	for (int j = 0; j < N_PARAMS; j++) *f += day * gsl_pow_uint(rk_data[j] - DATA[day][j], 2);
}

void fitness_exp(int day, double * rk_data, double * f) {
	for (int j = 0; j < N_PARAMS; j++) *f += exp(EXP_NU * day) * gsl_pow_uint(rk_data[j] - DATA[day][j], 2);
}

void fitness_max(int day, double * rk_data, double * f) {
	double ff = 0.0;
	for (int j = 0; j < N_PARAMS; j++) ff += gsl_pow_uint(rk_data[j] - DATA[day][j], 2);
	*f = GSL_MAX(ff, *f);
}

#define CoreModelDIM 8
void CoreModel(double t, double * x, unsigned dim, double * der, void * params){
	Parameters *par = (Parameters *) params; // To simplify the usage of Params (void pointer)
	double sigmae = par->sigma*x[1], gamma1i1 = par->gamma1*x[2], kappaA = par->kappa*x[3], alphai2 = par->alpha*x[5];
	der[0] = par->phi*x[2] + x[3] + (1-par->e1)*(x[4]+x[5]) + (1-par->eY)*x[6];
	der[0] = - par->beta * (x[0] * der[0])/POP_SIZE;
	der[1] = -der[0] - sigmae;
	der[2] = sigmae - gamma1i1;
	der[3] = (1-par->p)*gamma1i1 - kappaA - par->gamma2*x[3] ;
	der[4] = kappaA - par->gamma2*x[4];
	der[5] = par->p*gamma1i1 - par->gamma2*x[5] - alphai2;
	der[6] = alphai2 - (par->gamma2+par->delta)*x[6];
	der[7] = par->gamma2*(x[3] + x[4] + x[5] + x[6]);
}


int run_runge_putta(double * xt, void * ODE_pars, fitness_func func, double * fitness) {
	register int ndays;
	double t = 0.0, err, h = 1.e-3;

	double * rk_data;
	rk_data = (double *) malloc(N_PARAMS * sizeof(double));

	for (ndays = 1; ndays <= DAYS; ++ndays) {
		int status;
		while (t + h < ndays) {
			status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
			if (status) return status;
		}
		h = ndays - t;
		status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
		if (status) return status;

		rk_data[0] = xt[4];
		rk_data[1] = xt[5];
		rk_data[2] = xt[6];
		rk_data[3] = xt[7];
		rk_data[4] = POP_SIZE - (xt[0]+xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]+xt[7]);
		func(ndays, rk_data, fitness);
	}

	free(rk_data);

	return 0;
}


int compute_fitness(Genome * genome, fitness_func func) {
	IC * ic;
	ic = (IC *) malloc(sizeof(IC));
	Parameters * params;
	params = (Parameters *) malloc(sizeof(Parameters));
	phenotype_to_genotype(genome, ic, params);

	double xt[CoreModelDIM] = {POP_SIZE, ic->E, ic->I1, ic->A, DATA[0][0], DATA[0][1], DATA[0][2], DATA[0][3]};
	xt[0] -= (xt[1] + xt[2] + xt[3] + xt[4] + xt[5] + xt[6]);

	int status;
	double fitness = 0.0;
	if (xt[0] < 0 || (status = run_runge_putta(xt, params, func, &fitness))) {
		// printf("Oh shit, heere we go again! %d\n", status);
		genome->fitness = DBL_MAX;
		return 1;
	}

	genome->fitness = fitness;
	free(ic);
	free(params);
	return 0;
}
