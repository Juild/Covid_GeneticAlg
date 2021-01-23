# include "sim.h"
#define HMAX 1.0
#define HMIN 1.e-3
#define RKTOL 1.e-5

void phenotype_to_genotype(Genome * genome, IC * c1, Parameters * c2) {
	c1 -> E = crom2IC(genome -> c1[0]);
	c1 -> I1 = crom2IC(genome -> c1[1]);
	c1 -> A = crom2IC(genome -> c1[2]);

	c2 -> beta = crom2HSPar(genome -> c1[0]);
	c2 -> phi = crom2Par(genome -> c1[1]);
	c2 -> e1 = crom2LSPar(genome -> c1[2]);
	c2 -> eY = crom2LSPar(genome -> c1[3]);
	c2 -> sigma = crom2Par(genome -> c1[4]);
	c2 -> gamma1 = crom2Par(genome -> c1[5]);
	c2 -> gamma2 = crom2Par(genome -> c1[6]);
	c2 -> kappa = crom2LSPar(genome -> c1[7]);
	c2 -> p = crom2Par(genome -> c1[8]);
	c2 -> alpha = crom2HSPar(genome -> c1[9]);
	c2 -> delta = crom2HSPar(genome -> c1[10]);
}

double fitness_uniform(double ** solution) {
	int i, j;
	double f = 0.0;
	for (i = 0; i < DAYS; i++)
		for (j = 0; i < 5; i++)
			f += gsl_pow_uint(solution[i][j] - DATA[i][j], 2);

	return f;
}

double fitness_linear(double ** solution) {
	int i, j;
	double f = 0.0;
	for (i = 0; i < DAYS; i++)
		for (j = 0; i < 5; i++)
			f += i * gsl_pow_uint(solution[i][j] - DATA[i][j], 2);

	return f;
}

double fitness_exp(double ** solution) {
	int i, j;
	double f = 0.0;
	for (i = 0; i < DAYS; i++)
		for (j = 0; i < 5; i++)
			f += exp(EXP_NU * i) * gsl_pow_uint(solution[i][j] - DATA[i][j], 2);

	return f;
}

double fitness_max(double ** solution) {
	int i, j;
	double f, max_f = -1;
	for (i = 0; i < DAYS; i++) {
		f = 0.0;
		for (j = 0; i < 5; i++) f += gsl_pow_uint(solution[i][j] - DATA[i][j], 2);
		max_f = GSL_MAX(max_f, f);
	}

	return max_f;
}

#define CoreModelDIM 8
void CoreModel(double t, double * x, unsigned dim, double * der, void * params){
	Parameters *par = (Parameters *) params; // To simplify the usage of Params (void pointer)
	double sigmae = par->sigma*x[1], gamma1i1 = par->gamma1*x[2], kappaA = par->kappa*x[3], alphai2 = par->alpha*x[5];
	der[0] = par->phi*x[2] + x[3] + (1-par->e1)*(x[4]+x[5]) + (1-par->eY)*x[6];
	der[0] = - par->beta * (x[0] * der[0])/par->PopSize;
	der[1] = -der[0] - sigmae;
	der[2] = sigmae - gamma1i1;
	der[3] = (1-par->p)*gamma1i1 - kappaA - par->gamma2*x[3] ;
	der[4] = kappaA - par->gamma2*x[4];
	der[5] = par->p*gamma1i1 - par->gamma2*x[5] - alphai2;
	der[6] = alphai2 - (par->gamma2+par->delta)*x[6];
	der[7] = par->gamma2*(x[3] + x[4] + x[5] + x[6]);
}



int run_runge_putta(double *xt, void * ODE_pars, double ** out_data) {
	register unsigned ndays;
	double t = 0.0, err, h = 1.e-3;
	for (ndays = 1; ndays <= DAYS; ++ndays) {
		int status;
		while (t + h < ndays) {
			status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
			if (status) return status;
		}
		h = ndays - t;
		status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
		if (status) return status;

		out_data[ndays][0] = xt[4];
		out_data[ndays][1] = xt[5];
		out_data[ndays][2] = xt[6];
		out_data[ndays][3] = xt[7];
		out_data[ndays][4] = POP_SIZE - (xt[0]+xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]+xt[7]);
	}

	return 0;
}


int compute_fitness(Genome * genome, double (* fit_func) (double **)) {
	// fit_func(data);
	IC * ic;
	Parameters * params;
	phenotype_to_genotype(genome, ic, params);

	double xt[CoreModelDIM] = {POP_SIZE, ic->E, ic->I1, ic->A, DATA[0][0],DATA[0][1],DATA[0][2],DATA[0][3]};
	xt[0] -= (xt[1] + xt[2] + xt[3] + xt[4] + xt[5] + xt[6]);
	double rk_data[DAYS][5];

	CoreModel(double t, double * x, unsigned dim, double * der, void * params)

	run_runge_putta(double *xt, void * ODE_pars, double ** out_data)

	if (run_runge_kutta(xt, params, rk_data)) {
		printf("Oh shit, heere we go again!");
		genome->fitness = 999999999999999999999999999999999.9;
		return 1;
	}

	genome->fitness = fit_func(rk_data);
	return 0;
}
