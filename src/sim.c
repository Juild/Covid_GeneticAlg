# include "sim.h"

void genotype_to_phenotype(Genome * genome, double * c1) {
	c1[1] = crom2IC(genome -> c1[0]);
	c1[2] = crom2IC(genome -> c1[1]);
	c1[3] = crom2IC(genome -> c1[2]);
	c1[4] = DATA[0][0];
	c1[5] = DATA[0][1];
	c1[6] = DATA[0][2];
	c1[7] = DATA[0][3];
	c1[0] = POP_SIZE - (c1[1] + c1[2] + c1[3] + c1[4] + c1[5] + c1[6]);
}

void fitness_uniform(int day, double * rk_data, double * f) {
	*f += NORM_UNIFORM * (
		gsl_pow_2(rk_data[0] - DATA[day][0]) +
		gsl_pow_2(rk_data[1] - DATA[day][1]) +
		gsl_pow_2(rk_data[2] - DATA[day][2]) +
		gsl_pow_2(rk_data[3] - DATA[day][3]) +
		gsl_pow_2(rk_data[4] - DATA[day][4])
	);
}

void fitness_linear(int day, double * rk_data, double * f) {
	*f += NORM_LINEAR * day * (
		gsl_pow_2(rk_data[0] - DATA[day][0]) +
		gsl_pow_2(rk_data[1] - DATA[day][1]) +
		gsl_pow_2(rk_data[2] - DATA[day][2]) +
		gsl_pow_2(rk_data[3] - DATA[day][3]) +
		gsl_pow_2(rk_data[4] - DATA[day][4])
	);
}

void fitness_exp(int day, double * rk_data, double * f) {
	*f += NORM_EXP * exp(EXP_NU * day) * (
		gsl_pow_2(rk_data[0] - DATA[day][0]) +
		gsl_pow_2(rk_data[1] - DATA[day][1]) +
		gsl_pow_2(rk_data[2] - DATA[day][2]) +
		gsl_pow_2(rk_data[3] - DATA[day][3]) +
		gsl_pow_2(rk_data[4] - DATA[day][4])
	);
}

void fitness_max(int day, double * rk_data, double * f) {
	double ff = gsl_pow_2(rk_data[0] - DATA[day][0]) +
				gsl_pow_2(rk_data[1] - DATA[day][1]) +
				gsl_pow_2(rk_data[2] - DATA[day][2]) +
				gsl_pow_2(rk_data[3] - DATA[day][3]) +
				gsl_pow_2(rk_data[4] - DATA[day][4]);
	*f = GSL_MAX(ff, *f);
}

void CoreModel(double t, double * x, unsigned dim, double * der, void * params) {
	double *par = (double *) params; // To simplify the usage of Params (void pointer)
	double sigmae = par[4]*x[1], gamma1i1 = par[5]*x[2], kappaA = par[7]*x[3], alphai2 = par[9]*x[5];
	der[0] = par[1]*x[2] + x[3] + (1-par[2])*(x[4]+x[5]) + (1-par[3])*x[6];
	der[0] = - par[0] * (x[0] * der[0]) / (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]);
	der[1] = -der[0] - sigmae;
	der[2] = sigmae - gamma1i1;
	der[3] = (1-par[8])*gamma1i1 - kappaA - par[6]*x[3] ;
	der[4] = kappaA - par[6]*x[4];
	der[5] = par[8]*gamma1i1 - par[6]*x[5] - alphai2;
	der[6] = alphai2 - (par[6]+par[10])*x[6];
	der[7] = par[6]*(x[3] + x[4] + x[5] + x[6]);
}


int run_runge_putta(double * xt, void * ODE_pars, fitness_func func, double * fitness) {
	register int ndays;
	double t = 0.0, err, h = 1.e-3;

	double * rk_data;
	rk_data = (double *) malloc(N_PARAMS * sizeof(double));
	*fitness = 0.0;
	for (ndays = 1; ndays < DAYS; ++ndays) {
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

int store_trajectory(double * xt, void * ODE_pars, FILE *outfile) {
	register int ndays;
	double t = 0.0, err, h = 1.e-3;
	fprintf(outfile, "\nD,S,E,I_1,A,A_d,I_2,Y,R\n");

	double * rk_data;
	rk_data = (double *) malloc(N_PARAMS * sizeof(double));

	double fu = 0.0, fl = 0.0, fe = 0.0, fm = 0.0;

	for (ndays = 1; ndays < DAYS; ++ndays) {
		int status;
		while (t + h < ndays) {
			status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
			if (status) return status;
		}
		h = ndays - t;
		status = RKF78Sys(&t, xt, CoreModelDIM, &h, &err, HMIN, HMAX, RKTOL, ODE_pars, CoreModel);
		if (status) return status;

		fprintf(outfile, "%.16f,", POP_SIZE - (xt[0]+xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]+xt[7]));
		fprintf(outfile, "%.16f,", xt[0]);
		fprintf(outfile, "%.16f,", xt[1]);
		fprintf(outfile, "%.16f,", xt[2]);
		fprintf(outfile, "%.16f,", xt[3]);
		fprintf(outfile, "%.16f,", xt[4]); // A_d
		fprintf(outfile, "%.16f,", xt[5]);
		fprintf(outfile, "%.16f,", xt[6]);
		fprintf(outfile, "%.16f\n", xt[7]);

		rk_data[0] = xt[4];
		rk_data[1] = xt[5];
		rk_data[2] = xt[6];
		rk_data[3] = xt[7];
		rk_data[4] = POP_SIZE - (xt[0]+xt[1]+xt[2]+xt[3]+xt[4]+xt[5]+xt[6]+xt[7]);

		fitness_exp(ndays, rk_data, &fe);
		fitness_linear(ndays, rk_data, &fl);
		fitness_uniform(ndays, rk_data, &fu);
		fitness_max(ndays, rk_data, &fm);
	}

	fprintf(outfile, "F_linear: %.5f\n", fl);
	fprintf(outfile, "F_uniform: %.5f\n", fu);
	fprintf(outfile, "F_exp: %.5f\n", fe);
	fprintf(outfile, "F_max: %.5f\n", fm);

	printf("F_linear: %.5f\n", fl);
	printf("F_uniform: %.5f\n", fu);
	printf("F_exp: %.5f\n", fe);
	printf("F_max: %.5f\n", fm);

	free(rk_data);
	return 0;
}


int compute_fitness(Genome * genome, fitness_func func) {
	// fit_func(data);
	double * ic;
	ic = (double *) malloc(CoreModelDIM * sizeof(double));
	genotype_to_phenotype(genome, ic);

	int status;
	double fitness = 0.0;
	if (ic[0] < 0 || (status = run_runge_putta(ic, genome->c2, func, &fitness))) {
		// if the first ic is negative we can skip the calculation, we know that for sure it is unfeasible
		// printf("Oh shit, heere we go again! %d\n", status);
		genome->fitness = DBL_MAX;
		return 1;
	}
	genome->fitness = fitness;
	free(ic);
	return 0;
}
