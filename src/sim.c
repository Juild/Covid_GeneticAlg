# include "sim.h"

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
