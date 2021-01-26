# ifndef __HEADER_SIM__
# define __HEADER_SIM__

# include <float.h>
# include "RKF78.h"
# include "ga.h"
# include "utils.h"

# define POP_SIZE 1000000.000
# define CoreModelDIM 8

#define crom2IC(c) (((double) (c % 1000000000UL))/1000.0) // Initial conditions
#define crom2HSPar(c) (((double) (c % 1099511627776UL))/1099511627776.0) // High sensitivity
#define crom2Par(c) (((double) (c % 1048576U))/1048576.0) // Medium sensitivity
#define crom2LSPar(c) (((double) (c % 1024U))/1024.0) // Low sensitivity

typedef struct Parameters {

	double beta;
	double phi;
	double e1;
	double eY;
	double sigma;
	double gamma1;
	double gamma2;
	double kappa;
	double p;
	double alpha;
	double delta;

} Parameters;

#define HMAX 1.0
#define HMIN 1.e-3
#define RKTOL 1.e-5

void genotype_to_phenotype(Genome * genome, double * ic, Parameters * params);

/*
 * This function ideally should take care of calling genotype_to_phenotype,
 * calling evolve which runs the runge-kutta and computing the fitness using the
 * passed fit_func. This fitness should be added to the genome!!
 * returns the status, 0 if everything fine.
 */
int compute_fitness(Genome * genome, double (* fit_func) (double **));

/* Copy the code in the pdf!
 * @param x vector with initial conditions
 * @param params the set of Parameters that will be passed to the model
 * @param xt matrix with the evolution of the parameters, the dimension is DAYS x 5
 */
int evolve(double * x, void * params, double ** xt);

# define EXP_NU 0.05
double fitness_uniform(double ** solution);
double fitness_linear(double ** solution);
double fitness_exp(double ** solution);
double fitness_max(double ** solution);

# endif
