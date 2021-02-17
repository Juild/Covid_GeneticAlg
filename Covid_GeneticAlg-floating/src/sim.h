# ifndef __HEADER_SIM__
# define __HEADER_SIM__

# include <float.h>
# include "RKF78.h"
# include "ga.h"
# include "utils.h"

# define POP_SIZE 1000000.000
# define CoreModelDIM 8

#define crom2IC(c) (((double) (c))/1000.0) // Initial conditions

#define HMAX 1.0
#define HMIN 1.e-3
#define RKTOL 1.e-5

void genotype_to_phenotype(Genome * genome, double * ic);

/*
 * This is how a fitness function should be defined.
 * The first integer parameter is the day that corresponds to the array of data passed in the second parameter
 * whose size must be equal to `N_PARAMS`. The third aregument is a pointer to where the fitness should be stored.
 */
typedef void (* fitness_func)(int, double *, double *);

/*
 * This function ideally should take care of calling genotype_to_phenotype,
 * calling evolve which runs the runge-kutta and computing the fitness using the
 * passed fit_func. This fitness should be added to the genome!!
 * returns the status, 0 if everything fine.
 */
int compute_fitness(Genome * genome, fitness_func ff);

/* Copy the code in the pdf!
 * @param x vector with initial conditions
 * @param params the set of Parameters that will be passed to the model
 * @param ff the fitness function to use
 * @param fitness a pointer to where the fitness value calculated is stored
 */
int evolve(double * x, void * params, fitness_func ff, double * fitness);

int store_trajectory(double * xt, void * ODE_pars, FILE *outfile);

int run_runge_kutta(double * xt, void * ODE_pars, fitness_func func, double * fitness);

# define EXP_NU 0.05
# define NORM_UNIFORM 0.01
# define NORM_LINEAR 0.00019802			 // 1 over n(n + 1)/2
# define NORM_EXP 0.0003308427537690123  // evaluated for nu = 0.05, 1 over (1 - exp(nu*(N + 1)) / (1 - exp(nu)))
void fitness_uniform(int, double *, double *);
void fitness_linear(int, double *, double *);
void fitness_exp(int, double *, double *);
void fitness_max(int, double *, double *);

# endif
