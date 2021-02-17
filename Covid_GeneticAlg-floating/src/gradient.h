# ifndef __HEADER_GRADIENT__
# define __HEADER_GRADIENT__

# include <gsl/gsl_multimin.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_deriv.h>
# include "ga.h"
# include "sim.h"

# define GRADIENT_DIM 14
# define GRADIENT_MAX_ITERS 10

typedef struct {
    int i;
    gsl_vector * v;
    double * ic;
    double * params;
    fitness_func * ff;
} GradientParams;

int optimise_parameters(Genome * genome, fitness_func func);

# endif
