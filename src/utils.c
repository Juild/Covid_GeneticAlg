# include "utils.h"

void exit_error(const char *miss, int errcode)  {
    printf("\nERROR: %s.\nStopping...\n\n", miss); exit(errcode);
}

static gsl_rng * rng = NULL;

void init_rng() {
    rng = gsl_rng_alloc(gsl_rng_default);
    change_seed();
}

void change_seed() {
    unsigned seed = (unsigned) time(NULL);
    printf("Seed: %d\n", seed);
    gsl_rng_set(rng, seed);
}

void free_rng() {
    gsl_rng_free(rng);
}

double random_double() {
    return gsl_rng_uniform(rng);
}

unsigned long random_ulong() {
    return gsl_rng_get(rng);
}

unsigned int random_int(unsigned int n_max) {
    return gsl_rng_uniform_int(rng, n_max);
}

double random_gaussian(double sigma) {
    // return gsl_ran_gaussian(rng, sigma);
    return (random_double() * 2.0 - 1.0) * sigma;
}
