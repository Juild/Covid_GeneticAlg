#ifndef utils
#define utils
static const gsl_rng_type * T;
static gsl_rng * rng_generator;
gsl_rng_env_setup();
static T = gsl_rng_default;
static rng_generator = gsl_rng_alloc (T);

unsigned long random_num(){
  return gsl_rng_get(r);
}
double random_double(){
   return gsl_rng_uniform(r);
}
