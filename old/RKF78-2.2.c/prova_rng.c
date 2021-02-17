#include <stdio.h>
#include <gsl/gsl_rng.h>

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  int i, n = 10;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  printf("The larger number is: %lu\n", gsl_rng_max(r));
  for (i = 0; i < n; i++)
    {
      unsigned long u = gsl_rng_get(r);
      printf ("%lu\n", u);
    }

  gsl_rng_free (r);

  return 0;
}
