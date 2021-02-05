# ifndef __HEADER__GA__
# define __HEADER__GA__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include "utils.h"

# define UL_SIZE sizeof(unsigned long)

# define GENES_C1 3
# define GENES_C2 11

// we are going to use a two chromosome genome,
typedef struct Genome {
	unsigned long c1[GENES_C1];
	unsigned long c2[GENES_C2];
	double fitness;
} Genome;

void generate_genome(Genome * genome);

Genome * generate_population(int individuals);
int next_generation(
	Genome * parents, Genome * children,
	int n_elitism, int n_select, int n_cross, int n_new, double p_mutation, int mutation_bit
);

/*
 * Copies the information in the input genome to the output genome.
 * Note that the fitness function is also copied over.
 */
void copy_genome(Genome * in, Genome * out);

// a.k.a crossover
void tinder(Genome * population, int pop_size, Genome * last_individual);
// actual crossover
void crossover_genomes(
	Genome * gen1in,
	Genome * gen2in,
	Genome * gen1out,
	Genome * gen2out
);

void mutate_genome(Genome * genome, double p_mut);

/*
 * Implementation of elitism operator.
 * This selects the best of the best individuals to be passed over to the next generation.
 *
 * @param population the initial population with fitness values calculated
 * @param pop_size the total population size
 * @param number_elitism the number of best individuals to be selected
 * @param out the population to be filled with the best individuals
 *
 * @return the position of the best individual
 */


int extinction(int ek, Genome * population, Genome * survivors, int pop_size, int number_survivors);




void elitism(Genome * population, int pop_size, int number_elitism, Genome * out);

/*
 * Implementation of roulette wheel method.
 * This method chooses the individuals according to their fidelity, individuals with
 * lower fidelity have a lower change of being selected.
 *
 * @param population the initial population with fitness values calculated
 * @param pop_size the total population size
 * @param best_genomes the number of best individuals to choose
 * @param out the population to be filled with the best individuals
 */
int casting(Genome * population, int pop_size, int best_genomes, Genome * out);

/*
 * Add `n_new` individuals to the start of the passed genome array.
 */
void migration(
	Genome * genome,
	int n_new
);

# endif
