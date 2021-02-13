#include "ga.h"

void generate_genome(Genome * genome) {
	//generating initial states for chromosome 1
	genome -> c1[0] = random_int(500000000);
	genome -> c1[1] = random_int(500000000);
	genome -> c1[2] = random_int(500000000);
	//generating initial states for chromosome 2
	genome -> c2[0] = ((unsigned long) random_int(1024)) << 30;
	genome -> c2[1] = random_int(1024) << 10;
	genome -> c2[2] = random_int(1024);
	genome -> c2[3] = random_int(1024);
	genome -> c2[4] = random_int(1024) << 10;
	genome -> c2[5] = random_int(1024) << 10;
	genome -> c2[6] = random_int(1024) << 10;
	genome -> c2[7] = random_int(1024);
	genome -> c2[8] = random_int(1024) << 10;
	genome -> c2[9] = ((unsigned long) random_int(1024)) << 30;
	genome -> c2[10] = ((unsigned long) random_int(1024)) << 30;

	genome->fitness = -1;
}

Genome * generate_population(int individuals) {
	Genome * population;

	if( (population = (Genome *) malloc(individuals * sizeof(Genome)))== NULL )
		exit_error("when allocating memory for the population",13);

	for (int i = 0; i < individuals; i++) generate_genome(population + i);

	printf("Population of %d generated \n",individuals);
	return population;
}

void tinder(Genome * population, int pop_size, Genome * out) { //Superlike
	int ind1 = random_int(pop_size);
	int ind2;

	while ((ind2 = random_int(pop_size)) == ind1);

	crossover_genomes(population + ind1, population + ind2, out);
}

//genetic functions
int bitwise_mutation(unsigned long *f, double prob) {
	if (random_double() < prob) {
		*f = (*f)^(1U << random_int(UL_SIZE));
		return 1;
	}
	return 0;
}

int scaled_mutation(unsigned long *f, double prob, int max_bit) {
	if (random_double() < prob) {
		*f = (*f)^(1U << random_int(max_bit));
		return 1;
	}
	return 0;
}

void mutate_genome(Genome * genome, double prob) {
	if (scaled_mutation(genome->c1, prob, 30)
		+ scaled_mutation(genome->c1 + 1, prob, 30)
		+ scaled_mutation(genome->c1 + 2, prob, 30)
		+ scaled_mutation(genome->c2, prob, 30)
		+ scaled_mutation(genome->c2 + 1, prob, 10)
		+ scaled_mutation(genome->c2 + 2, prob, 6)
		+ scaled_mutation(genome->c2 + 3, prob, 6)
		+ scaled_mutation(genome->c2 + 4, prob, 10)
		+ scaled_mutation(genome->c2 + 5, prob, 10)
		+ scaled_mutation(genome->c2 + 6, prob, 10)
		+ scaled_mutation(genome->c2 + 7, prob, 6)
		+ scaled_mutation(genome->c2 + 8, prob, 10)
		+ scaled_mutation(genome->c2 + 9, prob, 30)
		+ scaled_mutation(genome->c2 + 10, prob, 30)) // gen `i` has been mutate, so reset fitness to default
			genome->fitness = -1;
}

void scaled_mutate_genome(Genome * genome, double prob, int max_bit) {
	if (scaled_mutation(genome->c1, prob, max_bit)
		+ scaled_mutation(genome->c1 + 1, prob, max_bit)
		+ scaled_mutation(genome->c1 + 2, prob, max_bit)
		+ scaled_mutation(genome->c2, prob, max_bit)
		+ scaled_mutation(genome->c2 + 1, prob, max_bit)
		+ scaled_mutation(genome->c2 + 2, prob, max_bit)
		+ scaled_mutation(genome->c2 + 3, prob, max_bit)
		+ scaled_mutation(genome->c2 + 4, prob, max_bit)
		+ scaled_mutation(genome->c2 + 5, prob, max_bit)
		+ scaled_mutation(genome->c2 + 6, prob, max_bit)
		+ scaled_mutation(genome->c2 + 7, prob, max_bit)
		+ scaled_mutation(genome->c2 + 8, prob, max_bit)
		+ scaled_mutation(genome->c2 + 9, prob, max_bit)
		+ scaled_mutation(genome->c2 + 10, prob, max_bit)) // gen `i` has been mutate, so reset fitness to default
			genome->fitness = -1;
}


void bitwise_crossover_uniform(
	unsigned int p1, unsigned int p2,
	unsigned int *f1, unsigned int *f2,
	double prob
) {
	unsigned char len = 8*sizeof(*f1);
	unsigned int mask = 0U;
	register unsigned char i;
	for(i=0; i < len; i++) if(random_double() < prob) mask = mask | (1U << i);
	*f1 = (p1 & mask) | (p2 & ~mask);
	*f2 = (p2 & mask) | (p1 & ~mask);
}


void crossover_genomes(
	Genome * gen1in,
	Genome * gen2in,
	Genome * out
) {
	// cross initial conditions
	int val = random_int(GENES_C1 + GENES_C2 - 1) + 1;
	if (val < GENES_C1) {
		// crossover in the first chromosome
		memcpy(out->c1, gen1in->c1, val * UL_SIZE);
		memcpy(out->c1, gen2in->c1, (GENES_C1 - val) * UL_SIZE);

		memcpy(out->c2, gen2in->c2, GENES_C2 * UL_SIZE);
	} else {
		// crossover in the 2nd chromosome
		memcpy(out->c1, gen1in->c1, GENES_C1 * UL_SIZE);

		val -= GENES_C1;
		memcpy(out->c2, gen1in->c2, val * UL_SIZE);
		memcpy(out->c2, gen2in->c2, (GENES_C2 - val) * UL_SIZE);
	}

	out->fitness = -1;
}

int extinction(int ek, Genome * population, Genome * survivors, int pop_size, int number_survivors){
	ek=0;
	int new_population =pop_size-number_survivors;

	elitism(population,pop_size, number_survivors,survivors);
	migration(survivors+number_survivors,new_population);

	return ek;

}

/*
 * Implementation of elitism operator.
 * This selects the best of the best individuals to be passed over to the next generation.
 *
 * @param population the initial population with fitness values calculated
 * @param pop_size the total population size
 * @param number_elitism the number of best individuals to be selected
 * @param out the population to be filled with the best individuals
 */
void elitism(Genome * population, int pop_size, int number_elitism, Genome * out) {
	int i, j;
	int best_index[number_elitism];

	// initialise the list of best individuals arbitrarialy to the first individuals of the previous population
	for (j = 0; j < number_elitism; j++) best_index[j] = j;

	// loop through the rest of them to find better candidates
	for (i = number_elitism; i < pop_size; i++)
		for (j = 0; j < number_elitism; j++)
			if (population[best_index[j]].fitness > population[i].fitness) {
				// we have found a better one at index i, so add it to the list of best individuals
				best_index[j] = i; break;
			}

	// copy the best individuals to the out population
	for (j = 0; j < number_elitism; j++) copy_genome(population + best_index[j], out + j);
}

/*
 * Implementation of roulette wheel method.
 * This method chooses the individuals according to their fidelity, individuals with
 * lower fidelity have a lower change of being selected.
 *
 * @param population the initial population with fitness values calculated
 * @param pop_size the total population size
 * @param best_genomes the number of best individuals to choose
 * @param out the population to be filled with the best individuals
 *
 * @return the index corresponding to the best individual
 */
int casting(Genome * population, int pop_size, int best_genomes, Genome * out) {
	// roulette wheeeeeeeel
	double sum_p, * p_cumsum, p, best_p;
	p_cumsum = (double *) malloc(pop_size * sizeof(double));

	int i, j, best_index = 0;
	p_cumsum[0] = best_p = 1.0 / population[0].fitness;
	for (i = 1; i < pop_size; i++) {
		p = 1.0 / population[i].fitness;
		p_cumsum[i] = p_cumsum[i - 1] + p;
		if (p > best_p) {
			best_p = p;
			best_index = i;
		}
	}
	//normalise
	sum_p = p_cumsum[pop_size - 1];
	for (i = 0; i < pop_size; i++) p_cumsum[i] /= sum_p;

	for (i = 0; i < best_genomes; i++) {
		p = random_double();
		for (j = 0; j < pop_size; j++) if (p <= p_cumsum[j]) break;
		copy_genome(population + j, out + i);
	}

	return best_index;
}

/* Does casting and elitism all at once to avoid double loops */
int elitist_casting(Genome * population, int pop_size, int best_genomes, int number_elitism, Genome * out) {
	// roulette wheeeeeeeel
	register int i, j;
	double p_cumsum[pop_size], p;

	int best_indices[number_elitism];

	// initialise the list of best individuals arbitrarialy to the first individuals of the previous population
	for (j = 0; j < number_elitism; j++) best_indices[j] = j;

	p_cumsum[0] = 1.0 / population[0].fitness;
	for (i = 1; i < pop_size; i++) {
		p = 1.0 / population[i].fitness;
		p_cumsum[i] = p_cumsum[i - 1] + p;
		for (j = 0; j < number_elitism; j++)
			if (population[best_indices[j]].fitness > population[i].fitness) {
				// we have found a better one at index i, so add it to the list of best individuals
				best_indices[j] = i; break;
			}
	}

	//normalise
	p = p_cumsum[pop_size - 1];
	for (i = 0; i < pop_size; i++) p_cumsum[i] /= p;

	// copy the best individuals to the out population
	for (j = 0; j < number_elitism; j++) copy_genome(population + best_indices[j], out + j);

	best_genomes += number_elitism; // just increase this number here and start after elitism to reduce operations
	for (i = number_elitism; i < best_genomes; i++) {
		p = random_double();
		for (j = 0; j < pop_size; j++) if (p < p_cumsum[j]) break;
		copy_genome(population + j, out + i);
	}

	// by how we calculated this array, the best is guaranteed to be in the first position
	return best_indices[0];
}

/*
 * Add `n_new` individuals to the start of the passed genome array.
 */
void migration(
	Genome * genome,
	int n_new
) {
	for (int i = 0; i < n_new; i++) generate_genome(genome + i);
}

void gaussian_migration(
	Genome * genome,
	int n_new,
	Genome * best
) {
	unsigned long mask = 1UL << 10;

	for (int i = 0; i < n_new; i++) {
		genome[i].c1[0] = (best->c1[0] ^ (mask << 20)) + random_int(1 << 20);
		genome[i].c1[1] = (best->c1[1] ^ (mask << 20)) + random_int(1 << 20);
		genome[i].c1[2] = (best->c1[2] ^ (mask << 20)) + random_int(1 << 20);
		genome[i].c2[0] = (best->c2[0] ^ (mask << 30)) + random_int(1 << 30);
		genome[i].c2[1] = (best->c2[1] ^ (mask << 10)) + random_int(1 << 10);
		genome[i].c2[2] = (best->c2[2] ^ (mask << 6)) + random_int(1 << 4);
		genome[i].c2[3] = (best->c2[3] ^ (mask << 6)) + random_int(1 << 4);
		genome[i].c2[4] = (best->c2[4] ^ (mask << 10)) + random_int(1 << 10);
		genome[i].c2[5] = (best->c2[5] ^ (mask << 10)) + random_int(1 << 10);
		genome[i].c2[6] = (best->c2[6] ^ (mask << 10)) + random_int(1 << 10);
		genome[i].c2[7] = (best->c2[7] ^ (mask << 6)) + random_int(1 << 4);
		genome[i].c2[8] = (best->c2[8] ^ (mask << 10)) + random_int(1 << 10);
		genome[i].c2[9] = (best->c2[9] ^ (mask << 30)) + random_int(1 << 30);
		genome[i].c2[10] = (best->c2[10] ^ (mask << 30)) + random_int(1 << 30);
		genome[i].fitness = -1;
	}
}

/*
 * Copies the information in the input genome to the output genome.
 * Note that the fitness function is also copied over.
 */
void copy_genome(Genome * in, Genome * out) {
	memcpy(out->c1, in->c1, GENES_C1 * UL_SIZE);
	memcpy(out->c2, in->c2, GENES_C2 * UL_SIZE);
	out->fitness = in->fitness;
}

int next_generation(
	Genome * parents, Genome * children,
	int n_elitism, int n_select, int n_cross, int n_new, double p_mutation, int mutation_bit
) {
	int pop_size = n_elitism + n_select + n_cross + n_new;
	int i;

	int best_individual;
	if (n_elitism > 0) // elitism takes the x bests and puts them to the new generation
		best_individual = elitist_casting(parents, pop_size, n_select, n_elitism, children);
	else // casting selects the best individuals randomly by fitness
	 	best_individual = casting(parents, pop_size, n_select, children + n_elitism);

	// crossover crosses the survivals to get new better individuals
	for (i = 0; i < n_cross; i++) tinder(children, n_elitism + n_select, children + n_elitism + n_select + i);

	// mutations mutate a bit some people so a bit of randomness is included
	// exclude elitist from mutation! this is rather artificial
	// TODO: parallel
	if (p_mutation > 0)
		for (i = n_elitism; i < pop_size - n_new; i++) mutate_genome(children + i, p_mutation);

	if (n_new > 0) {
		// gaussian_migration(children + (pop_size - n_new), n_new, parents + best_individual);
		migration(children + (pop_size - n_new), n_new);
	}

	return best_individual;
}
