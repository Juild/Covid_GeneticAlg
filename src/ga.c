#include "ga.h"

void generate_genome(Genome * genome) {
	//generating initial states for chromosome 1
	for (int i = 0; i < GENES_C1; i++) genome -> c1[i] = random_ulong();
	//generating initial states for chromosome 2
	for (int i = 0; i < GENES_C2; i++) genome -> c2[i] = random_ulong();
	genome->fitness = -1;
	printf("Individual generated");
}

Genome * generate_population(int individuals) {
	Genome * population;

	if( (population = (Genome *) malloc(individuals * sizeof(Genome)))== NULL )
		exit_error("when allocating memory for the population",13);

	for(int i = 0; i < individuals; i++) generate_genome(population + i);

	printf("Population of %d generated \n",individuals);
	return population;
}

void tinder(Genome * population, int pop_size, Genome * last_individual) { //Superlike
	int ind1 = random_int(pop_size);
	int ind2;

	while ((ind2 = random_int(pop_size)) == ind1);

	crossover_genomes(population + ind1, population + ind2, last_individual, last_individual + 1);
}

//genetic functions
int bitwise_mutation(unsigned long *f, double prob) {
	if (random_double() < prob) {
		*f = (*f)^(1U << ((unsigned char) random_double()*8*sizeof(*f)));
		return 1;
	}
	return 0;
}

void mutate_genome(Genome * genome, double prob) {
	int i;
	for (i = 0; i < GENES_C1; i++)
		if (bitwise_mutation(genome->c1 + i, prob)) // gen `i` has been mutate, so reset fitness to default
			genome->fitness = -1;
	for (i = 0; i < GENES_C2; i++)
		if (bitwise_mutation(genome->c2 + i, prob)) // gen `i` has been mutate, so reset fitness to default
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
	Genome * gen1out,
	Genome * gen2out
) {
	// cross initial conditions
	int i;
	int val = random_int(GENES_C1 + GENES_C2);
	for (i = 0; i < GENES_C1; i++) {
		if (i < val) { // copy as is
			gen1out->c1[i] = gen1in->c1[i];
			gen2out->c1[i] = gen2in->c1[i];
		} else {
			gen1out->c1[i] = gen2in->c1[i];
			gen2out->c1[i] = gen1in->c1[i];
		}
	}

	val -= GENES_C2;
	for (i = 0; i < GENES_C2; i++) {
		if (i < val) { // copy as is
			gen1out->c2[i] = gen1in->c2[i];
			gen2out->c2[i] = gen2in->c2[i];
		} else {
			gen1out->c2[i] = gen2in->c2[i];
			gen2out->c2[i] = gen1in->c2[i];
		}
	}

	gen1out->fitness = -1;
	gen2out->fitness = -1;
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
	double sum_p, * p_cumsum, p, best_p = -1;
	p_cumsum = (double *) malloc(pop_size * sizeof(double));

	int i, j, best_index;
	p_cumsum[0] = 1.0 / population[0].fitness;
	sum_p = p_cumsum[0];
	for (i = 1; i < pop_size; i++) {
		p = 1.0 / population[i].fitness;
		p_cumsum[i] = p_cumsum[i - 1] + p;
		sum_p += p;
		if (p > best_p) {
			best_p = p;
			best_index = j;
		}
	}
	//normalise
	for (i = 0; i < pop_size; i++) p_cumsum[i] /= sum_p;

	for (i = 0; i < best_genomes; i++) {
		p = random_double();
		for (j = 0; j < pop_size; j++) if (p <= p_cumsum[j]) break;
		copy_genome(population + j, out + i);
	}

	return best_index;
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

/*
 * Copies the information in the input genome to the output genome.
 * Note that the fitness function is also copied over.
 */
void copy_genome(Genome * in, Genome * out) {
	memcpy(out->c1, in->c1, GENES_C1 * sizeof(double));
	memcpy(out->c2, in->c2, GENES_C2 * sizeof(double));
	out->fitness = in->fitness;
}

int next_generation(
	Genome * parents, Genome * children,
	int n_elitism, int n_select, int n_cross, int n_new, double p_mutation
) {
	int pop_size = n_elitism + n_select + n_cross + n_new;
	int i;

	if (n_cross % 2 == 1) // ensure ǹumber of crossover is even
		exit_error("Please define the crossover number as even!", 33);

	// elitism takes the x bests and puts them to the new generation
	if (n_elitism > 0) elitism(parents, pop_size, n_elitism, children);

	// casting selects the best individuals randomly by fitness
	int best_individual = casting(parents, pop_size, n_select, children + n_elitism);

	// crossover crosses the survivals to get new better individuals
	Genome * last_individual = children + n_elitism + n_select;
	int n_cross_half = n_cross / 2;
	for (i = 0; i < n_cross_half; i++) {
		// at each iter, two individuals are added
		tinder(children, n_elitism + n_select, last_individual);
		last_individual += 2;
	}

	// mutations mutate a bit some people so a bit of randomness is included
	// exclude elitist from mutation! this is rather artificial
	if (p_mutation > 0) for (i = n_elitism; i < pop_size - n_new; i++) mutate_genome(children + i, p_mutation);

	if (n_new > 0) migration(children + (pop_size - n_new), n_new);

	return best_individual;
}
