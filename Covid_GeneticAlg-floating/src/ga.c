#include "ga.h"

/*
 *
 * Subroutines of the genetic algorithm calculating every next generation
 * with the evolution and genetic mechanisms
 *
 */




/*
 * Function generating a single genome in genotype form
 *
 * @param genome the genome of a single individual
 */

void generate_genome(Genome * genome) {
	//generating initial states for chromosome 1
	genome -> c1[0] = random_int(5000000);
	genome -> c1[1] = random_int(5000000);
	genome -> c1[2] = random_int(5000000);
	//generating initial states for chromosome 2
	genome -> c2[0] = random_double();
	genome -> c2[1] = random_double();
	genome -> c2[2] = random_double();
	genome -> c2[3] = random_double();
	genome -> c2[4] = random_double();
	genome -> c2[5] = random_double();
	genome -> c2[6] = random_double();
	genome -> c2[7] = random_double();
	genome -> c2[8] = random_double();
	genome -> c2[9] = random_double();
	genome -> c2[10] = random_double();

	genome->fitness = -1;
}


/*
 * Function generating the whole population
 * by doing the generate genome function
 *
 * @param individuals number of individuals of the population
 */
Genome * generate_population(int individuals) {
	Genome * population;

	if( (population = (Genome *) malloc(individuals * sizeof(Genome)))== NULL )
		exit_error("when allocating memory for the population",13);

	for (int i = 0; i < individuals; i++) generate_genome(population + i);

	printf("Population of %d generated \n",individuals);
	return population;
}

/*
 * Function deciding randomly which individuals will cross to return a new individual
 * is used only for the non elitist or fit individuals
 *
 * @param population the initial population with fitness values calculated
 * @param pop_size the total population size
 * @param out the temporal population to fill for the next generation
 */
void crossing(Genome * population, int pop_size, Genome * out) { //Superlike
	int ind1 = random_int(pop_size);
	int ind2;

	while ((ind2 = random_int(pop_size)) == ind1);

	crossover_genomes(population + ind1, population + ind2, out);
}


int bitwise_mutation(unsigned long *f, double prob) {
	//outdated function mutating the element given f
	if (random_double() < prob) {
		*f = (*f)^(1U << random_int(UL_SIZE));
		return 1;
	}
	return 0;
}

/*
 * Implementation of mutation operator with a scaling factor
 * that reduces the mutation in time
 *
 * @param f element to mutate
 * @param prob probability of the mutation
 * @param max_bit maximum bit that can be mutated
 */
int scaled_mutation(unsigned int *f, double prob, int max_bit) {
	if (random_double() < prob) {
		*f = (*f)^(1U << random_int(max_bit));
		return 1;
	}
	return 0;
}
/*
 * Implementation of mutation operator for floats
 * similar to the scaling but using a gaussian distribution
 *
 * @param f element to mutate
 * @param prob probability of the mutation
 * @param sigma stdeviation of the gaussian used for the mutation
 */
int float_mutation(double *f, double prob, double sigma) {
	if (random_double() < prob) {
		do { prob = random_gaussian(sigma); } while (prob + *f > 1 || prob + *f < 0);
		*f += prob;
		return 1;
	}
	return 0;
}
/*
 * Main function for the mutation operator
 *
 * @param genome genome of a single individual
 * @param prob probability of the mutation
 */
void mutate_genome(Genome * genome, double prob, double sigma) {
	if (scaled_mutation(genome->c1, prob, 30)
		+ scaled_mutation(genome->c1 + 1, prob, 30)
		+ scaled_mutation(genome->c1 + 2, prob, 30)
		+ float_mutation(genome->c2, prob, sigma)
		+ float_mutation(genome->c2 + 1, prob, sigma)
		+ float_mutation(genome->c2 + 2, prob, sigma)
		+ float_mutation(genome->c2 + 3, prob, sigma)
		+ float_mutation(genome->c2 + 4, prob, sigma)
		+ float_mutation(genome->c2 + 5, prob, sigma)
		+ float_mutation(genome->c2 + 6, prob, sigma)
		+ float_mutation(genome->c2 + 7, prob, sigma)
		+ float_mutation(genome->c2 + 8, prob, sigma)
		+ float_mutation(genome->c2 + 9, prob, sigma)
		+ float_mutation(genome->c2 + 10, prob, sigma)) // gen `i` has been mutate, so reset fitness to default
			genome->fitness = -1;
}




/*
 * Function for the crossover mechanism
 *
 * @param gen1in parent 1
 * @param gen2in parent 2
 * @param out children
 */
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
		memcpy(out->c1 + val, gen2in->c1 + val, (GENES_C1 - val) * UL_SIZE);

		memcpy(out->c2, gen2in->c2, GENES_C2 * UL_SIZE);
	} else {
		// crossover in the 2nd chromosome
		memcpy(out->c1, gen1in->c1, GENES_C1 * UL_SIZE);

		val -= GENES_C1;
		memcpy(out->c2, gen1in->c2, val * UL_SIZE);
		memcpy(out->c2 + val, gen2in->c2 + val, (GENES_C2 - val) * UL_SIZE);
	}

	out->fitness = -1;
}
/*
 * Implementation of extinction operator.
 * This selects a small group of individuals to be passed over to the next generation.
 *
 * @param population the initial population with fitness values calculated
 * @param pop_size the total population size
 * @param number_survivors the number of best individuals to be selected
 * @param survivors the new population which will be filled with the survivors and the new individuals
 * @param ek constant regulating the happening of extinctions
 */
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
	// roulette wheel
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
	// roulette wheel
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


/*
 * Copies the information in the input genome to the output genome.
 * Note that the fitness function is also copied over.
 */
void copy_genome(Genome * in, Genome * out) {
	memcpy(out->c1, in->c1, GENES_C1 * UL_SIZE);
	memcpy(out->c2, in->c2, GENES_C2 * UL_SIZE);
	out->fitness = in->fitness;
}

/*
 * Main function of the ga.c
 * calculates the next generation of the population with
 * all the operations of genetic algorithm
 *
 * @param parents previous step population
 * @param children next step population
 * @param n_elitism number of elitism individuals
 * @param n_select number of selected individuals
 * @param n_cross number of crossover individuals
 * @param n_new number of new individuals
 * @param p_mutation probability of mutating
 * @param mutation_bit maximum bit in where the mutation will happen
 */
int next_generation(
	Genome * parents, Genome * children,
	int n_elitism, int n_select, int n_cross, int n_new, double p_mutation, double sigma
) {
	int pop_size = n_elitism + n_select + n_cross + n_new;
	int i;

	int best_individual;
	if (n_elitism > 0) // elitism takes the x bests and puts them to the new generation
		best_individual = elitist_casting(parents, pop_size, n_select, n_elitism, children);
	else // casting selects the best individuals randomly by fitness
	 	best_individual = casting(parents, pop_size, n_select, children + n_elitism);

	// crossover crosses the survivals to get new better individuals
	for (i = 0; i < n_cross; i++) crossing(children, n_elitism + n_select, children + n_elitism + n_select + i);

	// mutations mutate a bit some people so a bit of randomness is included
	// exclude elitist from mutation! this is rather artificial
	// TODO: parallel
	if (p_mutation > 0)
		for (i = n_elitism; i < pop_size - n_new; i++) mutate_genome(children + i, p_mutation, sigma);

	if (n_new > 0) {
		// gaussian_migration(children + (pop_size - n_new), n_new, parents + best_individual);
		migration(children + (pop_size - n_new), n_new);
	}

	return best_individual;
}
