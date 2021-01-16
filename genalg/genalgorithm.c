#include <math.h>
#include "utils.h"
#include "functions.h"
#include <gsl/gsl_rng.h>


void generate_genome(Genome * genome) {
	//generating initial states for chromosome 1
	for (int i = 0; i < 3; i++) genome -> c1[i] = random_num();
	//generating initial states for chromosome 2
	for (int i = 0; i < 11; i++) genome -> c2[i] = random_num();
	genome->fitness = -1;
	printf("Individual generated");
}



void phenotype_to_genotype(Genome * genome, chromosome1 * c1, chromosome2 * c2) {
	c1 -> E = crom2IC(genome -> c1[0]);
	c1 -> I1 = crom2IC(genome -> c1[1]);
	c1 -> A = crom2IC(genome -> c1[2]);

	c2 -> beta = crom2HSPar(genome -> c1[0]);
	c2 -> phi = crom2Par(genome -> c1[1]);
	c2 -> e1 = crom2LSPar(genome -> c1[2]);
	c2 -> eY = crom2LSPar(genome -> c1[3]);
	c2 -> sigma = crom2Par(genome -> c1[4]);
	c2 -> gamma1 = crom2Par(genome -> c1[5]);
	c2 -> gamma2 = crom2Par(genome -> c1[6]);
	c2 -> kappa = crom2LSPar(genome -> c1[7]);
	c2 -> p = crom2Par(genome -> c1[8]);
	c2 -> alpha = crom2HSPar(genome -> c1[9]);
	c2 -> delta = crom2HSPar(genome -> c1[10]);
}


Genome * generate_population(unsigned long individuals) {
	printf("Generating population...\n")
	Genome * population;

	if( (population = (Genome *) malloc(individuals * sizeof(Genome)))== NULL )
		ExitError("when allocating memory for the population",13);

	for(int i = 0; i < individuals; i++)
		generate_genome(population + i);

	printf("Population of %d generated \n",individuals);
	return population;
}

void tinder(Genome * population, int pop_size, Genome * last_individual) { //Superlike

	int ind1 = random_num() % pop_size;
	int ind2;

	while ((ind2 = random_num() % pop_size) == ind1);
	crossover_genomes(population+ind1, population+ind2, last_individual + 1, last_individual + 2);
}

//genetic functions
void mutation(unsigned int *f, double prob) {

	if (uniform() < prob) *f = (*f)^(1U << ((unsigned char) uniform()*8*sizeof(*f)));
}


void UniformCrossover(

	unsigned int p1, unsigned int p2,
	unsigned int *f1, unsigned int *f2,
	double prob) {
	unsigned char len = 8*sizeof(*f);
	unsigned int mask = 0U;
	register unsigned char i;
	for(i=0; i < len; i++) if(uniform() < prob) mask = mask | (1U << i);
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
	int val = randomnumber % 2; // exclude last position as that means no change!
	for (i = 0; i < 3; i++) {
		if (i < val) { // copy as is
			gen1out->c1[i] = gen1in->c1[i];
			gen2out->c1[i] = gen2in->c1[i];
		} else {
			gen1out->c1[i] = gen2in->c1[i];
			gen2out->c1[i] = gen1in->c1[i];
		}
	}

	val = randomnumber % 10; // exclude last position as that means no change!
	for (i = 0; i < 11; i++) {
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

/*  PUTA MERDA SI JAJAJA TORNO AL MOVIL NO EM VAN ELS CASCOS BEUNO ADEU SI
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
	for (j = 0; j < number_elitism; j++) {
		// copy the first few individuals
		copy_genome(population + j, out + j);
	}
	for (i = number_elitism; i < pop_size; i++) {
		// loop through the rest of them to find better candidates
		for (j = 0; j < number_elitism; j++) {
			if (out[j].fitness > population[pop_size].fitness) {
				copy_genome(population + i, out + j);
			}
		}
	}
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
 */
void casting(Genome * population, int pop_size, int best_genomes, Genome * out) {
	// roulette wheeeeeeeel
	double sum_p, * p_cumsum, p;
	p_cumsum = (double *) malloc(pop_size * sizeof(double));

	int i, j;
	p_cumsum[0] = 1.0 / population[0]->fitness;
	sum_p = p_cumsum[0];
	for (i = 1; i < pop_size; i++) {
		p = 1.0 / population[i]->fitness;
		p_cumsum[i] = p_cumsum[i - 1] + p;
		sum_p += p;
	}
	//normalise
	for (i = 0; i < pop_size; i++) p_cumsum[i] /= sum_p;

	for (i = 0; i < best_genomes; i++) {
		p = random_float();
		for (j = 0; j < pop_size; j++) if (p <= p_cumsum[j]) break

		copy_genome(population + j, out + i);
	}
}

/*
 * Copies the information in the input genome to the output genome.
 * Note that the fitness function is also copied over.
 */
void copy_genome(Genome * in, Genome * out) {
	memcpy(out->c1, in->c1, 3 * sizeof(double));
	memcpy(out->c2, in->c2, 11 * sizeof(double));
	out->fitness = in->fitness;
}
///////////////////////////////////////////////////////////////////////////////
									//fitnessess
void fitness(){

}



double fitness_uniform(double ** solution) {
	int i, j;
	double f = 0.0;
	for (i = 0; i < DAYS; i++)
		for (j = 0; i < 5; i++)
			f += power(solution[i][j] - DATA[i][j], 2.0);

	return f;
}

double fitness_linear(double ** solution) {
	int i, j;
	double f = 0.0;
	for (i = 0; i < DAYS; i++)
		for (j = 0; i < 5; i++)
			f += i * power(solution[i][j] - DATA[i][j], 2.0);

	return f;
}

# define EXP_NU 0.05
double fitness_exp(double ** solution) {
	int i, j;
	double f = 0.0;
	for (i = 0; i < DAYS; i++)
		for (j = 0; i < 5; i++)
			f += exp(EXP_NU * i) * power(solution[i][j] - DATA[i][j], 2.0);

	return f;
}

double fitness_max(double ** solution) {
	int i, j;
	double f, max_f = -1;
	for (i = 0; i < DAYS; i++) {
		f = 0.0;
		for (j = 0; i < 5; i++) f += power(solution[i][j] - DATA[i][j], 2.0);
		if (f > max_f) max_f = f;
	}

	return max_f;
}

// other functions

void evolution() {

}

///////////////////////////////////////////////////////////////////////////////////////
																//Main Function//
//////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char * argv) {
	int individuals = 100;
	int maxiter = 10;
	int pop_size = 100;
	int termination=0;
	int iter=0;
	double fitness_threshold=10;

	int number_elitism = (int) (0.1 * pop_size);
	int number_selection = (int) (0.6 * pop_size);
	int number_crossover = pop_size - number_elitism - number_crossover;


	printf("Initializing population\n");
	Genome * population;
	Genome * temp_population;
	temp_population = (Genome *) malloc(individuals * sizeof(Genome));

	population = generate_population(individuals);



	printf("Entering genetic algorithm\n");
	do {

	 	// fitness calculates the fitness of every guy in the population

		// elitism takes the x bests and puts them to the new generation
		elitism(population, pop_size, number_elitism, temp_population);

		// casting  sorts the people surviving by order of fitness
		casting(population, pop_size, number_selection, temp_population + number_elitism);

		// crossover crosses the survivals to get new better individuals
		Genome * last_individual = temp_population + number_elitism + number_selection;
		for(i = 0; i < number_crossover / 2; i++) // at each iter, two individuals are added
			tinder(temp_population, number_elitism + number_selection, last_individual);
		// mutations mutate a bit some people so a bit of randomness is included
		for (i = number_elitism; i < pop_size; i++) {
			// exclude elitist from mutation! this is rather artificial
		}


		// termination condition
		if (fitness < fitness_threshold || iter > maxiter) {
			termination=1;
			if(fitness<fitness_threshold)
				printf("Exiting genetic algorithm by fitness threshold\n");
			if(iter>maxiter)
				printf("Exiting genetic algorithm by maxiter reached\n");
		}

	} while(termination == 0);


	printf("Exited genetic algorithm\n");
	printf("Fitness reached of %d, and total iterations of %d",fitness,iter);
	gsl_rng_free(rng_generator);
}
