# include <stdio.h>
# include "sim.h"
# include "ga.h"
# include "utils.h"

/*
 * Just for statistic purposes.
 */
void save_population(Genome * population, int individuals, const char * filename) {
	FILE *fp;

	if ((fp = fopen(filename, "w")) != 0)
		printf("Could not open file");

	for (int i = 0; i < individuals; i++) {
		fprintf(fp, "%.8f,%ld,%ld,%ld",
				population[i].fitness, population[i].c1[0], population[i].c1[1], population[i].c1[2]);
		fprintf(fp, "%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld",
				population[i].c2[0],
				population[i].c2[1],
				population[i].c2[2],
				population[i].c2[3],
				population[i].c2[4],
				population[i].c2[5],
				population[i].c2[6],
				population[i].c2[7],
				population[i].c2[8],
				population[i].c2[9],
				population[i].c2[10]);
	}

	fclose(fp);
}

///////////////////////////////////////////////////////////////////////////////////////
// Main Function //////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv) {
	int individuals = 100;
	int maxiter = 10; // ficar la possibilitat de donarho en runtime
	int pop_size = 100;
	int termination=0;
	int iter=0;
	double fitness_threshold = 1;

	int number_elitism = (int) (0.1 * individuals);
	int number_selection = (int) (0.5 * individuals);
	int number_crossover = individuals - number_elitism - number_selection;
	int best_individual;

	Genome * population;
	Genome * temp_population;
	temp_population = (Genome *) malloc(individuals * sizeof(Genome));

	init_rng();

    printf("Generating initial population\n");
	population = generate_population(individuals);

	int i;

	printf("Entering genetic algorithm\n");
	do {
	 	// fitness calculates the fitness of every guy in the population
        for (i = 0; i < individuals; i++) // exclude those simulation of repeated genes to speed up simulation!
            if (population[i].fitness < 0) compute_fitness(population + i, fitness_exp);

		// save_population(population, individuals, "step_" + i); pero esta bé

		// elitism takes the x bests and puts them to the new generation
		best_individual = next_generation(population, temp_population,
                                          number_elitism, number_selection, number_crossover, 0, 0.01);

		// termination condition
		if (population[best_individual].fitness < fitness_threshold || iter > maxiter) {
			termination = 1;
			if (population[best_individual].fitness < fitness_threshold)
				printf("Exiting genetic algorithm by fitness threshold\n");
			if (iter > maxiter)
				printf("Exiting genetic algorithm by maxiter reached\n");
		} else {
            // exchange pointers of parents and children populations
            Genome * tmp = population;
            population = temp_population;
            temp_population = tmp;
        }
	} while(termination == 0);

	printf("Exited genetic algorithm\n");
	printf("Fitness reached of %f, and total iterations of %d", population[best_individual].fitness, iter);
	free_rng();
}
