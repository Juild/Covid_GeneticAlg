# include "sim.h"
# include "ga.h"
# include "utils.h"

///////////////////////////////////////////////////////////////////////////////////////
// Main Function //////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv) {
	int individuals = 100;
	int maxiter = 10;
	int pop_size = 100;
	int termination=0;
	int iter=0;
	double fitness_threshold = 1;

	int number_elitism = (int) (0.1 * pop_size);
	int number_selection = (int) (0.6 * pop_size);
	int number_crossover = pop_size - number_elitism - number_selection;
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

		// elitism takes the x bests and puts them to the new generation
		best_individual = next_generation(population, temp_population,
                                          number_elitism, number_selection, number_crossover, 0.01);

		// termination condition
		if (population[best_individual].fitness < fitness_threshold || iter > maxiter) {
			termination=1;
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
