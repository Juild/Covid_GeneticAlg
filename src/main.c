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


void save_bestind(Genome * population, int bestindividual){
	FILE *fp;
	double * ic;
	ic = (double *) malloc(CoreModelDIM * sizeof(double));
	Parameters * pbest;
	pbest = (Parameters *) malloc(sizeof(Parameters));
	genotype_to_phenotype(population + bestindividual, ic, pbest);

	if ((fp = fopen("bestindividual.txt", "w")) != 0)
		printf("Could not open file");

	fprintf(fp, "fitness: %.8f ,E: %f  ,I_1: %f  ,A: %f \n", population[bestindividual].fitness, ic[1], ic[2], ic[3]);
	fprintf(fp, "beta: %.8f  ,phi: %.8f  ,epsilon_i: %.8f \nepsilon_Y: %.8f  ,sigma: %.8f  ,gamma_1: %.8f \ngamma_2: %.8f  ,kappa: %.8f  ,p: %.8f \nalpha: %.8f  ,delta: %.8f",
				pbest->beta,
				pbest->phi,
				pbest->e1,
				pbest->eY,
				pbest->sigma,
				pbest->gamma1,
				pbest->gamma2,
				pbest->kappa,
				pbest->p,
				pbest->alpha,
				pbest->delta);

	fclose(fp);


}
///////////////////////////////////////////////////////////////////////////////////////
// Main Function //////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv) {
	int individuals = 1000;
	int maxiter = 1000; // ficar la possibilitat de donarho en runtime
	int termination=0;
	int iter=0;
	double fitness_threshold = 1;

	int number_elitism = (int) (0.1 * individuals);
	int number_selection = (int) (0.5 * individuals);
	int number_migration = (int) (0.1 * individuals);
	int number_crossover = individuals - number_elitism - number_selection - number_migration;
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
                                          number_elitism, number_selection, number_crossover, number_migration, 0.01);

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
			++iter;
        }
        if(iter % (maxiter/100) == 0)printf("Generation %d with fitness %.8f\n",
        	iter,population[best_individual].fitness);
	} while(termination == 0);

	save_bestind(population,best_individual);
	printf("Exited genetic algorithm\n");
	printf("Fitness reached of %f, and total iterations of %d\n", population[best_individual].fitness, iter);
	free_rng();
}
