# include <stdio.h>
# include <omp.h>
# include "sim.h"
# include "ga.h"
# include "gradient.h"
# include "utils.h"
#include <omp.h>

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

	fprintf(fp, "fitness: %.16f ,E: %.16f  ,I_1: %.16f  ,A: %.16f \n", population[bestindividual].fitness, ic[1], ic[2], ic[3]);
	fprintf(fp, "beta: %.16f  ,phi: %.16f  ,epsilon_i: %.16f \nepsilon_Y: %.16f  ,sigma: %.16f  ,gamma_1: %.16f \ngamma_2: %.16f  ,kappa: %.16f  ,p: %.16f \nalpha: %.16f  ,delta: %.16f",
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

	store_trajectory(ic, pbest, fp);
	fclose(fp);


}
///////////////////////////////////////////////////////////////////////////////////////
// Main Function //////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv) {
	int individuals = 500;
	if(argc > 1) individuals = atoi(argv[1]);
	int maxiter = 2000; // ficar la possibilitat de donarho en runtime
	if(argc > 2) maxiter = atoi(argv[2]);
	printf("Initializing with %d individuals and %d maxiter\n", individuals, maxiter);
	int termination=0;
	int iter=0;
	double fitness_threshold = 1;

	int number_elitism = 2;
	int number_selection = 100;
	int number_crossover = 200;
	int number_migration = individuals - number_elitism - number_selection - number_crossover;
	int best_individual;

	Genome * population;
	Genome * temp_population;
	temp_population = (Genome *) malloc(individuals * sizeof(Genome));

	fitness_func ff;
	ff = fitness_exp;

	init_rng();

    printf("Generating initial population\n");
	population = generate_population(individuals);

	int i;

	printf("Entering genetic algorithm\n");
	do {
	 	// fitness calculates the fitness of every guy in the population
		#pragma omp parallel for
        for (i = 0; i < individuals; i++) // exclude those simulation of repeated genes to speed up simulation!
            if (population[i].fitness < 0) compute_fitness(population + i, ff); // TODO: parallel

		// save_population(population, individuals, "step_" + i); pero esta bé

		// elitism takes the x bests and puts them to the new generation

		// TODO: change the fitness by setting the value of ff to any of fitness_exp, fitness_max...
		best_individual = next_generation(population, temp_population,
										  number_elitism, number_selection, number_crossover, number_migration, 0.3);

		if(iter % (maxiter/100) == 0) {
			#pragma omp parallel for
			for (i = 0; i < number_elitism + number_selection; i ++)
				if (temp_population[i].fitness > 0) optimise_parameters(temp_population + i, ff);
			// copy_genome(population + best_individual, temp_population);
			printf("Generation %d with fitness %.8f\n", iter,population[best_individual].fitness);
		}

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
	} while(termination == 0);

	save_bestind(population,best_individual);
	printf("Exited genetic algorithm\n");
	printf("Fitness reached of %f, and total iterations of %d\n", population[best_individual].fitness, iter);
	free_rng();
}
