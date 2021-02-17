# include <stdio.h>
# include <omp.h>
# include "sim.h"
# include "ga.h"
# include "gradient.h"
# include "utils.h"
#include <omp.h>


/*
 * Function saving the best individual results
 *
 * @param population whole population
 * @param bestindividual position of the best individual
 */
void save_bestind(Genome * population, int bestindividual){
	FILE *fp;
	double * ic;
	ic = (double *) malloc(CoreModelDIM * sizeof(double));
	genotype_to_phenotype(population + bestindividual, ic);

	if ((fp = fopen("bestindividual.txt", "w")) != 0)
		printf("Could not open file");

	fprintf(fp, "fitness: %.16f\nE: %.16f\nI_1: %.16f\nA: %.16f\n", population[bestindividual].fitness, ic[1], ic[2], ic[3]);
	fprintf(fp, "beta: %.16f\nphi: %.16f\nepsilon_i: %.16f\nepsilon_Y: %.16f\nsigma: %.16f\ngamma_1: %.16f\ngamma_2: %.16f\nkappa: %.16f\np: %.16f\nalpha: %.16f\ndelta: %.16f",
				population->c2[0],
				population->c2[1],
				population->c2[2],
				population->c2[3],
				population->c2[4],
				population->c2[5],
				population->c2[6],
				population->c2[7],
				population->c2[8],
				population->c2[9],
				population->c2[10]);

	store_trajectory(ic, population->c2, fp);
	fclose(fp);
	free(ic);
}
/*
 *Function printing out the genome
 */
void printf_genome(Genome * g) {
	printf("%.8f,%d,%d,%d", g->fitness, g->c1[0], g->c1[1], g->c1[2]);
	printf("%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",
			g->c2[0],
			g->c2[1],
			g->c2[2],
			g->c2[3],
			g->c2[4],
			g->c2[5],
			g->c2[6],
			g->c2[7],
			g->c2[8],
			g->c2[9],
			g->c2[10]);
}

///////////////////////////////////////////////////////////////////////////////////////
// Main Function //////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////


int main(int argc, char ** argv) {

	//set up of the genetic algorithm
	int individuals = 200;
	if(argc > 1) individuals = atoi(argv[1]);
	int maxiter = 20000; // ficar la possibilitat de donarho en runtime
	if(argc > 2) maxiter = atoi(argv[2]);
	printf("Initializing with %d individuals and %d maxiter\n", individuals, maxiter);
	int iter = 0;


	// parameters for the next generation calculation
	int number_elitism = 2;
	int number_selection = 40;
	int number_crossover = 70;
	int number_migration = individuals - number_elitism - number_selection - number_crossover;
	int best_individual;
	double scale_factor = 0.1;


	//parameters for the extinctions, similarly has the parameters for the next generation calculation after extinction
	int cooldown = 200;
	unsigned int extinction_period = 500;
	int number_survivors = 2;
	int extinc_migration = 0;
	int extinc_selection = 50;
	int extinc_cross     = individuals - number_survivors - extinc_migration - extinc_selection;

	Genome * population;
	Genome * temp_population;
	temp_population = (Genome *) malloc(individuals * sizeof(Genome));

	//fitness used
	fitness_func ff;
	ff = fitness_uniform;

	//some auxiliar variables
	int ek = 0;
	int recovery = cooldown + 1 ;
	double fitness_temp=1.f;
	float epsilon = 1.0;

	//random initializer
	init_rng();

    printf("Generating initial population\n");
	population = generate_population(individuals);

	int i;

	printf("Entering genetic algorithm\n");
	while (iter < maxiter) {
	 	// fitness calculates the fitness of every guy in the population
		#pragma omp parallel for
        for (i = 0; i < individuals; i++) {// exclude those simulation of repeated genes to speed up simulation!
            if (population[i].fitness < 0) compute_fitness(population + i, ff); // TODO: parallel
		}

		// every some steps the individuals are optimized with a deterministic optimisation
		if (iter % 1000 == 0) {
			#pragma omp parallel for
			for (i = 0; i < 100; i ++)
				if (population[i].fitness > 0) optimise_parameters(population + i, ff);
	
		}

		//extinction routines
		if (recovery < cooldown) {
				best_individual = next_generation(population, temp_population,
	                                          	number_survivors, extinc_selection, extinc_cross, extinc_migration,
												0.3, scale_factor * (1 - ((double)iter) / ((double)maxiter)));
				recovery++;
				fitness_temp=population[best_individual].fitness;

			} else {
	        	int rdn=random_int(extinction_period);
	        	if (rdn < ek) {
	        		recovery=0;
					change_seed();
	        		ek = extinction( ek, population, temp_population, individuals, number_survivors);
	        		printf("An extinction has occurred\n");
	        	} else {
	        		//if no extinction happens a normal next generation is done
					best_individual = next_generation(population, temp_population,
	                                          		number_elitism, number_selection, number_crossover, number_migration,
													0.5, scale_factor * (1 - ((double)iter) / ((double)maxiter)));
					
					//if the fitness is not improving the chances of extinction rise
					if ((abs(population[best_individual].fitness -fitness_temp) < epsilon )
						||((population[best_individual].fitness -fitness_temp)==0)) ek++;
					fitness_temp=population[best_individual].fitness;
				}
			}


		if (iter % (maxiter/100) == 0)
			printf("Generation %d with fitness %.8f\n", iter, population[best_individual].fitness);

        // exchange pointers of parents and children populations
        Genome * tmp = population;
        population = temp_population;
        temp_population = tmp;
		++iter;
	}

	printf_genome(temp_population + best_individual);
	save_bestind(temp_population, best_individual);
	printf("Exited genetic algorithm\n");
	printf("Fitness reached of %f, and total iterations of %d\n", temp_population[best_individual].fitness, iter);
	free_rng();
}
