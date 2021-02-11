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

	fprintf(fp, "fitness: %.16f\nE: %.16f\nI_1: %.16f\nA: %.16f\n", population[bestindividual].fitness, ic[1], ic[2], ic[3]);
	fprintf(fp, "beta: %.16f\nphi: %.16f\nepsilon_i: %.16f\nepsilon_Y: %.16f\nsigma: %.16f\ngamma_1: %.16f\ngamma_2: %.16f\nkappa: %.16f\np: %.16f\nalpha: %.16f\ndelta: %.16f",
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
	free(ic);
	free(pbest);
}

void printf_genome(Genome * g) {
	printf("%.8f,%ld,%ld,%ld", g->fitness, g->c1[0], g->c1[1], g->c1[2]);
	printf("%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",
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
	int individuals = 300;
	if(argc > 1) individuals = atoi(argv[1]);
	int maxiter = 10000; // ficar la possibilitat de donarho en runtime
	if(argc > 2) maxiter = atoi(argv[2]);
	printf("Initializing with %d individuals and %d maxiter\n", individuals, maxiter);
	int iter = 0;
	int termination = 0;

	int number_elitism = 2;
	int number_selection = 20;
	int number_crossover = 70;
	int number_migration = individuals - number_elitism - number_selection - number_crossover;
	int best_individual;
	int mutation_bit = UL_SIZE;

	int cooldown = 200;
	unsigned int extinction_period = 500;
	int number_survivors = 2;
	int extinc_migration = 0;
	int extinc_selection = 50;
	int extinc_cross     = individuals - number_survivors - extinc_migration - extinc_selection;

	Genome * population;
	Genome * temp_population;
	temp_population = (Genome *) malloc(individuals * sizeof(Genome));

	fitness_func ff;
	ff = fitness_exp;

	int ek = 0;
	int recovery = cooldown + 1 ;
	double fitness_temp=1.f;
	float epsilon = 1.0;

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

		if (iter % 200 == 0) {
			#pragma omp parallel for
			for (i = 0; i < individuals; i ++)
				if (population[i].fitness > 0) optimise_parameters(population + i, ff);
			// copy_genome(population + best_individual, temp_population);
		}

		mutation_bit = 1 + (int) (UL_SIZE - 1) * ((1.0 - ((float)iter) / ((float)maxiter)));

		if (recovery < cooldown) {
				best_individual = next_generation(population, temp_population,
	                                          	number_survivors, extinc_selection, extinc_cross, extinc_migration,
												0.1, mutation_bit);
				recovery++;
				//if(iter % (maxiter/100) == 0)printf("Entering cooldown if\n");
				fitness_temp=population[best_individual].fitness;

			} else {
	        	int rdn=random_int(extinction_period);
	        	if (rdn < ek) {
	        		recovery=0;
					change_seed();
	        		ek = extinction( ek, population, temp_population, individuals, number_survivors);
	        		printf("An extinction has occurred\n");
	        	} else {
					best_individual = next_generation(population, temp_population,
	                                          		number_elitism, number_selection, number_crossover, number_migration,
													0.3, mutation_bit);
					//if(iter % (maxiter/100) == 0)printf("normal behaviour, values %.8f,%.8f\n",population[best_individual].fitness,fitness_temp);

					if ((abs(population[best_individual].fitness -fitness_temp) < epsilon )
						||((population[best_individual].fitness -fitness_temp)==0)) ek++;
					fitness_temp=population[best_individual].fitness;
				}
			}

/*
		best_individual = next_generation(population, temp_population,
										number_elitism, number_selection, number_crossover, number_migration,
										0.1, mutation_bit);
*/
		// mutation_bit = UL_SIZE;
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
