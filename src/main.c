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

	int cooldown = 200;
	unsigned int extinction_period = 500;


	int individuals = 500;
	int maxiter = 2000; // ficar la possibilitat de donarho en runtime
	int termination=0;
	int iter=0;
	double fitness_threshold = 1;

	int number_elitism = 2;
	int number_selection = 20;
	int number_migration = 70;
	int number_crossover = individuals - number_elitism - number_selection - number_migration;
	int best_individual;

	int number_survivors = 4;
	int extinc_migration = 0;
	int extinc_selection = 200;
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
	do {
	 	// fitness calculates the fitness of every guy in the population
        for (i = 0; i < individuals; i++) // exclude those simulation of repeated genes to speed up simulation!
            if (population[i].fitness < 0) compute_fitness(population + i, ff); // TODO: parallel

		// save_population(population, individuals, "step_" + i); pero esta bé

		// elitism takes the x bests and puts them to the new generation

		// TODO: change the fitness by setting the value of ff to any of fitness_exp, fitness_max...

		if (recovery < cooldown) {
			best_individual = next_generation(population, temp_population,
                                          	number_survivors, extinc_selection, extinc_cross, extinc_migration, 0.1);
			recovery++;
			//if(iter % (maxiter/100) == 0)printf("Entering cooldown if\n");
			fitness_temp=population[best_individual].fitness;

		} else {
        	int rdn=random_int(extinction_period);
        	if (rdn < ek) {
        		recovery=0;
        		ek = extinction( ek, population, temp_population, individuals, number_survivors);
        		printf("An extinction has occurred\n");
        		init_rng();
        	} else {
				best_individual = next_generation(population, temp_population,
                                          		number_elitism, number_selection, number_crossover, number_migration, 0.3);
				//if(iter % (maxiter/100) == 0)printf("normal behaviour, values %.8f,%.8f\n",population[best_individual].fitness,fitness_temp);

				if ((abs(population[best_individual].fitness -fitness_temp) < epsilon )||((population[best_individual].fitness -fitness_temp)==0)) ek++;
				fitness_temp=population[best_individual].fitness;
			}
		}

		if(iter % (maxiter/100) == 0)
			printf("Generation %d with fitness %.8f and ek%d\n",
					iter,population[best_individual].fitness,ek);
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
