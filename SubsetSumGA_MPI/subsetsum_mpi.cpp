/**
* This program is a distributed memory implementation of SubsetSum problem using Message Passing Interface.
* I have parallelized this in rather unusual domain decomposition contrary to existing GA paralleization 
* techniques employed which is population wise domain decomposition. Instead of going for this, I used
* Chromosome wise domain decomposition and column wise division of population.
* 
* 
* Author : Jahanzeb Maqbool
* Ajou University, Dept. of Computer Engineering, South Korea.
*/

#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <fstream>
#include <mpi.h>

#define root 0
#define target 25000
#define NUM_GEN 500
#define POPULATION_SIZE 40
#define _DEBUG true

int ReadsetSizeLength (char* inputFile) 
{
    FILE* file;	
    int rv ;
    int setSize = 0;		
    // Read input data from file
    file = fopen(inputFile, "r");
    if(file == NULL) {
		printf("ERROR: Unable to open file `%s'.\n", inputFile);
		exit(1);
    }
    rv = fscanf(file, "%i", &setSize);
    if(rv != 1) {
		printf("ERROR: Unable to read from file `%s'.\n", inputFile);
		fclose(file);
		exit(1);
    }
	fclose (file);
	return setSize; 
}

void ReadRandValues (char* inputFile, int setSize, int* randomSet)
{
	FILE* file;	
    int i ;
    int rv ;
    file = fopen(inputFile, "r");
    if(file == NULL) {
		printf("ERROR: Unable to open file `%s'.\n", inputFile);
		exit(1);
    }
	int garbageVariable;
    int err = fscanf(file, "%d", &garbageVariable);
	if (_DEBUG) printf ("garbage <ignore> %d err: %d\n", garbageVariable, err); 
	for ( i = 0; i < setSize; ++i ) {
		rv = fscanf(file, "%d", &randomSet[i]);
    }
	rv = fclose(file);

    if(rv != 0) {
		printf("ERROR: Unable to close file `%s'.\n", inputFile);
 	    exit(1);
    }
}

void Print(int* arr, int len) {
	printf ("{");
	for (int i = 0; i < len; i++) {
		printf ("%d", arr[i]);
		if (i != len-1) printf (", ");
	}
	printf ("}\n");
}
void PrintL(long* arr, int len) {
	printf ("{");
	for (int i = 0; i < len; i++) {
		printf ("%li", arr[i]);
		if (i != len-1) printf (", ");
	}
	printf ("}\n");
}


void InitializePopulation (int* population, int setSize) {

	for (int i = 0; i < POPULATION_SIZE; i++) {
		for (int j = 0; j < setSize; j++) {
			population [i*setSize+j] = 0.5 > ((float) rand()) / (float) RAND_MAX ? 1 : 0;
			if (_DEBUG) printf ("%d, ", population [i*setSize+j]);
		}
		if (_DEBUG) printf ("\n");
	}
	
}

void DistributePopulation (int* population, int* subPopulation, int setSize, int chunkSize, int proc) {
	
	int i, j;
	int r      = 0;
	int c      = 0;
	int limit  = 0;
	int offset = 0;
	
	if (_DEBUG) printf ("========\n");
	for (int p = 0; p < proc; p++) {
		limit += chunkSize;
		r = 0;
		for (i = 0; i < POPULATION_SIZE; i++, r++) {
			c = 0;
			for (j = offset; j < limit; j++, c++) {
				subPopulation [r*chunkSize+c] = population [i*setSize+j];
			}
		}
		offset = limit;
		if (_DEBUG) printf ("========\n");
		if (_DEBUG) Print (subPopulation, chunkSize*POPULATION_SIZE);
		
		MPI_Send(subPopulation, chunkSize*POPULATION_SIZE, MPI_INT, p, 101, MPI_COMM_WORLD);
	}
}

void SubPopulationSummation (int* subPopulation, int chunkSize, int* subSet, long* sumVec_p) {
	
	for (int i = 0; i < POPULATION_SIZE; i++) sumVec_p [i] = 0;
	for (int r = 0; r < POPULATION_SIZE; r++) {
		for (int c = 0; c < chunkSize; c++) {
			if (subPopulation[r*chunkSize+c] == 1) 
				sumVec_p [r] += subPopulation[r*chunkSize+c] * subSet[c];
		}
		
	}
}


void GetSingleChromosome (int* population, int* chrom1, int setSize, int rand) {

	int start = rand*setSize;
	int end   = start + setSize; 
	int idx = 0;
	for (int a = start; a < end; a++, idx++) chrom1 [idx] = population [a];				
	if (_DEBUG) printf ("random %d ",rand); 
}

long GetFitness (long geneVecSum) {

	// if target is greater than or equal to sum, then fitness=diff, otherwise sum itself
	// otherwise
	
	return (target >= geneVecSum)  ? target-geneVecSum : geneVecSum;
/*
	long sum = geneVecSum;
	long s = ((target - sum) >= 0) ? 1 : 0;
	return s*(target - sum) + (1 - s)*sum;		
*/
}

// Uniform crossover : Winners produce offsprings
void ProduceOffsprings (int* winner1, int* winner2, int* child1, int* child2, int setSize) {
	
	for (int i = 0; i < setSize; i++) {
		bool b = ((float) rand()) / (float) RAND_MAX >= 0.5;
		
		if (_DEBUG) printf ("boolean %d\n", b);
		if (b) {
			child1[i] = winner1[i];
			child2[i] = winner2[i];
		} else {
			child1[i] = winner2[i];
			child2[i] = winner1[i];
		}
	}
}

// Mutation : Single bit flip (gene flipping)
void PerformMutation (int* chrom, int setSize) {
	int i = rand () % setSize;
	chrom[i] = chrom[i] == 1 ? 0 : 1; // flip
}

long GetSingleChromosomeSummation (int* chrom, int* randomSet, int setSize) {
	long sum = 0;
	for (int c = 0; c < setSize; c++) {
		if (chrom[c] == 1) 
			sum += chrom[c] * randomSet[c];
	}
	return sum;
}

// Only for final result testing.
void GetAllChromSummation (int* population, int* randomSet, long* sumVec, int setSize) {
	long sum = 0;
	for (int r = 0; r < POPULATION_SIZE; r++) {
		for (int c = 0; c < setSize; c++) {
			if (population[r*setSize+c] == 1)  {
				sum += population[r*setSize+c] * randomSet[c];
			}
		}
		sumVec [r] = sum;
		sum = 0;
	}
}
 // Least fit chromosomes will be replaced. Let the best ones survive.
void ReplaceChromInPop (int* population, int* chrom, int setSize, int idx) {

	for (int i = 0; i < setSize; i++) {
		population [idx*setSize+i] = chrom [i];
	}
}
 
void ProduceNextGen (int* population, long* sumVec_p, int* randomSet, int setSize) {
	
	int idx = 0;
	
	while (idx < POPULATION_SIZE) {
	
		int i = rand () % POPULATION_SIZE;
		
		// 4-distinct integers to get distinct chromosomes
		int j, k, l;
		j = k = l = i;
		while (j==i) j = rand () % POPULATION_SIZE;
		while (k==i || k==j) k = rand () % POPULATION_SIZE;
		while (l==i || l==j || k==l) l = rand () % POPULATION_SIZE;
		if (_DEBUG) printf ("i, j, k, l : %d, %d, %d, %d\n", i, j, k, l);
		
		int* chrom1 = (int*) malloc (setSize*sizeof(int));
		int* chrom2 = (int*) malloc (setSize*sizeof(int));
		int* chrom3 = (int*) malloc (setSize*sizeof(int));
		int* chrom4 = (int*) malloc (setSize*sizeof(int));
		
		GetSingleChromosome (population, chrom1, setSize, i);
		GetSingleChromosome (population, chrom2, setSize, j);
		GetSingleChromosome (population, chrom3, setSize, k);
		GetSingleChromosome (population, chrom4, setSize, l);
		
		if (_DEBUG) Print (chrom1, setSize); 
		if (_DEBUG) Print (chrom2, setSize); 
		if (_DEBUG) Print (chrom3, setSize); 
		if (_DEBUG) Print (chrom4, setSize);
		
		if (_DEBUG) printf ("f_i, f_j, f_k, f_l : %li, %li, %li, %li \n", GetFitness (sumVec_p [i]),
				            GetFitness (sumVec_p [j]), GetFitness (sumVec_p [k]), GetFitness (sumVec_p [l]));
		
		// Selection through tournament selection
		int* winner1; 
		int* winner2; 
		if (GetFitness (sumVec_p [i]) < GetFitness (sumVec_p [j])) 
			winner1 = chrom1;		
		else 
			winner1 = chrom2;
		if (GetFitness (sumVec_p [k]) < GetFitness (sumVec_p [l])) 
			winner2 = chrom3;		
		else 
			winner2 = chrom4;
				
		if (_DEBUG) Print (winner1, setSize);
		if (_DEBUG) Print (winner2, setSize);
		
		int* child1 = (int*) malloc (setSize*sizeof (int));
		int* child2 = (int*) malloc (setSize*sizeof (int));
		// uniform crossover : winners produce offsprings
		ProduceOffsprings (winner1, winner2, child1, child2, setSize);
		if (_DEBUG) Print (child1, setSize);
		if (_DEBUG) Print (child2, setSize);
		
		// Mutation : Single Random Bit inversion
		double mutatePercent = 0.01;
		bool m1 = ((float) rand()) / (float) RAND_MAX <= mutatePercent;
		bool m2 = ((float) rand()) / (float) RAND_MAX <= mutatePercent;
		if (m1) PerformMutation (child1, setSize); 
		if (m2) PerformMutation (child2, setSize); 
		
		// Evaluate Children
		long sum_c1 = GetSingleChromosomeSummation (child1, randomSet, setSize);
		long sum_c2 = GetSingleChromosomeSummation (child2, randomSet, setSize);
		long sum_w1 = GetSingleChromosomeSummation (winner1, randomSet, setSize);
		long sum_w2 = GetSingleChromosomeSummation (winner2, randomSet, setSize);
		
		bool isChild1Good = GetFitness (sum_c1) < GetFitness (sum_w1);
		bool isChild2Good = GetFitness (sum_c2) < GetFitness (sum_w2);	
		
		//replace the chrom at (newPopSize)-index of population with two children
		if (isChild1Good) ReplaceChromInPop (population, child1, setSize, idx); 
		else ReplaceChromInPop (population, winner1, setSize, idx); 
		idx++;
		if (isChild2Good) ReplaceChromInPop (population, child2, setSize, idx); 
		else ReplaceChromInPop (population, winner2, setSize, idx); 
		idx++;

		if (_DEBUG) printf ("IDX %d\n", idx);
		
		free (child1); free (child2); //free (winner1); free (winner2); 
		free (chrom1); free (chrom2); free (chrom3); free (chrom4);   
	}

}

int main (int argc, char *argv[])
{
	char * inputFile  = argv [1]; 
	int  * values;
	int  * randomSet;
	int  * population;
	int  * subPopulation;
	int  * subSet;
	long * sumVec_p;
	int setSize;
	int chunkSize;
	int offset;
	int limit;
	int rank, proc;
	
	// execution time vars.
	double elapsedTime_p;
		
	// time seed for random generation
	time_t seconds;
	time(&seconds);
	srand((unsigned int) seconds);
	 
	// MPI initializations. 
	MPI_Init (&argc, &argv);	
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &proc);	 
	MPI_Status status;
	
	// all proc need to know the set size	
	setSize = ReadsetSizeLength (inputFile);
	if (_DEBUG) printf ("i am <%d> out of [%d] and I have read set size = %d\n", rank, proc, setSize);
	
	// DEBUG : I doubt if its necessary or not. For the time being its ommited.
	//MPI_Barrier (MPI_COMM_WORLD);
	
	// The chromosome length is same as the set length.
	chunkSize       = setSize/proc; 
	subPopulation   = (int*) malloc((chunkSize*POPULATION_SIZE)*sizeof(int));
	subSet          = (int*) malloc((chunkSize)*sizeof(int));
	sumVec_p        = (long*) malloc((POPULATION_SIZE)*sizeof(long));
	// Only root significant initializaitons.
	if (rank == 0) {
		randomSet     = (int*) malloc(setSize*sizeof(int));
		population    = (int*) malloc((setSize*POPULATION_SIZE)*sizeof(int));
		ReadRandValues (inputFile, setSize, randomSet);
		InitializePopulation (population, setSize);
		if (_DEBUG) Print (randomSet, setSize);
	}
	// Scatter the Random set into subSets
	MPI_Scatter(randomSet, chunkSize, MPI_INT, 
               subSet, chunkSize, MPI_INT, root, 
               MPI_COMM_WORLD);
	
	if (_DEBUG) Print (subSet, chunkSize);
	
	elapsedTime_p = MPI_Wtime ();
	double tTime = 0.0;
	for (int gen = 0; gen < NUM_GEN; gen++) {
		if (_DEBUG) printf ("Generation = %d\n", gen);
		if (rank == 0) {
			DistributePopulation (population, subPopulation, setSize, chunkSize, proc);
		}	
		MPI_Recv(subPopulation, chunkSize*POPULATION_SIZE, MPI_INT, root, 101, MPI_COMM_WORLD, &status);
		
		// every process calculates its local subpopulation's fitness summation
		SubPopulationSummation (subPopulation, chunkSize, subSet, sumVec_p);
		
		/* Ring communication pattern to rotate the sum of all the processors and gather it
		 * to root with the sum of all the parts (subsums) to make it a larger sum.
		 */ 
		
		if (proc > 1) {
			
			if (rank == root) {
			
				long* tempVec = (long*) malloc((POPULATION_SIZE)*sizeof(long)); 
				for (int i = 0; i < POPULATION_SIZE; i++){
					tempVec [i] = sumVec_p [i];
				}
				// blocking receive.
				MPI_Recv(sumVec_p, POPULATION_SIZE, MPI_LONG, root+1, 100, MPI_COMM_WORLD, &status);	
				// add both vectors
				for (int i = 0; i < POPULATION_SIZE; i++) {
					sumVec_p [i] += tempVec [i];
				}
				// free the temp vector
				free (tempVec);	
			}
			else if (rank == proc-1) {
				MPI_Send(sumVec_p, POPULATION_SIZE, MPI_LONG, rank-1, 100, MPI_COMM_WORLD);
			}
			// intermediate processors.
			else {
				
				long* tempVec = (long*) malloc((POPULATION_SIZE)*sizeof(long)); 
				for (int i = 0; i < POPULATION_SIZE; i++){
					tempVec [i] = sumVec_p [i];
				}
				// Recv the incoming sum vector
				MPI_Recv(sumVec_p, POPULATION_SIZE, MPI_LONG, rank+1, 100, MPI_COMM_WORLD, &status);	
				// add both vectors
				for (int i = 0; i < POPULATION_SIZE; i++) {
					sumVec_p [i] += tempVec [i];
				}
				// free the temp vector
				free (tempVec);
				// send to the previous node.
				MPI_Send(sumVec_p, POPULATION_SIZE, MPI_LONG, rank-1, 100, MPI_COMM_WORLD);
			}
		}
		
		if (rank == 0) {
			double t1 = MPI_Wtime ();
			ProduceNextGen (population, sumVec_p, randomSet, setSize);
			tTime += MPI_Wtime ()-t1;
		} 
		
		// synchronize before proceeding to the next generation.
		// DEBUG: I doube if its necessary or not. For the time being its ommitted.
		//MPI_Barrier (MPI_COMM_WORLD);
	}
	
	// Get each processor's elapsed time
	elapsedTime_p = MPI_Wtime () - elapsedTime_p;
	// Reduced time - On Head node
	double* recvbuf;
	double* sendbuf = (double*) malloc (1*sizeof(double));
	if (rank == 0) 
		recvbuf = (double*) malloc (1*sizeof(double));
	
	sendbuf [0] = elapsedTime_p;
	
	MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, 
               MPI_SUM, root, MPI_COMM_WORLD);
	
	if (_DEBUG) printf ("Total time taken on <%d> : %f\n", rank, elapsedTime_p);
	
	//if (_DEBUG) {
		if (rank == 0) {
			
			long* sumVec = (long*) malloc (POPULATION_SIZE*sizeof(long));
			
			GetAllChromSummation (population, randomSet, sumVec, setSize);
			
			printf ("Target = %d, Num of Generations = %d\n", target, NUM_GEN);
			for (int i = 0; i < POPULATION_SIZE; i++) {
				for (int j = 0; j < setSize; j++) {
					if (_DEBUG) printf ("%d ", population [i*setSize+j]);
				}
					printf (": sum = %li, fitness = %li\n", sumVec[i], GetFitness (sumVec[i]));
			}
			printf ("%f\n", tTime);
			free (sumVec);
		}	
	//}
	
	
	free (sendbuf);
	if (rank == 0) {
		//if (_DEBUG) 
		printf ("Average Execution Time of all processors: %f s\n", recvbuf[0] / proc);
		
		free (recvbuf);
		free (population);
		free (randomSet);
	}
   free (subPopulation);
   free (subSet);
   free (sumVec_p);
  
  MPI_Finalize();
  return 0;
}
