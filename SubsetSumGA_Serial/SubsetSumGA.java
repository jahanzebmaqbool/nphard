package com.ga.subsetsum;

/**
* This code implements subset sum problem using Java.
* This code is a modification to a general template by Kunuk Nykjaer.
* But this code varies alot from the general template as it is specifically
* solving subset sum problem.

* @author: jahanzeb
*/
import java.util.*;

public class SubsetSumGA {

	static final int NUM_GENERATIONS = 50000;
	static final int TARGET = 25000;
	static final int SET_SIZE = 1000;
	static final int CHROMOSOME_SIZE = SET_SIZE;
	static long BEGIN;
	static final boolean _DEBUG = true;
	LinkedList<Chromosome> population = new LinkedList<Chromosome>();
	final Random rand;

	final int populationSize = 40;
	int[] set;

	public SubsetSumGA() {

		set = new int[SET_SIZE];
		rand = new Random();

		for (int i = 0; i < SET_SIZE; i++) {
			set[i] = rand.nextInt(SET_SIZE);
			System.out.print(set[i] + ",");
		}

		for (int i = 0; i < populationSize; i++) {
			Chromosome c = new Chromosome(set, TARGET);
			c.random();
			population.add(c);
		}

//		Collections.sort(population); // sort method
//		System.out.println("Init population sorted");
		print();
	}

	void print() {
//		System.out.println(population.size());
		System.out.println("-- print");
		for (Chromosome c : population) {
			System.out.println(c);
		}
	}

	/**
	 * Selection strategy: Tournament method and steady state find 4 random in
	 * population not same let 2 fight, and 2 fight the winners makes 2 children
	 */
	void produceNextGen() {
		LinkedList<Chromosome> newpopulation = new LinkedList<Chromosome>();

		// while (newpopulation.size() < populationSize
		// * (1.0 - (parentUsePercent / 100.0))) {
		while (newpopulation.size() < populationSize) {
			int size = population.size();
			int i = rand.nextInt(size);
			int j, k, l;
			j = k = l = i;
			while (j == i)
				j = rand.nextInt(size);
			while (k == i || k == j)
				k = rand.nextInt(size);
			while (l == i || l == j || k == l)
				l = rand.nextInt(size);

			Chromosome c1 = population.get(i);
			Chromosome c2 = population.get(j);
			Chromosome c3 = population.get(k);
			Chromosome c4 = population.get(l);

			Chromosome w1, w2;

			if (c1.fitness() < c2.fitness())
				w1 = c1;
			else
				w1 = c2;
			if (c3.fitness() < c4.fitness())
				w2 = c3;
			else
				w2 = c4;

			Chromosome child1, child2;

			// Method uniform crossover
			Chromosome[] childs = newChilds(w1, w2);
			child1 = childs[0];
			child2 = childs[1];

			double mutatePercent = 0.01;
			boolean m1 = rand.nextFloat() <= mutatePercent;
			boolean m2 = rand.nextFloat() <= mutatePercent;

			if (m1)
				mutate(child1);
			if (m2)
				mutate(child2);
			
			// here fitness function results in opposite manner
			boolean isChild1Good = child1.fitness() < w1.fitness();
			boolean isChild2Good = child2.fitness() < w2.fitness();


			newpopulation.add(isChild1Good ? child1 : w1);
			newpopulation.add(isChild2Good ? child2 : w2);
		}

		population = newpopulation;

//		print();
//		try {
//			Thread.sleep(1000);
//		} catch (InterruptedException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
	}

//	// one-point crossover random pivot
//	Chromosome newChild(Chromosome c1, Chromosome c2, int pivot) {
//		Chromosome child = new Chromosome(set, TARGET);
//
//		for (int i = 0; i < pivot; i++) {
//			child.genotype[i] = c1.genotype[i];
//		}
//		for (int j = pivot; j < CHROMOSOME_SIZE; j++) {
//			child.genotype[j] = c2.genotype[j];
//		}
//
//		return child;
//	}

	// Uniform crossover
	Chromosome[] newChilds(Chromosome c1, Chromosome c2) {
		Chromosome child1 = new Chromosome(set, TARGET);
		Chromosome child2 = new Chromosome(set, TARGET);

		for (int i = 0; i < CHROMOSOME_SIZE; i++) {
			boolean b = rand.nextFloat() >= 0.5;
			if (b) {
				child1.genotype[i] = c1.genotype[i];
				child2.genotype[i] = c2.genotype[i];
			} else {
				child1.genotype[i] = c2.genotype[i];
				child2.genotype[i] = c1.genotype[i];
			}
		}

		return new Chromosome[] { child1, child2 };
	}

	void mutate(Chromosome c) {
		int i = rand.nextInt(SET_SIZE);
		// c.genotype[i] = !c.genotype[i]; // flip
		c.genotype[i] = c.genotype[i] == 1 ? 0 : 1; // flip
	}

	public static void main(String[] args) {
		BEGIN = System.currentTimeMillis();

		SubsetSumGA ga = new SubsetSumGA();
		ga.run();

		long END = System.currentTimeMillis();
		System.out.println("Time: " + (END - BEGIN) / 1000.0 + " sec.");
	}

	void run() {
		
		int count = 0;
		while (count < NUM_GENERATIONS) {
			count++;
			produceNextGen();
		}

		System.out.println("\nResult");
		print();
	}

	public class Chromosome {
		public int cSize;
		public int[] genotype;
		public int[] set;
		public int target;

		public Chromosome(int[] aSet, int aTarget) {
			this.set = aSet;
			this.target = aTarget;
			genotype = new int[aSet.length];

		}

		void random() {
			for (int i = 0; i < genotype.length; i++) {
				genotype[i] = 0.5 > rand.nextFloat() ? 1 : 0;
			}
		}

		private String gene() {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < genotype.length; i++) {
				sb.append(genotype[i]);
			}
			return sb.toString();
		}

		int summation() {
			int sum = 0;
			for (int i = 0; i < genotype.length; i++) {
				if (genotype[i] == 1)
					sum += genotype[i] * this.set[i];
			}
			return sum;
		}

		public int fitness () {
			
			int sum = this.summation();
			int s = ((target - sum) >= 0) ? 1 : 0;
			
			return s*(target - sum) + (1 - s)*sum;		
		}
		
		@Override
		public String toString() {
			return "gene=" + gene() + " sum=" + summation();
		}

	}
}
