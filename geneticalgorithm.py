import numpy as np
import random

GA_POPSIZE = 10
GA_MAXITER = 1000
GA_ELITRATE = 0.10
GA_MUTATIONRATE = 0.25


class GAStruct:
    def __init__(self, d1, d2, fitness):
        self.d1 = d1
        self.d2 = d2
        self.fitness = fitness


class GeneticAlgorithm:
    population = []
    buffer = []

    def runsimulation(self):
        self.init_population()
        for i in range(0, GA_MAXITER):
            self.calc_fitness()
            self.sort_by_fitness()
            print("Best: " + str(self.population[0].fitness) + " and " + " d1= " + str(self.population[0].d1) + " d2= " + str(self.population[0].d2))
            self.mate()
            self.swap()

    def init_population(self):
        for i in range(0, GA_POPSIZE):
            gastruct = GAStruct(0.0, 0.0, 0.0)
            gastruct.d1 = random.uniform(0.2, 10)
            #gastruct.d2 = random.uniform(0.23, 0.27)
            gastruct.d2 = gastruct.d1
            self.population.append(gastruct)

    def calc_fitness(self):
        from geoexample import Geo
        for poulationobj in self.population:
            geo = Geo(poulationobj.d1, poulationobj.d2, 0.0000001)
            poulationobj.fitness = geo.getfitness()

    def sort_by_fitness(self):
        self.population = sorted(self.population, key=lambda pop: pop.fitness)

    def mate(self):
        self.buffer = self.population
        esize = GA_POPSIZE * GA_ELITRATE
        for i in range(int(esize), GA_POPSIZE):
            i1 = random.randint(0, GA_POPSIZE / 2)
            i2 = random.randint(0, GA_POPSIZE / 2)

            self.buffer[i].d1 = self.population[i1].d1
            self.buffer[i].d2 = self.buffer[i].d1
            if random.uniform(0, 1) < GA_MUTATIONRATE:
                self.mutate(self.buffer[i])

    def mutate(self, member):
        #randomposition = random.randint(0, 1)
        #if randomposition == 0:
        #     member.d1 = random.uniform(0.23, 0.27)
        # else:
        #     member.d2 = random.uniform(0.23, 0.27)
        member.d1 = random.uniform(0.2, 10)
        member.d2 = member.d1

    def swap(self):
        self.population = self.buffer




if __name__ == '__main__':
    print("hello world")
    genetic = GeneticAlgorithm()
    #genetic.runsimulation()
    import geoexample
    #geoexample.get1dplot()
    geoexample.get2dplot()