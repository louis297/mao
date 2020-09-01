import random
from operator import itemgetter
from typing import List
from datetime import datetime
from math import exp, sqrt

class Gene:
    """
    This is a class to represent individual(Gene) in GA algorithom
    each object of this class have two attribute: data, size
    """

    def __init__(self, **data):
        self.__dict__.update(data)
        self.size = len(data['data'])  # length of gene


class GA:
    """
    This is a class of GA algorithm.
    """

    def __init__(self, parameter):
        """
        Initialize the pop of GA algorithom and evaluate the pop by computing its' fitness value.
        The data structure of pop is composed of several individuals which has the form like that:

        {'Gene':a object of class Gene, 'fitness': 1.02(for example)}
        Representation of Gene is a list: [b s0 u0 sita0 s1 u1 sita1 s2 u2 sita2]

        """
        # parameter = [CXPB, MUTPB, NGEN, POP_SIZE, [low, up], g1]

        self.CXPB = parameter[0]
        self.MUTPB = parameter[1]
        self.NGEN = parameter[2]
        self.POP_SIZE = parameter[3]
        self.bound = parameter[4]
        self.g1 = parameter[5]
        self.s_inds = ''

        # TODO: Read pop data from existing file
        # filename pattern: ga_mao_1_<CXPB>_<MUTPB>_<NGEN>_<POP_SIZE>.txt
        pop = []
        for i in range(self.POP_SIZE):
            geneinfo = []
            for pos in range(len(self.bound[0])):
                geneinfo.append(random.uniform(self.bound[0][pos], self.bound[1][pos]))  # initialise popluation

            fitness = self.evaluate(geneinfo)  # evaluate each chromosome
            pop.append({'Gene': Gene(data=geneinfo), 'fitness': fitness})  # store the chromosome and its fitness

        self.pop = pop
        self.bestindividual = self.selectBest(self.pop)  # store the best chromosome in the population

    def evaluate(self, geneinfo):
        """
        fitness function
        """
        # y = sum(x**(2) for i, x in enumerate(geneinfo)) / (geneinfo[2] ** 3 + geneinfo[3] ** 4)
        # para = [a, b, Xi1, Zeta1, Lambda1, Xi2, Zeta2, Lambda2, g1]
        # energy = self.calculate_from_mao(geneinfo)
        energy = self.calculate_from_louis(geneinfo)
        # return y for maximum fitness, 1/y for minimum fitness
        # the fitness relates to the chance of selection, do not simply use something like 10000000000-y
        return 2-energy

    def calculate_from_mao(self, geneinfo):
        a, Xi1, Zeta1, Lambda1, Xi2, Zeta2, Lambda2 = geneinfo
        omega = 1
        t = 1.1
        # the range of g1 is from 0 to 4.
        g = (sqrt(2) * self.g1) / omega

        b = 1 / 2 * (-((4 * a * exp(
            -((g ** 2 * Zeta2 * Lambda1 ** 2 * Xi1) / (2 * (Zeta2 + Xi1))) - (g ** 2 * Zeta1 * Lambda2 ** 2 * Xi2) / (
                    2 * (Zeta1 + Xi2))) * Zeta1 ** (1 / 4) * Zeta2 ** (1 / 4) * Xi1 ** (1 / 4) * Xi2 ** (1 / 4)) / (
                               sqrt(Zeta2 + Xi1) * sqrt(Zeta1 + Xi2))) + sqrt(2) * sqrt(1 - 2 * a ** 2 + (
                8 * a ** 2 * exp(
            -((g ** 2 * Zeta2 * Lambda1 ** 2 * Xi1) / (Zeta2 + Xi1)) - (g ** 2 * Zeta1 * Lambda2 ** 2 * Xi2) / (
                    Zeta1 + Xi2)) * sqrt(Zeta1) * sqrt(Zeta2) * sqrt(Xi1) * sqrt(Xi2)) / ((Zeta2 + Xi1) * (
                Zeta1 + Xi2))))

        result = (a ** 2 * (Zeta1 + Xi1 + Zeta1 * Xi1 * (Zeta1 + 2 * g ** 2 * (-2 + Lambda1) * Lambda1 + Xi1))) / (
                2 * Zeta1 * Xi1) \
                 + (b ** 2 * (Zeta2 + Xi2 + Zeta2 * Xi2 * (Zeta2 + 2 * g ** 2 * Lambda2 ** 2 + Xi2))) / (
                         2 * Zeta2 * Xi2) \
                 + 2 * t * ((2 * a ** 2 * exp(-((g ** 2 * Zeta1 * Lambda1 ** 2 * Xi1) / (Zeta1 + Xi1))) * sqrt(
            Zeta1) * sqrt(Xi1)) / (Zeta1 + Xi1) + (2 * b ** 2 * exp(
            -((g ** 2 * Zeta2 * Lambda2 ** 2 * Xi2) / (Zeta2 + Xi2))) * sqrt(Zeta2) * sqrt(Xi2)) / (Zeta2 + Xi2) + (
                                    2 * a * b * exp(
                                -((g ** 2 * (Lambda1 - Lambda2) ** 2 * Xi1 * Xi2) / (2 * (Xi1 + Xi2)))) * Zeta1 ** (
                                            1 / 4) * Zeta2 ** (1 / 4) * Xi1 ** (1 / 4) * Xi2 ** (1 / 4)) / (
                                    sqrt(Zeta1 + Zeta2) * sqrt(Xi1 + Xi2)) + (2 * a * b * exp(
            -((g ** 2 * (-Lambda1 + Lambda2) ** 2 * Xi1 * Xi2) / (2 * (Xi1 + Xi2)))) * Zeta1 ** (1 / 4) * Zeta2 ** (
                                                                                      1 / 4) * Xi1 ** (
                                                                                      1 / 4) * Xi2 ** (
                                                                                      1 / 4)) / (
                                    sqrt(Zeta1 + Zeta2) * sqrt(Xi1 + Xi2))) \
                 - 1 / ((Zeta2 + Xi1) ** (5 / 2) * (Zeta1 + Xi2) ** (5 / 2)) * 3 * a * b * exp(
            -((g ** 2 * Zeta2 * Lambda1 ** 2 * Xi1) / (2 * (Zeta2 + Xi1))) - (g ** 2 * Zeta1 * Lambda2 ** 2 * Xi2) / (
                    2 * (Zeta1 + Xi2))) * Zeta1 ** (1 / 4) * Zeta2 ** (1 / 4) * Xi1 ** (1 / 4) * Xi2 ** (1 / 4) * (
                         Xi2 * (-Xi1 * (Xi1 + Xi2 + g ** 2 * (
                         -2 + Lambda1) * Lambda1 * Xi1 * Xi2 + g ** 2 * Lambda2 ** 2 * Xi1 * Xi2) - Zeta2 ** 2 * (
                                        1 + g ** 2 * Lambda2 ** 2 * Xi2 + Xi1 * Xi2 - g ** 2 * Lambda1 ** 2 * Xi1 ** 2 * Xi2) - Zeta2 * (
                                        Xi2 + Xi1 ** 2 * Xi2 + Xi1 * (
                                        2 - 2 * g ** 2 * Lambda1 * Xi2 + 2 * g ** 2 * Lambda2 ** 2 * Xi2))) + Zeta1 * (
                                 Zeta2 ** 2 * (
                                 -1 - 2 * Xi1 * Xi2 + 2 * g ** 2 * Lambda1 ** 2 * Xi1 ** 2 * Xi2 - Xi2 ** 2) - 2 * Zeta2 * (
                                         Xi2 + Xi1 ** 2 * Xi2 + Xi1 * (
                                         1 - 2 * g ** 2 * Lambda1 * Xi2 + Xi2 ** 2)) - Xi1 * (
                                         2 * Xi2 + Xi1 * (1 + 2 * g ** 2 * (
                                         -2 + Lambda1) * Lambda1 * Xi2 + Xi2 ** 2))) + Zeta1 ** 2 * (
                                 Zeta2 ** 2 * (-Xi1 + g ** 2 * Lambda1 ** 2 * Xi1 ** 2 + Xi2 * (
                                 -1 + g ** 2 * Lambda2 ** 2 * Xi2)) + Zeta2 * (
                                         -1 + 2 * g ** 2 * Lambda1 * Xi1 - Xi1 ** 2 + 2 * Xi1 * Xi2 * (
                                         -1 + g ** 2 * Lambda2 ** 2 * Xi2)) - Xi1 * (
                                         1 + g ** 2 * (-2 + Lambda1) * Lambda1 * Xi1 + Xi1 * (
                                         Xi2 - g ** 2 * Lambda2 ** 2 * Xi2 ** 2)))) \
                 + 1 / ((Zeta2 + Xi1) ** (5 / 2) * (Zeta1 + Xi2) ** (5 / 2)) * a * b * exp(
            -((g ** 2 * Zeta2 * Lambda1 ** 2 * Xi1) / (2 * (Zeta2 + Xi1))) - (g ** 2 * Zeta1 * Lambda2 ** 2 * Xi2) / (
                    2 * (Zeta1 + Xi2))) * Zeta1 ** (1 / 4) * Zeta2 ** (1 / 4) * Xi1 ** (1 / 4) * Xi2 ** (1 / 4) * (
                         Xi2 * (Xi1 * (Xi1 + Xi2 + g ** 2 * (
                         -2 + Lambda1) * Lambda1 * Xi1 * Xi2 + g ** 2 * Lambda2 ** 2 * Xi1 * Xi2) + Zeta2 ** 2 * (
                                        1 + g ** 2 * Lambda2 ** 2 * Xi2 + Xi1 * Xi2 - g ** 2 * Lambda1 ** 2 * Xi1 ** 2 * Xi2) + Zeta2 * (
                                        Xi2 + Xi1 ** 2 * Xi2 + Xi1 * (
                                        2 - 2 * g ** 2 * Lambda1 * Xi2 + 2 * g ** 2 * Lambda2 ** 2 * Xi2))) + Zeta1 * (
                                 Zeta2 ** 2 * (
                                 1 + 2 * Xi1 * Xi2 - 2 * g ** 2 * Lambda1 ** 2 * Xi1 ** 2 * Xi2 + Xi2 ** 2) + 2 * Zeta2 * (
                                         Xi2 + Xi1 ** 2 * Xi2 + Xi1 * (
                                         1 - 2 * g ** 2 * Lambda1 * Xi2 + Xi2 ** 2)) + Xi1 * (
                                         2 * Xi2 + Xi1 * (1 + 2 * g ** 2 * (
                                         -2 + Lambda1) * Lambda1 * Xi2 + Xi2 ** 2))) + Zeta1 ** 2 * (
                                 Zeta2 ** 2 * (
                                 Xi1 - g ** 2 * Lambda1 ** 2 * Xi1 ** 2 + Xi2 - g ** 2 * Lambda2 ** 2 * Xi2 ** 2) + Zeta2 * (
                                         1 - 2 * g ** 2 * Lambda1 * Xi1 + Xi1 ** 2 - 2 * Xi1 * Xi2 * (
                                         -1 + g ** 2 * Lambda2 ** 2 * Xi2)) + Xi1 * (
                                         1 + g ** 2 * (-2 + Lambda1) * Lambda1 * Xi1 + Xi1 * (
                                         Xi2 - g ** 2 * Lambda2 ** 2 * Xi2 ** 2)))) \
                 - 1

        energy = omega * result


        return energy

    def calculate_from_louis(self, geneinfo):
        a, xi1, zeta1, lambda1, xi2, zeta2, lambda2 = geneinfo
        omega = 1
        t = 1.1
        # the range of g1 is from 0 to 4.
        g = (sqrt(2) * self.g1) / omega

        xi12zeta12 = xi1 * xi2 * zeta1 * zeta2
        xi12zeta12qsqrt = xi12zeta12 ** 0.25

        b_1_1_1 = g * g * zeta2 * lambda1 * lambda1 * xi1
        b_1_1_2 = zeta2 + xi1
        b_1_2_1 = g * g * zeta1 * lambda2 * lambda2 * xi2
        b_1_2_2 = zeta1 + xi2
        b_1_1 = b_1_1_1 / 2 / b_1_1_2
        b_1_2 = b_1_2_1 / 2 / b_1_2_2
        b_1 = 4 * a * exp(-b_1_1 - b_1_2) * xi12zeta12qsqrt / (sqrt(b_1_1_2 * b_1_2_2))
        b_2_1 = b_1_1_1 / b_1_1_2
        b_2_2 = b_1_2_1 / b_1_2_2
        b_2 = 1 - 2 * a * a + (8 * a * a * exp(-b_2_1 - b_2_2) * sqrt(xi12zeta12)) / (b_1_1_2 * b_1_2_2)

        b = 0.5 * (-b_1 + sqrt(2 * b_2))

        e_1 = a * a * (zeta1 + xi1 + zeta1*xi1*(zeta1 + 2*g*g * (-2+lambda1) * lambda1 + xi1)) / 2 / zeta1 / xi1
        e_2 = b * b * (zeta2 + xi2 + zeta2*xi2*(zeta2 + 2*g*g * lambda2 * lambda2 + xi2)) / 2 / zeta2 / xi2
        e_t1 = 2*a*a * exp(-g*g*zeta1*lambda1*lambda1*xi1/(zeta1+xi1)) * sqrt(zeta1*xi1) / (zeta1+xi1)
        e_t2 = 2*b*b * exp(-g*g*zeta2*lambda2*lambda2*xi2/(zeta2+xi2)) * sqrt(zeta2*xi2) / (zeta2+xi2)
        e_t3 = 2*a*b * exp(-(g*(lambda1-lambda2))**2 * xi1*xi2 /2/(xi1+xi2)) * xi12zeta12qsqrt / sqrt((zeta1+zeta2) * (xi1+xi2))
        e_t4 = 2*a*b * exp(-(g*(-lambda1+lambda2))**2 *xi1*xi2 /2/(xi1+xi2)) * xi12zeta12qsqrt / sqrt((zeta1+zeta2) * (xi1+xi2))
        e_t = 2*t*(e_t1 + e_t2 + e_t3 + e_t4)

        e_3_1 = 1 / ((zeta2+xi1) * (zeta1+xi2)) ** 2.5 * a*b * exp(-b_1_1 - b_1_2) * xi12zeta12qsqrt
        e_3_2_1 = xi1 + xi2 + g*g*(-2+lambda1)*lambda1*xi1*xi2 + g*g*lambda2*lambda2*xi1*xi2
        e_3_2_2 = 1 + g*g*lambda2*lambda2*xi2 + xi1*xi2 - (g*lambda1*xi1)**2 * xi2
        e_3_2_3 = xi2 + xi1*xi1*xi2 + xi1*2*(1 - g*g*lambda1*xi2 + g*g*lambda2*lambda2*xi2)
        e_3_2_a = xi2 * (-xi1*e_3_2_1 - zeta2*zeta2*e_3_2_2 - zeta2*e_3_2_3)

        e_3_2_4 = -1-2*xi1*xi2 + 2*(g*lambda1*xi1)**2 * xi2 - xi2*xi2
        e_3_2_5 = xi2 + xi1*xi1*xi2 + xi1 * (1 - 2*g*g*lambda1*xi2 + xi2*xi2)
        e_3_2_6 = 2*xi2 + xi1 * (1 + 2*g*g*(-2+lambda1)*lambda1*xi2 + xi2*xi2)
        e_3_2_b = zeta1 * (zeta2*zeta2*e_3_2_4 - 2*zeta2*e_3_2_5 - xi1*e_3_2_6)

        e_3_2_7 = -xi1 + (g*lambda1*xi1) ** 2 + xi2 * (-1+g*g*lambda2*lambda2*xi2)
        e_3_2_8 = -1 + 2*g*g*lambda1*xi1 - xi1*xi1 + 2*xi1*xi2 * (-1 + g*g*lambda2*lambda2*xi2)
        e_3_2_9 = 1 + g*g*(-2+lambda1)*lambda1*xi1 + xi1*(xi2 - (g*lambda2*xi2)**2)
        e_3_2_c = zeta1*zeta1 * (zeta2*zeta2*e_3_2_7 + zeta2*e_3_2_8 - xi1*e_3_2_9)
        e_3_2 = e_3_2_a + e_3_2_b + e_3_2_c
        e_3 = 4 * e_3_1 * e_3_2

        result = e_1 + e_2 + e_t - e_3 - 1
        energy = omega * result

        return energy

    def selectBest(self, pop):
        """
        select the best individual from pop
        """
        self.s_inds = sorted(pop, key=itemgetter("fitness"), reverse=True)  # from large to small, return a pop

        return self.s_inds[0]

    def selection(self):
        """
        select some good individuals from pop, note that good individuals have greater probability to be chosen
        for example: a fitness list like that:[5, 4, 3, 2, 1], sum is 15,
        [-----|----|---|--|-]
        012345|6789|101112|1314|15
        we randomly choose a value in [0, 15],
        it belongs to first scale with greatest probability
        """
        self.s_inds = sorted(self.pop, key=itemgetter("fitness"), reverse=True)  # sort the pop by the reference of fitness
        sum_fits = sum(ind['fitness'] for ind in self.pop)  # sum up the fitness of the whole pop

        chosen = []
        for i in range(self.POP_SIZE):
            u = random.random() * sum_fits  # randomly produce a num in the range of [0, sum_fits], as threshold
            sum_ = 0
            for ind in self.s_inds:
                sum_ += ind['fitness']  # sum up the fitness
                if sum_ >= u:
                    # when the sum of fitness is bigger than u, choose the one, which means u is in the range of
                    # [sum(1,2,...,n-1),sum(1,2,...,n)] and is time to choose the one ,namely n-th individual in the pop
                    chosen.append(ind)
                    break

        chosen = sorted(chosen, key=itemgetter("fitness"), reverse=False)
        return chosen

    def crossoperate(self, offspring):
        """
        cross operation
        here we use two points crossoperate
        for example: gene1: [5, 2, 4, 7], gene2: [3, 6, 9, 2], if pos1=1, pos2=2
        5 | 2 | 4  7
        3 | 6 | 9  2
        =
        3 | 2 | 9  2
        5 | 6 | 4  7
        """
        dim = len(offspring[0]['Gene'].data)

        geninfo1 = offspring[0]['Gene'].data  # Gene's data of first offspring chosen from the selected pop
        geninfo2 = offspring[1]['Gene'].data  # Gene's data of second offspring chosen from the selected pop

        if dim == 1:
            pos1 = 1
            pos2 = 1
        else:
            pos1 = random.randrange(1, dim)  # select a position in the range from 0 to dim-1,
            pos2 = random.randrange(1, dim)

        newoff1 = Gene(data=[])  # offspring1 produced by cross operation
        newoff2 = Gene(data=[])  # offspring2 produced by cross operation
        temp1 = []
        temp2 = []
        for i in range(dim):
            if min(pos1, pos2) <= i < max(pos1, pos2):
                temp2.append(geninfo2[i])
                temp1.append(geninfo1[i])
            else:
                temp2.append(geninfo1[i])
                temp1.append(geninfo2[i])
        newoff1.data = temp1
        newoff2.data = temp2

        return newoff1, newoff2

    def mutation(self, crossoff):
        """
        mutation operation
        """
        dim = len(crossoff.data)

        if dim == 1:
            pos = 0
        else:
            pos = random.randrange(0, dim)  # chose a position in crossoff to perform mutation.

        crossoff.data[pos] = random.uniform(self.bound[0][pos], self.bound[1][pos])
        return crossoff

    def GA_main(self):
        """
        main frame work of GA
        """
        print("Start of evolution")

        # Begin the evolution
        for g in range(self.NGEN+1):
            # Apply selection based on their converted fitness
            selectpop = self.selection()

            next_offspring = []
            while len(next_offspring) != self.POP_SIZE:
                # Apply crossover and mutation on the offspring

                # Select two individuals
                offspring = [selectpop.pop() for _ in range(2)]

                if random.random() < self.CXPB:  # cross two individuals with probability CXPB
                    cross_offspring1, cross_offspring2 = self.crossoperate(offspring)
                    if random.random() < self.MUTPB:  # mutate an individual with probability MUTPB
                        cross_offspring1 = self.mutation(cross_offspring1)
                        cross_offspring2 = self.mutation(cross_offspring2)
                    fit_cross_offspring1 = self.evaluate(cross_offspring1.data)  # Evaluate the individuals
                    fit_cross_offspring2 = self.evaluate(cross_offspring2.data)
                    next_offspring.append({'Gene': cross_offspring1, 'fitness': fit_cross_offspring1})
                    next_offspring.append({'Gene': cross_offspring2, 'fitness': fit_cross_offspring2})
                else:
                    next_offspring.extend(offspring)

            # The population is entirely replaced by the offspring
            self.pop = next_offspring

            # Gather all the fitnesses in one list and print the stats
            fits = [ind['fitness'] for ind in self.pop]

            best_ind = self.selectBest(self.pop)

            if best_ind['fitness'] > self.bestindividual['fitness']:
                self.bestindividual = best_ind

            if g % 500 == 0:
                print("############### Generation {} : {} ###############".format(g, self.g1))

                print("Best individual found is {}, {}".format(self.bestindividual['Gene'].data,
                                                               2-self.bestindividual['fitness']))
                print("  Max fitness of current pop: {}".format(2-max(fits)))

        print("------ End of (successful) evolution ------")
        return 2-self.bestindividual['fitness']