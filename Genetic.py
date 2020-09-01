# Author: Louis
# Genetic Algorithms demo

import random
from typing import List
from datetime import datetime
from math import exp, sqrt
from gaclass import GA
import matplotlib.pyplot as plt

# from sympy import Symbol
# from sympy import exp,sqrt
# a=Symbol("a")
# b=Symbol("b")
# g=Symbol("g")
# g1=Symbol("g1")
# t=Symbol("t")
# omega=Symbol("omega")
# Xi1=Symbol("Xi1")
# Zeta1=Symbol("Zeta1")
# Lambda1=Symbol("Lambda1")
# Xi2=Symbol("Xi2")
# Zeta2=Symbol("Zeta2")
# Lambda2=Symbol("Lambda2")


##parameter range:
#0<a<sqrt(1/2)
#0<b<sqrt(1/2
#Xi1>0
#Xi2>0
#Zeta1>0
#Zeta2>0


def main():
    random.seed(datetime.now())
    # config
    CXPB, MUTPB, NGEN, POP_SIZE = 0.8, 0.1, 10000, 100  # POP_SIZE must be even number
    # up = [30, 30, 30, 30]  # upper range for variables
    # low = [1, 1, 1, 1]  # lower range for variables
    # para = [a, Xi1, Zeta1, Lambda1, Xi2, Zeta2, Lambda2]
    low = (0, 0, 0, -1, 0, 0, -1)
    high = (sqrt(1/2), 1, 1, 1, 1, 1, 1)
    bound = [low, high]
    # g1 = 0.8
    # bound = read_parameters()
    g_array = []
    e_array = []
    for i in range(30):
        g1 = 0.01 + 0.1 * i
        parameter = [CXPB, MUTPB, NGEN, POP_SIZE, bound, g1]
        run = GA(parameter)
        e = run.GA_main()
        g_array.append(g1)
        e_array.append(e)
    plt.plot(g_array, e_array)
    print(g_array, e_array)
    plt.show()

def read_parameters() -> (tuple, tuple):
    filename = 'ga1_const.txt'
    with open(filename, 'r') as infile:
        low: List[float] = []
        high: List[float] = []
        # \r\n for windows, \n for unix like OS, \r for old mac OS
        all_content = infile.read().replace('\r\n', '\n').replace('\r', '\n')
        content_list = all_content.strip().split('\n')
        for line in content_list:
            line = line.strip()
            # comment in file
            if line[0] == '#' or not line:
                continue
            line_content = line.split(',')
            low.append(float(line_content[1].strip()))
            high.append(float(line_content[2].strip()))
        return tuple(low), tuple(high)


if __name__ == "__main__":
    main()
