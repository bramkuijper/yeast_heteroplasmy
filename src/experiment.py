# mimicking the simple experiment by 
# Kazancioglu & Arnqvist
# assuming the action of drift only

import os, re, random, sys

import numpy as np

class Experiment:

    Females = []
    Males = []
    Offspring = []

    N = 500
    sr = 0.5
    mean_fec = 2.8 
    experiment_id = 0

    def __init__(self, experiment_id):

        self.experiment_id = experiment_id

        for i in range(0,self.N/2):

            # initialize the population
            self.Females.append(random.uniform(0,1) < 0.2)
            self.Males.append(random.uniform(0,1) < 0.2)

        # set number of males and females

        self.theGeneration = 1

        self.write_data()

        for i in range(2, 11):
            self.theGeneration = i
            self.generation()
            
            self.write_data()

        self.sample()

    # single experimental generation
    def generation(self):

        self.Offspring = []

        for female_i in self.Females:

            clutch = np.random.poisson(self.mean_fec)

            for j in range(0, clutch):

                # inherit the mitochondrion to the offspring
                self.Offspring.append(female_i)

        if len(self.Offspring) < self.N:
            print("too few offspring: Noff = ", str(len(self.Offspring)) + ", N: ", str(self.N))
            sys.exit(1)

        self.Females = []
        self.Males = []

        # sample N offspring out off all the offspring
        self.Offspring = random.sample(self.Offspring,self.N)

        # sample adults from offspring
        for i in range(0, self.N):

            if random.uniform(0,1) < 0.5:

                # make a female
                self.Females.append(self.Offspring[i])

            else:

                # make a male
                self.Males.append(self.Offspring[i])

        self.write_data()

    def write_data(self):

        # total number of type 1 mitochondria obtained by summing over all of them
        self.n1 = sum(self.Females + self.Males)

        self.n2 = len(self.Females) + len(self.Males) - self.n1

        proportion = float(self.n1) / (len(self.Females) + len(self.Females))

        # starting 0 means that this the actual average, not the sample
        print("0;" + str(self.experiment_id) + ";" + str(self.theGeneration) + ";" + str(self.n1) + ";" + str(self.n2) + ";" + str(proportion))

    # now sample mitochondria from the experiment
    def sample(self):

        total = self.Males + self.Females

        total_sample = sum(random.sample(total,10))

        proportion = float(total_sample) / 10

        # starting 1 means we're taking a sample
        print("1;" + str(self.experiment_id) + ";" + str(self.theGeneration) + ";" + str(total_sample) + ";" + str(10 -total_sample) + ";" + str(proportion))

nExperiment = 180

print("type;exp_id;generation;n1;n2;proportion")

for experiment_i in range(0,180):

    theExperiment = Experiment(experiment_i)

