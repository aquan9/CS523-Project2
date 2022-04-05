#Author: Andres Quan
#Brief: A program to simulate B-cell mutation using a genetic algorithm tuned to similar properties

import pygad
import random



#Generate the epitope we're trying to match to

#uncomment for debug
correct_epitope = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,]

#uncomment for random value
#correct_epitope = []
#for i in range(0, 128):
#    correct_epitope.append(random.randint(0,1))

print("Correct Epitope: " + str(correct_epitope))

def fitness_function(solution, solution_idx):
    #print("Fitness function: " + str(solution))
    #print("Sum solution: " + str(sum(solution)))
    fitness = 0
    for i in range(0, 128):
        if solution[i] == correct_epitope[i]:
            fitness += 1
    return fitness

#Using binary values by using:
#gene_type=int
#int_range_low=0
#int_range_high=0

#Getting 2^128 possible Epitopes by setting:
#num_genes=128

#B cells differentially survive based on fitness achived by setting:
#parent_selection_type="rank"

#Only 3 B cell types per germinal center at a time 

ga_instance = pygad.GA(num_generations=20000,
                       num_parents_mating=2,
                       sol_per_pop=9,
                       num_genes=128,
                       keep_parents=1,
                       parent_selection_type="rank",
                       fitness_func=fitness_function,
                       #mutation_by_replacement=True,
                       mutation_type = "inversion",
                       mutation_percent_genes = 1,
                       crossover_type = "scattered",
                       init_range_low=0,
                       init_range_high=2,
                       stop_criteria="reach_128",
                       gene_type=int,
                        )

print(ga_instance.initial_population)
print(ga_instance.initial_population.shape)



ga_instance.run()

solution, solution_fitness, solution_idx = ga_instance.best_solution()
#print("Completed generations: " + str(ga_instance.generations_completed))
ga_instance.plot_fitness()


time_to_complete_128 = []
for i in range(0,100):
    ga_instance2 =  pygad.GA(num_generations=20000,
                       num_parents_mating=2,
                       sol_per_pop=9,
                       num_genes=128,
                       keep_parents=1,
                       parent_selection_type="rank",
                       fitness_func=fitness_function,
                       #mutation_by_replacement=True,
                       mutation_type = "inversion",
                       mutation_percent_genes = 1,
                       crossover_type = "scattered",
                       init_range_low=0,
                       init_range_high=2,
                       stop_criteria="reach_128",
                       gene_type=int,
                        )

    ga_instance2.run()
    print("Completed generations: " + str(ga_instance2.generations_completed))
    time_to_complete_128.append(int(ga_instance2.generations_completed))

#90% matching
time_to_complete_115 = []
for i in range(0,100):
    ga_instance2 =  pygad.GA(num_generations=20000,
                       num_parents_mating=2,
                       sol_per_pop=9,
                       num_genes=128,
                       keep_parents=1,
                       parent_selection_type="rank",
                       fitness_func=fitness_function,
                       #mutation_by_replacement=True,
                       mutation_type = "inversion",
                       mutation_percent_genes = 1,
                       crossover_type = "scattered",
                       init_range_low=0,
                       init_range_high=2,
                       stop_criteria="reach_115",
                       gene_type=int,
                        )

    ga_instance2.run()
    print("Completed generations: " + str(ga_instance2.generations_completed))
    time_to_complete_115.append(int(ga_instance2.generations_completed))



#80% matching
time_to_complete_102 = []
for i in range(0,100):
    ga_instance2 =  pygad.GA(num_generations=20000,
                       num_parents_mating=2,
                       sol_per_pop=9,
                       num_genes=128,
                       keep_parents=1,
                       parent_selection_type="rank",
                       fitness_func=fitness_function,
                       #mutation_by_replacement=True,
                       mutation_type = "inversion",
                       mutation_percent_genes = 1,
                       crossover_type = "scattered",
                       init_range_low=0,
                       init_range_high=2,
                       stop_criteria="reach_102",
                       gene_type=int,
                        )

    ga_instance2.run()
    print("Completed generations: " + str(ga_instance2.generations_completed))
    time_to_complete_102.append(int(ga_instance2.generations_completed))

print("Total generations completed in each instance:")
print("128 genes matching (100%): ")
print(time_to_complete_128)
print("115 genes matching (90%): ")
print(time_to_complete_115)
print("102 genes matching (80%): ")
print(time_to_complete_102)

#Get the average
average_time_to_complete_128 = sum(time_to_complete_128) / 100
average_time_to_complete_115 = sum(time_to_complete_115) / 100
average_time_to_complete_102 = sum(time_to_complete_102) / 100

print("Average number of generations to match 128 genes: " + str(average_time_to_complete_128))
print("Average number of generations to match 115 genes: " + str(average_time_to_complete_115))
print("Average number of generations to match 102 genes: " + str(average_time_to_complete_102))


#---------------------------------------------------------------------

#Time to play the antigen elimination game!
#Assume gene matchup is directly correlated to probability of antibody neutralizing antigen
time_per_generation = 8 #hours
initial_antigen_population = 100000000
antibody_production_rate = 2000 * 3600 #per hour
initial_b_cells = 3

time_to_eliminate_128 = (average_time_to_complete_128 * time_per_generation) + ((initial_antigen_population / (antibody_production_rate * initial_b_cells)) * 128/128)
time_to_eliminate_115 = (average_time_to_complete_115 * time_per_generation) + ((initial_antigen_population / (antibody_production_rate * initial_b_cells)) * 115/128)
time_to_eliminate_102 = (average_time_to_complete_102 * time_per_generation) + ((initial_antigen_population / (antibody_production_rate * initial_b_cells)) * 102/128)

print("Time to eliminate all cells in hours (wait until 128 matching genes): " + str(time_to_eliminate_128))
print("Time to eliminate all cells in hours (wait until 115 matching genes): " + str(time_to_eliminate_115))
print("Time to eliminate all cells in hours (wait until 102 matching genes): " + str(time_to_eliminate_102))

