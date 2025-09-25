#include "genetic_algorithm.h"
#include "problems.h"
#include "genetic_algorithm.tpp" // Include template implementations
int main() {
    using namespace genetic_algorithms;
    using namespace example_problems;

    std::vector<knapsack_item> items = {
        {2.0, 10.0}, {3.0, 15.0}, {5.0, 25.0}, {7.0, 40.0}, {1.0, 5.0},
        {4.0, 20.0}, {6.0, 30.0}, {8.0, 45.0}, {9.0, 50.0}, {2.5, 12.0}
    };
    double max_cap = 20.0;
    
    knapsack_problem knapsack(items, max_cap);
    genetic_algorithm<knapsack_solution> ga(
        100,  // population size
        0.8,  // crossover rate
        0.1,  // mutation rate
        2,    // elitism count
        500,  // max generations
        [&knapsack](const knapsack_solution& sol) { 
            return knapsack.fitness_function(sol);  
        },
        [&knapsack](const knapsack_solution& p1, const knapsack_solution& p2) {
            return knapsack.crossover_function(p1, p2);
        },
        [&knapsack](knapsack_solution& sol) {
            knapsack.mutation_function(sol);
        },
        [&knapsack]() {
            return knapsack.initialize_function();
        },
        genetic_algorithm<knapsack_solution>::selection_method::TOURNAMENT
    );

    auto best_solution = ga.run(true);
    return 0;
}
