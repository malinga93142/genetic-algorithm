#ifndef PROBLEMS_H
#define PROBLEMS_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>

namespace example_problems {

struct knapsack_item {
    double weight;
    double value;
};

using knapsack_solution = std::vector<bool>;

class knapsack_problem {
private:
    std::vector<knapsack_item> items;
    double max_capacity;

public:
    knapsack_problem(const std::vector<knapsack_item> &item_list, double capacity);
    double fitness_function(const knapsack_solution& solution) const;
    knapsack_solution crossover_function(const knapsack_solution& parent1, const knapsack_solution& parent2) const;
    void mutation_function(knapsack_solution& solution) const;
    knapsack_solution initialize_function() const;
};

struct city {
    double x;
    double y;
};

using tsp_solution = std::vector<size_t>;

class tsp_problem {
private:
    std::vector<city> cities;

public:
    tsp_problem(const std::vector<city>& _cities);
    double fitness_function(const tsp_solution& solution) const;
    tsp_solution crossover_function(const tsp_solution& p1, const tsp_solution& p2);
    void mutation_function(tsp_solution& solution) const;
    tsp_solution initialization_function() const;
};

} // namespace example_problems

#endif // PROBLEMS_H
