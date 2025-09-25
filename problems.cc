#include "problems.h"
#include <random>
#include <algorithm>
#include <numeric>

namespace example_problems {

// knapsack_problem implementations
knapsack_problem::knapsack_problem(const std::vector<knapsack_item> &item_list, double capacity)
    : items(item_list), max_capacity(capacity) {}

double knapsack_problem::fitness_function(const knapsack_solution& solution) const {
    double total_value = 0.0;
    double total_weight = 0.0;

    for (size_t i = 0; i < solution.size(); ++i) {
        if (solution[i]) {
            total_value += items[i].value;
            total_weight += items[i].weight;
        }
    }
    // penalize solutions that exceed capacity
    if (total_weight > max_capacity) {
        return 0.0;
    }
    return total_value;
}

knapsack_solution knapsack_problem::crossover_function(const knapsack_solution& parent1, const knapsack_solution& parent2) const {
    knapsack_solution offspring = parent1;

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, parent1.size()-1);

    // single-point crossover
    size_t crossover_point = dist(rng);
    for (size_t i = crossover_point; i < offspring.size(); ++i) {
        offspring[i] = parent2[i];
    }
    return offspring;
}

void knapsack_problem::mutation_function(knapsack_solution& solution) const {
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> dist(0, solution.size()-1);

    // flip random bit
    size_t mutation_point = dist(rng);
    solution[mutation_point] = !solution[mutation_point];
}

knapsack_solution knapsack_problem::initialize_function() const {
    knapsack_solution solution(items.size());
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, 1);

    for (size_t i = 0; i < items.size(); ++i) {
        solution[i] = (dist(rng) == 1);
    }
    return solution;
}

// tsp_problem implementations
tsp_problem::tsp_problem(const std::vector<city>& _cities) : cities(_cities) {}

double tsp_problem::fitness_function(const tsp_solution& solution) const {
    double total_dist = 0.0;
    for (size_t i = 0; i < solution.size(); ++i) {
        const city& from = cities[solution[i]];
        const city& to = cities[solution[(i+1) % solution.size()]];
        double dx = from.x - to.x;
        double dy = from.y - to.y;
        total_dist += std::sqrt(dx*dx + dy*dy);
    }
    return 1.0 / total_dist;
}

tsp_solution tsp_problem::crossover_function(const tsp_solution& p1, const tsp_solution& p2) {
    size_t n = p1.size();
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<size_t> dist(0, n-1);
    
    size_t start = dist(rng), end = dist(rng);
    if (start > end) std::swap(start, end);
    
    tsp_solution offspring(n, n); // Initialize with value 'n' (invalid index)
    std::vector<bool> in_offspring(n, false);
    
    // copy substring from p1
    for (size_t i = start; i <= end; ++i) {
        offspring[i] = p1[i];
        in_offspring[p1[i]] = true;
    }
    
    // fill remaining from p2
    size_t curr = (end + 1) % n;
    for (size_t i = 0; i < n; ++i) {
        size_t index = (end + 1 + i) % n;
        if (!in_offspring[p2[index]]) {
            offspring[curr] = p2[index];
            curr = (curr + 1) % n;
        }
    }
    return offspring;
}

void tsp_problem::mutation_function(tsp_solution& solution) const {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<size_t> dist(0, solution.size()-1);
    size_t i = dist(rng), j = dist(rng);
    std::swap(solution[i], solution[j]);
}

tsp_solution tsp_problem::initialization_function() const {
    tsp_solution sol(cities.size());
    std::iota(sol.begin(), sol.end(), 0);
    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(sol.begin(), sol.end(), rng);
    return sol;
}

} // namespace example_problems
