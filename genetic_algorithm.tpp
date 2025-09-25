#ifndef GENETIC_ALGORITHM_TPP
#define GENETIC_ALGORITHM_TPP

#include "genetic_algorithm.h"

namespace genetic_algorithms {

// individual struct implementations
template<typename T>
individual<T>::individual(const T& chromosome, double fit) 
    : chromosome(chromosome), fitness(fit) {}

template<typename T>
individual<T>::individual() {}

template<typename T>
bool individual<T>::operator<(const individual &other) const {
    return fitness > other.fitness;
}

// genetic_algorithm class implementations
template <typename T>
genetic_algorithm<T>::genetic_algorithm(
    size_t pop_size, double cross_rate, double mut_rate,
    size_t elitism, size_t max_gen,
    std::function<double(const T&)> fit_func,
    std::function<T(const T&, const T&)> cross_func,
    std::function<void(T&)> mut_func,
    std::function<T()> init_func,
    selection_method __method)
    : population_size(pop_size), crossover_rate(cross_rate),
      mutation_rate(mut_rate), elitism_count(elitism),
      max_generations(max_gen), fitness_function(fit_func),
      crossover_function(cross_func), mutation_function(mut_func),
      initialization_function(init_func), uniform_dist(0.0, 1.0), method(__method) {
    
    std::random_device rd;
    rng.seed(rd());

    if (population_size == 0) {
        throw std::invalid_argument("population size must be greater than 0");
    }
}

template <typename T>
void genetic_algorithm<T>::initialize_population() {
    population.clear();
    for (size_t i = 0; i < population_size; ++i) {
        T chromosome = initialization_function();
        double fitness = fitness_function(chromosome);
        population.emplace_back(chromosome, fitness);
    }
    std::sort(population.begin(), population.end());
}

template <typename T>
individual<T> genetic_algorithm<T>::tournament_selection(size_t tournament_size) {
    if (tournament_size < 1) tournament_size = 1;
    if (tournament_size > population_size) tournament_size = population_size;

    std::vector<individual<T>> tournament;
    std::uniform_int_distribution<size_t> dist(0, population_size-1);

    for (size_t i = 0; i < tournament_size; ++i) {
        size_t random_index = dist(rng);
        tournament.push_back(population[random_index]);
    }

    return *std::max_element(tournament.begin(), tournament.end(),
        [](const individual<T>& a, const individual<T>& b) {
            return a.fitness < b.fitness;
        });
}

template <typename T>
std::vector<individual<T>> genetic_algorithm<T>::stochastic_universal_sampling() {
    std::vector<individual<T>> selected;

    // calculate total fitness
    double total_fitness = 0.0;
    for (const auto& indv : population) {
        total_fitness += indv.fitness;
    }

    if (total_fitness <= 0.0) {
        std::uniform_int_distribution<size_t> dist(0, population_size-1);
        for (size_t i = 0; i < population_size; ++i) {
            selected.push_back(population[dist(rng)]);
        }
        return selected;
    }

    // calculate selection probabilities
    std::vector<double> cumulative_probabilities;
    double cumulative = 0.0;

    for (const auto& indv : population) {
        double probability = indv.fitness / total_fitness;
        cumulative += probability;
        cumulative_probabilities.push_back(cumulative);
    }

    // SUS selection
    double pointer_distance = 1.0 / population_size;
    double start_point = uniform_dist(rng) * pointer_distance;

    for (size_t i = 0; i < population_size; ++i) {
        double pointer = start_point + i * pointer_distance;

        for (size_t j = 0; j < population_size; ++j) {
            if (pointer <= cumulative_probabilities[j]) {
                selected.push_back(population[j]);
                break;
            }
        }
    }
    return selected;
}

template <typename T>
individual<T> genetic_algorithm<T>::roulette_wheel_selection() {
    double total_fitness = 0.0;
    for (const auto &indv : population) {
        total_fitness += indv.fitness;
    }
    if (total_fitness <= 0.0) {
        std::uniform_int_distribution<size_t> dist(0, population_size-1);
        return population[dist(rng)];
    }
    double r = uniform_dist(rng) * total_fitness;
    double cumulative = 0.0;
    for (const auto& indv : population) {
        cumulative += indv.fitness;
        if (cumulative >= r) {
            return indv;
        }
    }
    return population.back();
}

template <typename T>
individual<T> genetic_algorithm<T>::select_parent() {
    switch(method) {
        case selection_method::SUS: {
            auto selected = stochastic_universal_sampling();
            std::uniform_int_distribution<size_t> dist(0, selected.size() - 1);
            return selected[dist(rng)];
        }
        case selection_method::ROULETTE_WHEEL: {
            return roulette_wheel_selection();
        }
        case selection_method::TOURNAMENT:
        default:
            return tournament_selection();
    }
}

template <typename T>
void genetic_algorithm<T>::evolve() {
    std::vector<individual<T>> new_population;

    // elitism
    for (size_t i = 0; i < elitism_count && i < population_size; ++i) {
        new_population.push_back(population[i]);
    }

    // If using SUS, pre-select the entire mating pool
    std::vector<individual<T>> mating_pool;
    if (method == selection_method::SUS) {
        mating_pool = stochastic_universal_sampling();
    }

    // create rest of the population
    while (new_population.size() < population_size) {
        individual<T> parent1, parent2;
        if (method == selection_method::SUS) {
            std::uniform_int_distribution<size_t> dist(0, mating_pool.size()-1);
            parent1 = mating_pool[dist(rng)];
            parent2 = mating_pool[dist(rng)];
        } else {
            parent1 = select_parent();
            parent2 = select_parent();
        }

        // crossover
        T offspring_chromosome;
        if (uniform_dist(rng) < crossover_rate) {
            offspring_chromosome = crossover_function(parent1.chromosome, parent2.chromosome);
        } else {
            // no crossover, inherit from one parent
            offspring_chromosome = (uniform_dist(rng) < 0.5) ? parent1.chromosome : parent2.chromosome;
        }

        // mutation
        if (uniform_dist(rng) < mutation_rate) {
            mutation_function(offspring_chromosome);
        }

        double offspring_fitness = fitness_function(offspring_chromosome);
        new_population.emplace_back(offspring_chromosome, offspring_fitness);
    }

    population = std::move(new_population);
    std::sort(population.begin(), population.end());
}

template <typename T>
individual<T> genetic_algorithm<T>::run(bool verbose) {
    initialize_population();
    if (verbose) {
        std::cout << "Generation 0 - Best fitness: " << population[0].fitness << std::endl;
    }

    for (size_t generation = 1; generation <= max_generations; ++generation) {
        evolve();
        if (verbose && (generation % 100 == 0 || generation == max_generations)) {
            std::cout << "Generation " << generation << " - Best fitness: " << population[0].fitness << std::endl;
        }
    }
    return population[0];
}

template <typename T>
const std::vector<individual<T>>& genetic_algorithm<T>::get_population() const {
    return population;
}

template <typename T>
individual<T> genetic_algorithm<T>::best_individual() const {
    return population[0];
}

template <typename T>
void genetic_algorithm<T>::set_selection_method(selection_method __method) {
    method = __method;
}

template <typename T>
typename genetic_algorithm<T>::selection_method genetic_algorithm<T>::get_selection_method() const {
    return method;
}

} // namespace genetic_algorithms

#endif // GENETIC_ALGORITHM_TPP
