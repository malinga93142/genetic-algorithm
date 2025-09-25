#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H
#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace genetic_algorithms{
  template <typename T>
    struct individual{
      T chromosome;
      double fitness;
      individual(const T& chromosome, double fit);
      individual();
      bool operator <(const individual &other) const;
    };

  template <typename T>
  class genetic_algorithm{
    public:
      enum class selection_method{
        TOURNAMENT,
        SUS,
        ROULETTE_WHEEL
      };
      genetic_algorithm(size_t pop_size, double cross_rate, double mut_rate,
          size_t elitism, size_t max_gen,
          std::function<double(const T&)> fit_func,
          std::function<T(const T&, const T&)> cross_func,
          std::function<void(T&)> mut_func,
          std::function<T()> init_func,
          selection_method __method = selection_method::TOURNAMENT);

      void initialize_population();
      individual<T> tournament_selection(size_t tournament_size = 3);
      std::vector<individual<T>> stochastic_universal_sampling();
      individual<T> roulette_wheel_selection();
      individual<T> select_parent();
      void evolve();
      individual<T> run(bool verbose = false);

      // Getters and Setters
      const std::vector<individual<T>> &get_population() const;
      individual<T> best_individual() const;
      void set_selection_method(selection_method __method);
      selection_method get_selection_method() const;
    private:
      size_t population_size;
      double crossover_rate;
      double mutation_rate;
      size_t elitism_count;
      size_t max_generations;
      selection_method method;


      // function pointers for problem-specific operations
      std::function<double(T&)> fitness_function;
      std::function<T(const T&, const T&)> crossover_function;
      std::function<void(T&)> mutation_function;
      std::function<T()> initialization_function;

      // population
      std::vector<individual<T>> population;
      
      // random number generation
      std::mt19937 rng;
      std::uniform_real_distribution<double> uniform_dist;

  };
} // namespace genetic_algorithms
#include "genetic_algorithm.tpp"
#endif  // GENETIC_ALGORITHM_H
