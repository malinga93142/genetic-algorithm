#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>


namespace genetic_algorithms{

template <typename T>
  class genetic_algorithm;

template<typename T>
  struct individual{
    T chromosome;
    double fitness;
    
    individual(const T& chromosome, double fit): chromosome(chromosome), fitness(fit) {}
    individual(){}
    bool operator<(const individual &other) const{
      return fitness > other.fitness;
    }
  };
template <typename T>
  class genetic_algorithm{
    public:
      enum class selection_method{
        TOURNAMENT,
        SUS, // Stochastic Universal Sampling
        ROULETTE_WHEEL
      };
    private:
      // configuration parameters
      size_t population_size;
      double crossover_rate;
      double mutation_rate;
      size_t elitism_count;
      size_t max_generations;
      selection_method method;

      // function pointers for problem-specific operations
      std::function<double(const T&)> fitness_function;
      std::function<T(const T&, const T&)> crossover_function;
      std::function<void(T&)> mutation_function;
      std::function<T()> initialization_function;
      
      // population
      std::vector<individual<T>> population;

      // Random number generation
      std::mt19937 rng;
      std::uniform_real_distribution<double> uniform_dist;
    public:
      genetic_algorithm(size_t pop_size, double cross_rate, double mut_rate,
          size_t elitism, size_t max_gen,
          std::function<double(const T&)> fit_func,
          std::function<T(const T&, const T&)> cross_func,
          std::function<void(T&)> mut_func,
          std::function<T()> init_func,
          selection_method __method = selection_method::TOURNAMENT)
        :
          population_size(pop_size), crossover_rate(cross_rate),
          mutation_rate(mut_rate), elitism_count(elitism),
          max_generations(max_gen), fitness_function(fit_func),
          crossover_function(cross_func), mutation_function(mut_func),
          initialization_function(init_func), uniform_dist(0.0, 1.0), method(__method)
    {
      std::random_device rd;
      rng.seed(rd());

      if (population_size == 0){
        throw std::invalid_argument("population size must be greater than 0");
      }
    }

      void initialize_population(){
        population.clear();
        for (size_t i = 0; i < population_size; ++i){
          T chromosome = initialization_function();
          double fitness = fitness_function(chromosome);
          population.emplace_back(chromosome, fitness);
        }
        std::sort(population.begin(), population.end());
      }

      individual<T> tournament_selection(size_t tournament_size = 3){
        if (tournament_size < 1) tournament_size = 1;
        if (tournament_size > population_size)  tournament_size = population_size;

        std::vector<individual<T>> tournament;
        std::uniform_int_distribution<size_t> dist(0, population_size-1);

        for (size_t i = 0; i < tournament_size; ++i){
          size_t random_index = dist(rng);
          tournament.push_back(population[random_index]);
        }

        return *std::max_element(tournament.begin(), tournament.end(),
            [](const individual<T>& a, const individual<T>& b){
              return a.fitness < b.fitness;
            });
      }

      // stochastic universal sampling (SUS)
      std::vector<individual<T>> stochastic_universal_sampling(){
        std::vector<individual<T>> selected;

        // calculate total fitness
        double total_fitness = 0.0f;
        for (const auto& indv: population){
          total_fitness += indv.fitness;
        }

        if (total_fitness <= 0.0){
          std::uniform_int_distribution<size_t> dist(0, population_size-1);
          for (size_t i = 0; i < population_size; ++i){
            selected.push_back(population[dist(rng)]);
          }
          return selected;
        }
        // calculate selection probabilities

        std::vector<double> probabilities;
        std::vector<double> cumulative_probabilities;
        double cumulative = 0.0;

        for (const auto& indv: population){
          double probability = indv.fitness / total_fitness;
          probabilities.push_back(probability);
          cumulative += probability;
          cumulative_probabilities.push_back(cumulative);
        }

        // SUS selection

        double pointer_distance = 1.0 / population_size;
        double start_point = uniform_dist(rng)*pointer_distance;

        for (size_t i = 0; i < population_size; ++i){
          double pointer = start_point + i * pointer_distance;

          for (size_t j = 0; j < population_size; ++j){
            if (pointer <= cumulative_probabilities[j]) {
              selected.push_back(population[j]);
              break;
            }
          }
        }
        return selected;
      }

      individual<T> roulette_wheel_selection(){
        double total_fitness = 0.0f;
        for (const auto &indv: population){
          total_fitness += indv.fitness;
        }
        if (total_fitness <= 0.0f){
          std::uniform_int_distribution<size_t> dist(0, population_size-1);
          return population[dist(rng)];
        }
        double r = uniform_dist(rng) * total_fitness;
        double cumulative = 0.0f;
        for (const auto& indv: population){
          cumulative += indv.fitness;
          if (cumulative >= r){
            return indv;
          }
        }
        return population.back();
      }
      individual<T> select_parent(){
        switch(method){
          case selection_method::SUS:
            {
              auto selected = stochastic_universal_sampling();
              std::uniform_int_distribution<size_t> dist(0, selected.size() - 1);
              return selected[dist(rng)];
            }
          case selection_method::ROULETTE_WHEEL:
            {
              return roulette_wheel_selection();
            }
          case selection_method::TOURNAMENT:
          default:
            return tournament_selection();
        }
      }
      void evolve(){
        std::vector<individual<T>> new_population;

        // elitism

        for (size_t i = 0; i < elitism_count && i < population_size; ++i){
          new_population.push_back(population[i]);
        }

        // If using SUS, pre-select the entire mating pool
        //
        std::vector<individual<T>> mating_pool;
        if (method == selection_method::SUS){
          mating_pool = stochastic_universal_sampling();
        }
        // create rest of the pop
        while (new_population.size() < population_size){
          individual<T> parent1, parent2;
          if (method == selection_method::SUS){
            std::uniform_int_distribution<size_t> dist(0, mating_pool.size()-1);
            parent1 = mating_pool[dist(rng)];
            parent2 = mating_pool[dist(rng)];
          }
          else{
            parent1 = select_parent();
            parent2 = select_parent();
          }
          // crossover
          T offspring_chromosome;
          if (uniform_dist(rng) < crossover_rate){
            offspring_chromosome = crossover_function(parent1.chromosome, parent2.chromosome);
          }
          else{
            // no crossover, inherit from one parent

            offspring_chromosome = (uniform_dist(rng) < 0.5)  ? parent1.chromosome  : parent2.chromosome;
          }

          // mutation

          if (uniform_dist(rng) < mutation_rate){
            mutation_function(offspring_chromosome);
          }
          double offspring_fitness = fitness_function(offspring_chromosome);
          new_population.emplace_back(offspring_chromosome, offspring_fitness);
        }
        population = std::move(new_population);
        std::sort(population.begin(), population.end());
      }

      
      individual<T> run(bool verbose=false){
        initialize_population();
        if (verbose){
          std::cout << "Generation 0 - Best fitness: " << population[0].fitness << std::endl;
        }

        for (size_t generation = 1; generation <= max_generations; ++generation){
          evolve();
          if (verbose && (generation % 100 == 0 || generation == max_generations)){
            std::cout << "Generation " << generation << " - Best fitness: " << population[0].fitness << std::endl;
          }
        }
        return population[0];
      }

      // get current population size
      const std::vector<individual<T>>& get_population() const{
        return population;
      }

      // get best individual
      individual<T> best_individual() const{
        return population[0];
      }
      void set_selection_method(selection_method __method){
        method = __method;
      }
      selection_method get_selection_method() const{
        return method;
      }
  };
} // namespace genetic_algorithm

 
namespace example_problems{
  struct knapsack_item{
    double weight;
    double value;
  };
  using knapsack_solution = std::vector<bool>;

  class knapsack_problem{
    private:
      std::vector<knapsack_item> items;
      double max_capacity;
    public:
      knapsack_problem(const std::vector<knapsack_item> &item_list, double capacity)
      :
        items(item_list), max_capacity(capacity)
      {}

      double fitness_function(const knapsack_solution& solution) const{
        double total_value = 0.0;
        double total_weight= 0.0;

        for (size_t i = 0; i < solution.size(); ++i)  {
          if (solution[i]){
            total_value += items[i].value;
            total_weight += items[i].weight;
          }
        }
        // penalize solutions that exceed capacity
        if (total_weight > max_capacity){
          return 0.0;
        }
        return total_value;
      }

      knapsack_solution crossover_function(const knapsack_solution& parent1, const knapsack_solution& parent2) const{
        knapsack_solution offspring = parent1;

        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<size_t> dist(0, parent1.size()-1);

        // single-point crossover
        size_t crossover_point = dist(rng);
        for (size_t i = crossover_point; i < offspring.size(); ++i){
          offspring[i] = parent2[i];
        }
        return offspring;
      }

      void mutation_function(knapsack_solution& solution) const{
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<size_t> dist(0, solution.size()-1);

        // flip random bit
        size_t mutation_point = dist(rng);
        solution[mutation_point] = !solution[mutation_point];
      }
      knapsack_solution initialize_function() const{
        knapsack_solution solution(items.size());
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, 1);

        for (size_t i = 0; i < items.size(); ++i){
          solution[i] = (dist(rng) == 1);
        }
        return solution;
      }
  };
 }  // namespace example_problem
 


int main(){
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
      500,
      [&knapsack](const knapsack_solution& sol) { return knapsack.fitness_function(sol);  },
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
}

