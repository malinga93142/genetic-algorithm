#include "genetic_algorithm.h"
#include "problems.h"
#include "genetic_algorithm.tpp"  // template implementations
#include <iostream>
#include <iomanip>


int main(){
  using namespace genetic_algorithms;
  using namespace example_problems;

  std::vector<city> cities = {
    {0.0, 0.0}, // city 0
    {1.0, 2.0}, // city 1
    {3.0, 1.0}, // city 2
    {2.0, 3.0}, // city 3
    {4.0, 4.0}, // city 4
    {5.0, 2.0}, // city 5
    {6.0, 5.0}, // city 6
    {3.0, 6.0}, // city 7
    {1.0, 5.0}, // city 8
    {2.0, 7.0}
  };

  tsp_problem tsp(cities);

  std::cout << "=== Traveling Salesman Problem Genetic Algorithm ==="<<std::endl;
  std::cout << "Number of cities: " << cities.size() << std::endl;
  std::cout << "===================================================" <<std::endl;

  genetic_algorithm<tsp_solution> ga(
        200,
        0.9,
        0.05,
        5,
        1000,
        [&tsp](const tsp_solution& sol){
          return tsp.fitness_function(sol);
        },
        [&tsp](const tsp_solution& p1, const tsp_solution& p2){
          return tsp.crossover_function(p1, p2);
        },
        [&tsp](tsp_solution& sol){
          return tsp.mutation_function(sol);
        },
        [&tsp](){
          return tsp.initialization_function();
        },
        genetic_algorithm<tsp_solution>::selection_method::TOURNAMENT
      );

  std::cout << "Running genetic algorithm..." << std::endl;
  auto best_solution = ga.run(true);

  // display results
  std::cout << "\n=== FINAL RESULTS ===" << std::endl;
  std::cout << "Best fitness: " << best_solution.fitness << std::endl;
  std::cout << "Best distance:" << std::fixed << std::setprecision(2) << (1.0/best_solution.fitness) << std::endl;


  std::cout << "Best route: ";
  for (size_t i = 0; i < best_solution.chromosome.size(); ++i){
    std::cout << best_solution.chromosome[i];
    if (i < best_solution.chromosome.size()-1){
      std::cout << " -> ";
    }
  }
  std::cout << " -> " << best_solution.chromosome[0] << std::endl;

  double total_distance = 0.0;
  std::cout << "\nRoute details:" << std::endl;
  for (size_t i = 0; i < best_solution.chromosome.size(); ++i){
    size_t from_city  = best_solution.chromosome[i];
    size_t to_city    = best_solution.chromosome[(i+1) % best_solution.chromosome.size()];

    const city &from  = cities[from_city];
    const city &to    = cities[to_city];

    double dx = from.x - to.x;
    double dy = from.y - to.y;
    double segment_distance = std::sqrt(dx*dx + dy*dy);
    total_distance += segment_distance;
    std::cout << "  city " << from_city << " (" << from.x << ", " << from.y
      << ") to city " << to_city << " (" << to.x << ", " << to.y
      << ") distance: " << std::fixed << std::setprecision(2) << segment_distance << std::endl;
  }
  std::cout << "total distance: " << std::fixed << std::setprecision(2) << total_distance << std::endl;
}
