//
// Created by andre on 8/6/22.
//

#ifndef NSGAII__NSGA_H_
#define NSGAII__NSGA_H_

#include "Rand.h"
#include "Individual.h"
#include "Problem.h"
#include "lib/WFG/WFG.h"

class NSGA {
 public:
  float seed;
  Rand random_struct;
  Problem problem;
  WFG wfg_evaluator;
  std::vector<Individual> parent_population, child_population;
  NSGA() = default;
  void test_problem(Individual &individuo);
  void evaluate_population(std::vector<Individual> &population);
  void initialize_population(std::vector<Individual> &population);
  void mutation(std::vector<Individual> &new_population);
  void crossover(const Individual &parent_1, const Individual &parent_2, Individual &child_1, Individual &child_2);
  int check_dominance(const Individual &individuo_1, const Individual &individuo_2);
  Individual tournament(const Individual &individuo_1, const Individual &individuo_2);
  void selection(const std::vector<Individual> &old_population, std::vector<Individual> &new_population);
  void crowding_distance_assignment(std::vector<size_t> &F, std::vector<Individual> &population);
  std::vector<std::vector<size_t>> fast_non_dominated_sort(std::vector<Individual> &population);
  std::vector<Individual> merge(const std::vector<Individual> &parents, const std::vector<Individual> &offspring);
  void fill_new_parents(std::vector<Individual> &mixed_population);
};

#endif //NSGAII__NSGA_H_
