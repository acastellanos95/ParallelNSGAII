//
// Created by andre on 7/28/22.
//

#ifndef NSGAII__PROBLEM_H_
#define NSGAII__PROBLEM_H_

/* Estructura que mantiene información del problema */
struct Problem {
  double pcross_real;
  double pmut_real;
  double eta_c;
  double eta_m;
  size_t number_of_real_variables;
  size_t number_of_objectives;
  size_t number_of_constraints;
  size_t population_size;
  size_t number_of_generations;
  size_t problem_index;
  size_t migration_interval;
  std::vector<std::pair<double, double>> x_ranges;
};

#endif //NSGAII__PROBLEM_H_
