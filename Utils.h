//
// Created by andre on 7/28/22.
//

#ifndef NSGAII__UTILS_H_
#define NSGAII__UTILS_H_

#include <vector>
#include <random>
#include <algorithm>
#include "Individual.h"
#include "Problem.h"



/// Evaluar problema elegido zdt2
/// \param individuo
void test_problem(Individual &individuo){
  double f1, f2, g, h;
  int i;
  f1 = individuo.x[0];
  g = 0.0;
  for (i=1; i<30; i++)
  {
    g += individuo.x[i];
  }
  g = 9.0*g/29.0;
  g += 1.0;
  h = 1.0 - pow((f1/g),2.0);
  f2 = g*h;
  individuo.objectives[0] = f1;
  individuo.objectives[1] = f2;
}

/// Evaluamos los atributos de cada individuo
/// \param population
void evaluate_population(std::vector<Individual> &population){
  for(auto &individual: population){
    test_problem(individual);

    // Valor de restricciones violadas
    individual.constraint_violation = 0.0;
    for(auto &constraint_value: individual.constraints){
      if(constraint_value < 0.0){
        individual.constraint_violation += constraint_value;
      }
    }
  }
}

/// Inicializamos a la población
/// \param population
void initialize_population(std::vector<Individual> &population){
  for(auto &individual: population){
    individual.x = std::vector<double>(global_problem.number_of_real_variables, 0.0);
    individual.objectives = std::vector<double>(global_problem.number_of_objectives, 0.0);
    individual.constraints = std::vector<double>(global_problem.number_of_constraints, 0.0);
    for(size_t index_of_xi = 0; index_of_xi < individual.x.size(); ++index_of_xi){
      individual.x[index_of_xi] = rndreal(global_problem.x_ranges[index_of_xi].first, global_problem.x_ranges[index_of_xi].second);
    }
  }
}

/// Mutación polinomial
/// \param new_population
void mutation(std::vector<Individual> &new_population) {
  // Iteramos la población
  for (auto &individual: new_population) {
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;

    for (size_t indexOfXi = 0; indexOfXi < individual.x.size(); indexOfXi++) {

      // Si la moneda cargada es verdadera entonces hacemos mutación polinomial
      if (flip(global_problem.pmut_real)) {
        y = individual.x[indexOfXi];
        yl = global_problem.x_ranges[indexOfXi].first;
        yu = global_problem.x_ranges[indexOfXi].second;
        delta1 = (y - yl) / (yu - yl);
        delta2 = (yu - y) / (yu - yl);
        rnd = randomperc();
        mut_pow = 1.0 / (global_problem.eta_m + 1.0);
        if (rnd <= 0.5) {
          xy = 1.0 - delta1;
          val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (global_problem.eta_m + 1.0)));
          deltaq = std::pow(val, mut_pow) - 1.0;
        } else {
          xy = 1.0 - delta2;
          val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow(xy, (global_problem.eta_m + 1.0)));
          deltaq = 1.0 - (std::pow(val, mut_pow));
        }
        y = y + deltaq * (yu - yl);
        if (y < yl)
          y = yl;
        if (y > yu)
          y = yu;
        individual.x[indexOfXi] = y;
      }
    }
  }
}

/// Cruza SBX
/// \param parent_1
/// \param parent_2
/// \param child_1
/// \param child_2
void crossover(const Individual &parent_1, const Individual &parent_2, Individual &child_1, Individual &child_2) {
  double rand;
  double y1, y2, yl, yu;
  double c1, c2;
  double alpha, beta, betaq;

  // Igualamos para deshacernos de condiciones innecesarias en el código de Deb
  child_1.x = parent_1.x;
  child_2.x = parent_2.x;

  // Lanzamos moneda cargada para saber si hacemos cruza o no
  if (flip(global_problem.pcross_real)) {
    for (size_t indexOfXi = 0; indexOfXi < parent_1.x.size(); indexOfXi++) {
      // Igualamos para deshacernos de dos condiciones innecesarias en el código de Deb
      child_1.x[indexOfXi] = parent_1.x[indexOfXi];
      child_2.x[indexOfXi] = parent_2.x[indexOfXi];

      // Lancemos una moneda justa para saber si modificamos la variable x_i
      if (flip(0.5)) {

        // Solo si no son iguales realizamos la cruza
        if (std::abs(parent_1.x[indexOfXi] - parent_2.x[indexOfXi]) > EPS) {
          if (parent_1.x[indexOfXi] < parent_2.x[indexOfXi]) {
            y1 = parent_1.x[indexOfXi];
            y2 = parent_2.x[indexOfXi];
          } else {
            y1 = parent_2.x[indexOfXi];
            y2 = parent_1.x[indexOfXi];
          }
          yl = global_problem.x_ranges[indexOfXi].first;
          yu = global_problem.x_ranges[indexOfXi].second;
          rand = randomperc();
          beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
          alpha = 2.0 - pow(beta, -(global_problem.eta_c + 1.0));
          if (rand <= (1.0 / alpha)) {
            betaq = std::pow((rand * alpha), (1.0 / (global_problem.eta_c + 1.0)));
          } else {
            betaq = std::pow((1.0 / (2.0 - rand * alpha)), (1.0 / (global_problem.eta_c + 1.0)));
          }
          c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
          beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
          alpha = 2.0 - pow(beta, -(global_problem.eta_c + 1.0));
          if (rand <= (1.0 / alpha)) {
            betaq = std::pow((rand * alpha), (1.0 / (global_problem.eta_c + 1.0)));
          } else {
            betaq = std::pow((1.0 / (2.0 - rand * alpha)), (1.0 / (global_problem.eta_c + 1.0)));
          }
          c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
          if (c1 < yl)
            c1 = yl;
          if (c2 < yl)
            c2 = yl;
          if (c1 > yu)
            c1 = yu;
          if (c2 > yu)
            c2 = yu;
          if (flip(0.5)) {
            child_1.x[indexOfXi] = c2;
            child_2.x[indexOfXi] = c1;
          } else {
            child_1.x[indexOfXi] = c1;
            child_2.x[indexOfXi] = c2;
          }
        }
      }
    }
  }
}

/// Rutina que determina dominancia entre dos individuos (ver sección VI del artículo) TODO: evaluar si normalizaremos
/// \param individuo_1
/// \param individuo_2
/// \return
int check_dominance(const Individual &individuo_1, const Individual &individuo_2) {
  /*
   * 1 si individuo_1 domina a individuo_2
   * -1 si individuo_2 domina a individuo_1
   * 0 si ambos individuos son no dominados
   * */

  // Si ambos individuos violaron las restricciones compararemos por quien violó más
  if (individuo_1.constraint_violation < 0.0 && individuo_2.constraint_violation < 0.0) {
    // Individuo 1 violó menos restricciones
    if (individuo_1.constraint_violation > individuo_2.constraint_violation) {
      return (1);
    } else {
      // Individuo 2 violó menos restricciones
      if (individuo_1.constraint_violation < individuo_2.constraint_violation) {
        return (-1);
      }
        // Ambos violaron igualmente las restricciones
      else {
        return (0);
      }
    }
  } else {
    // Si solo el individuo 1 violó restricciones
    if (individuo_1.constraint_violation < 0.0 && individuo_2.constraint_violation == 0.0) {
      return (-1);
    } else {
      // Si solo el individuo 2 violó restricciones
      if (individuo_1.constraint_violation == 0.0 && individuo_2.constraint_violation < 0.0) {
        return (1);
      }
        // Nadie violó restricciones
      else {
        int flag1;
        int flag2;
        flag1 = 0;
        flag2 = 0;

        for (size_t index_objective_function_values = 0;
             index_objective_function_values < individuo_1.objectives.size(); ++index_objective_function_values) {
          // Individuo 2 en alguna función objetivo es mayor que el individuo 1
          if (individuo_1.objectives[index_objective_function_values]
              < individuo_2.objectives[index_objective_function_values]) {
            flag1 = 1;
          } else {
            // Individuo 1 en alguna función objetivo es mayor que el individuo 2
            if (individuo_1.objectives[index_objective_function_values]
                > individuo_2.objectives[index_objective_function_values]) {
              flag2 = 1;
            }
          }
        }

        // Individuo 2 en alguna(s) función(es) objetivo es mayor que el individuo 1 pero nunca el individuo 1 es mayor que el individuo 2 en todas las funciones objetivo
        // individuo_1 domina a individuo_2
        if (flag1 == 1 && flag2 == 0) {
          return (1);
        } else {
          // Individuo 1 en alguna(s) función(es) objetivo es mayor que el individuo 2 pero nunca el individuo 2 es mayor que el individuo 1 en todas las funciones objetivo
          // individuo_2 domina a individuo_1
          if (flag1 == 0 && flag2 == 1) {
            return (-1);
          }
            // Ambos son no dominantes
          else {
            return (0);
          }
        }
      }
    }
  }

}

/// Torneo entre dos individuos
/// \param individuo_1
/// \param individuo_2
/// \return
Individual tournament(const Individual &individuo_1, const Individual &individuo_2) {

  /*
   * 1 si individuo_1 domina a individuo_2
   * -1 si individuo_2 domina a individuo_1
   * 0 si ambos individuos son no dominados
   * */
  int flag = check_dominance(individuo_1, individuo_2);

  if (flag == 1) {
    return individuo_1;
  }

  if (flag == -1) {
    return individuo_2;
  }

  if (individuo_1.crowding_distance > individuo_2.crowding_distance) {
    return individuo_1;
  }

  if (individuo_2.crowding_distance > individuo_1.crowding_distance) {
    return individuo_2;
  }

  if (flip(0.5)) {
    return individuo_1;
  } else {
    return individuo_2;
  }
}

/// Seleccionamos como en NSGAII con torneo binario y hacemos la cruza
/// \param old_population
/// \param new_population
void selection(const std::vector<Individual> &old_population, std::vector<Individual> &new_population) {

  /* Solo por si las dudas si la nueva población no está inicializada hay que inicializar */
  if (new_population.size() != old_population.size())
    new_population = old_population;

  // Arreglos de índices
  std::vector<size_t> populationIndexes_1(old_population.size(), 0);
  std::vector<size_t> populationIndexes_2(old_population.size(), 0);
  for (size_t index_population = 0; index_population < old_population.size(); ++index_population) {
    populationIndexes_1[index_population] = index_population;
    populationIndexes_2[index_population] = index_population;
  }

  // Usando el Rseed determinamos la semilla del generador de aleatorios
  auto seed = (unsigned) round(4294967295.0 * Rseed);
  std::mt19937 g(seed);

  // Barajear los arreglos
  std::shuffle(populationIndexes_1.begin(), populationIndexes_1.end(), g);
  std::shuffle(populationIndexes_2.begin(), populationIndexes_2.end(), g);

  // Selección por torneo binario determinístico
  for (size_t index_population = 0; index_population < old_population.size(); index_population += 4) {
    auto parent_1 = tournament(old_population[populationIndexes_1[index_population]],
                               old_population[populationIndexes_1[index_population + 1]]);
    auto parent_2 = tournament(old_population[populationIndexes_1[index_population + 2]],
                               old_population[populationIndexes_1[index_population + 3]]);
    crossover(parent_1, parent_2, new_population[index_population], new_population[index_population + 1]);
    parent_1 = tournament(old_population[populationIndexes_2[index_population]],
                          old_population[populationIndexes_2[index_population + 1]]);
    parent_2 = tournament(old_population[populationIndexes_2[index_population + 2]],
                          old_population[populationIndexes_2[index_population + 3]]);
    crossover(parent_1, parent_2, new_population[index_population + 2], new_population[index_population + 3]);
  }
}

#endif //NSGAII__UTILS_H_
