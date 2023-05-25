//
// Created by andre on 8/6/22.
//

#include <random>
#include <algorithm>
#include "NSGA.h"
#include "Problem.h"

/// Evaluar individuo según el problema elegido
/// \param individuo
void NSGA::test_problem(Individual &individuo) {
  // Seleccionamos según el problema seleccionado
  if(problem.problem_index == 1){
    double gx = 10.0;

    for(size_t i = 2; i < individuo.x.size(); ++i){
      gx += ((individuo.x[i] - 0.5)*(individuo.x[i] - 0.5) - std::cos(20*M_PI*(individuo.x[i] - 0.5)));
    }

    gx *= 100.0;

    individuo.objectives[0] = 0.5 * individuo.x[0] * individuo.x[1] * (1 + gx);
    individuo.objectives[1] = 0.5 * individuo.x[0] * (1.0 - individuo.x[1]) * (1 + gx);
    individuo.objectives[2] = 0.5 * (1.0 - individuo.x[0]) * (1 + gx);
  } else if(problem.problem_index == 2){
    double gx = 0.0;

    for(size_t i = 2; i < individuo.x.size(); ++i){
      gx += (individuo.x[i] - 0.5)*(individuo.x[i] - 0.5);
    }

    individuo.objectives[0] = std::cos(0.5 * M_PI * individuo.x[0]) * std::cos(0.5 * M_PI * individuo.x[1]) * (1 + gx);
    individuo.objectives[1] = std::cos(0.5 * M_PI * individuo.x[0]) * std::sin(0.5 * M_PI * individuo.x[1]) * (1 + gx);
    individuo.objectives[2] = std::sin(0.5 * M_PI * individuo.x[0]) * (1 + gx);
  } else if(problem.problem_index == 3){
    double gx = 0.0;

    for(size_t i = 2; i < individuo.x.size(); ++i){
      gx += (individuo.x[i] - 0.5)*(individuo.x[i] - 0.5);
    }

    individuo.objectives[0] = std::cos(0.5 * M_PI * std::pow(individuo.x[0], 100.0)) * std::cos(0.5 * M_PI * std::pow(individuo.x[1], 100.0)) * (1 + gx);
    individuo.objectives[1] = std::cos(0.5 * M_PI * std::pow(individuo.x[0], 100.0)) * std::sin(0.5 * M_PI * std::pow(individuo.x[1], 100.0)) * (1 + gx);
    individuo.objectives[2] = std::sin(0.5 * M_PI * std::pow(individuo.x[0], 100.0)) * (1 + gx);
  } else if(problem.problem_index == 4){

    wfg_evaluator.evaluate_wfg1(individuo.x, individuo.objectives);

  } else if(problem.problem_index == 5){

    wfg_evaluator.evaluate_wfg2(individuo.x, individuo.objectives);

  } else if(problem.problem_index == 6) {
    double xAlpha = individuo.x[0];
    double xHA = individuo.x[1];
    double xOA = individuo.x[2];
    double xOPTT = individuo.x[3];

    // f1 (TF_max)
    individuo.objectives[0] = 0.692 + (0.477 * xAlpha) - (0.687 * xHA) - (0.080 * xOA) - (0.0650 * xOPTT) - (0.167 * xAlpha * xAlpha) - (0.0129 * xHA * xAlpha) + (0.0796 * xHA * xHA) - (0.0634 * xOA * xAlpha) - (0.0257 * xOA * xHA) + (0.0877 * xOA * xOA) - (0.0521 * xOPTT * xAlpha) + (0.00156 * xOPTT * xHA) + (0.00198 * xOPTT * xOA) + (0.0184 * xOPTT * xOPTT);

    // f2 (X_cc)
    individuo.objectives[1] = 0.153 - (0.322 * xAlpha) + (0.396 * xHA) + (0.424 * xOA) + (0.0226 * xOPTT) + (0.175 * xAlpha * xAlpha) + (0.0185 * xHA * xAlpha) - (0.0701 * xHA * xHA) - (0.251 * xOA * xAlpha) + (0.179 * xOA * xHA) + (0.0150 * xOA * xOA) + (0.0134 * xOPTT * xAlpha) + (0.0296 * xOPTT * xHA) + (0.0752 * xOPTT * xOA) + (0.0192 * xOPTT * xOPTT);

    // f3 (TT_max)
    individuo.objectives[2] = 0.370 - (0.205 * xAlpha) + (0.0307 * xHA) + (0.108 * xOA) + (1.019 * xOPTT) - (0.135 * xAlpha * xAlpha) + (0.0141 * xHA * xAlpha) + (0.0998 * xHA * xHA) + (0.208 * xOA * xAlpha) - (0.0301 * xOA * xHA) - (0.226 * xOA * xOA) + (0.353 * xOPTT * xAlpha) - (0.0497 * xOPTT * xOA) - (0.423 * xOPTT * xOPTT) + (0.202 * xHA * xAlpha * xAlpha) - (0.281 * xOA * xAlpha * xAlpha) - (0.342 * xHA * xHA * xAlpha) - (0.245 * xHA * xHA * xOA) + (0.281 * xOA * xOA * xHA) - (0.184 * xOPTT * xOPTT * xAlpha) - (0.281 * xHA * xAlpha * xOA);
  }
}

/// Evaluamos los atributos de cada individuo
/// \param population
void NSGA::evaluate_population(std::vector<Individual> &population) {
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

/// Inicializamos a la población TODO: inicializar objetivos, restricciones
/// \param population
void NSGA::initialize_population(std::vector<Individual> &population) {
  for(auto &individual: population){
    individual.x = std::vector<double>(problem.number_of_real_variables, 0.0);
    individual.objectives = std::vector<double>(problem.number_of_objectives, 0.0);
    individual.constraints = std::vector<double>(problem.number_of_constraints, 0.0);
    for(size_t index_of_xi = 0; index_of_xi < individual.x.size(); ++index_of_xi){
      individual.x[index_of_xi] = random_struct.rndreal(problem.x_ranges[index_of_xi].first, problem.x_ranges[index_of_xi].second);
    }
  }
}

/// Mutación polinomial
/// \param new_population
void NSGA::mutation(std::vector<Individual> &new_population) {
  // Iteramos la población
  for (auto &individual: new_population) {
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;

    for (size_t indexOfXi = 0; indexOfXi < individual.x.size(); indexOfXi++) {

      // Si la moneda cargada es verdadera entonces hacemos mutación polinomial
      if (random_struct.flip(problem.pmut_real)) {
        y = individual.x[indexOfXi];
        yl = problem.x_ranges[indexOfXi].first;
        yu = problem.x_ranges[indexOfXi].second;
        delta1 = (y - yl) / (yu - yl);
        delta2 = (yu - y) / (yu - yl);
        rnd = random_struct.randomperc();
        mut_pow = 1.0 / (problem.eta_m + 1.0);
        if (rnd <= 0.5) {
          xy = 1.0 - delta1;
          val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (problem.eta_m + 1.0)));
          deltaq = std::pow(val, mut_pow) - 1.0;
        } else {
          xy = 1.0 - delta2;
          val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow(xy, (problem.eta_m + 1.0)));
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
void NSGA::crossover(const Individual &parent_1, const Individual &parent_2, Individual &child_1, Individual &child_2) {
  double rand;
  double y1, y2, yl, yu;
  double c1, c2;
  double alpha, beta, betaq;

  // Igualamos para deshacernos de condiciones innecesarias en el código de Deb
  child_1.x = parent_1.x;
  child_2.x = parent_2.x;

  // Lanzamos moneda cargada para saber si hacemos cruza o no
  if (random_struct.flip(problem.pcross_real)) {
    for (size_t indexOfXi = 0; indexOfXi < parent_1.x.size(); indexOfXi++) {
      // Igualamos para deshacernos de dos condiciones innecesarias en el código de Deb
      child_1.x[indexOfXi] = parent_1.x[indexOfXi];
      child_2.x[indexOfXi] = parent_2.x[indexOfXi];

      // Lancemos una moneda justa para saber si modificamos la variable x_i
      if (random_struct.flip(0.5)) {

        // Solo si no son iguales realizamos la cruza
        if (std::abs(parent_1.x[indexOfXi] - parent_2.x[indexOfXi]) > 1.0e-14) {
          if (parent_1.x[indexOfXi] < parent_2.x[indexOfXi]) {
            y1 = parent_1.x[indexOfXi];
            y2 = parent_2.x[indexOfXi];
          } else {
            y1 = parent_2.x[indexOfXi];
            y2 = parent_1.x[indexOfXi];
          }
          yl = problem.x_ranges[indexOfXi].first;
          yu = problem.x_ranges[indexOfXi].second;
          rand = random_struct.randomperc();
          beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
          alpha = 2.0 - pow(beta, -(problem.eta_c + 1.0));
          if (rand <= (1.0 / alpha)) {
            betaq = std::pow((rand * alpha), (1.0 / (problem.eta_c + 1.0)));
          } else {
            betaq = std::pow((1.0 / (2.0 - rand * alpha)), (1.0 / (problem.eta_c + 1.0)));
          }
          c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
          beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
          alpha = 2.0 - pow(beta, -(problem.eta_c + 1.0));
          if (rand <= (1.0 / alpha)) {
            betaq = std::pow((rand * alpha), (1.0 / (problem.eta_c + 1.0)));
          } else {
            betaq = std::pow((1.0 / (2.0 - rand * alpha)), (1.0 / (problem.eta_c + 1.0)));
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
          if (random_struct.flip(0.5)) {
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
int NSGA::check_dominance(const Individual &individuo_1, const Individual &individuo_2) {
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
Individual NSGA::tournament(const Individual &individuo_1, const Individual &individuo_2) {
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

  if (random_struct.flip(0.5)) {
    return individuo_1;
  } else {
    return individuo_2;
  }
}

/// Seleccionamos como en NSGAII con torneo binario y hacemos la cruza
/// \param old_population
/// \param new_population
void NSGA::selection(const std::vector<Individual> &old_population, std::vector<Individual> &new_population) {
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
  auto Gseed = (unsigned) round(4294967295.0 * this->seed);
  std::mt19937 g(Gseed);

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

/// Asignamos el crowding distance
/// \param F
/// \param population
void NSGA::crowding_distance_assignment(std::vector<size_t> &F, std::vector<Individual> &population) {
  size_t r = F.size();

  if(r == 0)
    throw std::runtime_error("Un frente de pareto local que no tiene puntos?!");

  // Inicializando a 0 los que están en esta frontera
  for(auto &i: F){
    population[i].crowding_distance = 0;
  }

  for(size_t index_objective = 0; index_objective < population[0].objectives.size(); ++index_objective){
    // Ordenar por objetivo este frente local
    std::sort(F.begin(), F.end(), [&population, &index_objective](size_t i, size_t j){ return population[i].objectives[index_objective] < population[j].objectives[index_objective];});
    population[F[0]].crowding_distance = population[F[r-1]].crowding_distance = 1.0e14;

    // Solo si hay más de dos puntos entonces podemos dar crowding distance
    // Solo si hay más de dos puntos entonces podemos dar crowding distance
    if(r > 3){
      // del 1,...,r-2
      for(size_t index_of_front = 1; index_of_front < r-1; ++index_of_front){
        // Si f_max y f_min son iguales la distancia se aumenta en nada
        if(std::abs(population[F[r-1]].objectives[index_objective] - population[F[0]].objectives[index_objective]) < 1.0e-14){
          population[F[index_of_front]].crowding_distance += 0.0;
          // En caso contrario podemos usar la fórmula en el artículo
        } else {
          population[F[index_of_front]].crowding_distance += (population[F[index_of_front + 1]].objectives[index_objective] - population[F[index_of_front - 1]].objectives[index_objective]) / (population[F[r-1]].objectives[index_objective] - population[F[0]].objectives[index_objective]);
        }
      }
    }
  }
}

/// Ordenamos los índices por frentes no dominados y regresamos estos grupos indexados
/// \param population
/// \return
std::vector<std::vector<size_t>> NSGA::fast_non_dominated_sort(std::vector<Individual> &population) {
  // Arreglo de conjuntos no dominados
  std::vector<std::vector<size_t>> F;
  // Conjunto i-ésimo no dominado
  std::vector<size_t> F_i;
  // Población final
  std::vector<Individual> final_population;

  // para p en P
  for(size_t index_of_p_in_P = 0; index_of_p_in_P < population.size(); ++index_of_p_in_P){
    // Conjunto de soluciones (índices) que domina p
    population[index_of_p_in_P].solutions_dominated.clear();
    // Inicializar contador de dominación sobre p
    population[index_of_p_in_P].np = 0;

    //para q en P
    for(size_t index_of_q_in_P = 0; index_of_q_in_P < population.size(); ++index_of_q_in_P){
      // Dominancia entre ambos
      int dominance_flag = check_dominance(population[index_of_p_in_P], population[index_of_q_in_P]);
      if(dominance_flag == 1)
        population[index_of_p_in_P].solutions_dominated.push_back(index_of_q_in_P);
      else if(dominance_flag == -1)
        ++population[index_of_p_in_P].np;
    }

    // p pertenece a la primera frente de pareto
    if(population[index_of_p_in_P].np == 0){
      population[index_of_p_in_P].rank = 1;
      F_i.push_back(index_of_p_in_P);
    }
  }

  // Asignamos rangos
  int rank = 1;
  while(!F_i.empty()){
    F.push_back(F_i);
    crowding_distance_assignment(F_i, population);
    std::vector<size_t> F_iplusone;
    // Para cada p en F_i
    for(auto p: F_i){
      // Para cada q en S_p
      for(auto q: population[p].solutions_dominated){
        --population[q].np;

        // Ninguno domina a q
        if(population[q].np == 0){
          population[q].rank = rank;
          F_iplusone.push_back(q);
        }
      }
    }

    // Aumentamos rango y nueva frontera de pareto
    ++rank;
    F_i = F_iplusone;
  }

  return F;
}

/// Mezcla padres e hijos
/// \param parents
/// \param offspring
/// \return
std::vector<Individual> NSGA::merge(const std::vector<Individual> &parents, const std::vector<Individual> &offspring) {
  // R en el artículo
  std::vector<Individual> mixed_population = parents;

  // Mezclar
  mixed_population.insert(mixed_population.end(), offspring.begin(), offspring.end());

  return mixed_population;
}

/// Llenamos los nuevos padres como pide NSGA-II aquí está el elitismo
/// \param mixed_population
void NSGA::fill_new_parents(std::vector<Individual> &mixed_population) {
  auto F = fast_non_dominated_sort(mixed_population);

  // Ciclo para llenar P_t+1 en el artículo
  std::vector<Individual> final_parent_population;
  size_t index_of_F= 0;
  while(final_parent_population.size() + F[index_of_F].size() <= problem.population_size){
    crowding_distance_assignment(F[index_of_F], mixed_population);
    for(auto &p: F[index_of_F]){
      final_parent_population.push_back(mixed_population[p]);
    }
    ++index_of_F;
  }

  // A este último frente le falta su crowding_distance (error en el artículo)
  crowding_distance_assignment(F[index_of_F], mixed_population);
  // Ordenamos por crowding distance en orden ascendente
  std::sort(F[index_of_F].begin(), F[index_of_F].end(), [&mixed_population](size_t i, size_t j){
    return mixed_population[i].crowding_distance < mixed_population[j].crowding_distance;
  });

  // elegimos los últimos hasta llenar al tamaño deseado
  auto itF_l = F[index_of_F].rbegin();
  while(final_parent_population.size() != problem.population_size){
    final_parent_population.push_back(mixed_population[*itF_l]);
    ++itF_l;
  }

  // Última verificación
  if(final_parent_population.size() != problem.population_size)
    throw std::runtime_error("Algo ocurrió y no se llenó al tamaño deseado");

  mixed_population = final_parent_population;
}
