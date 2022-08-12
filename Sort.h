//
// Created by andre on 7/30/22.
//

#ifndef NSGAII__SORT_H_
#define NSGAII__SORT_H_

/// Asignamos el crowding distance
/// \param F
/// \param population
void crowding_distance_assignment(std::vector<size_t> &F, std::vector<Individual> &population){
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
    population[F[0]].crowding_distance = population[F[r-1]].crowding_distance = INF;

    // Solo si hay más de dos puntos entonces podemos dar crowding distance
    if(r > 3){
      // del 1,...,r-2
      for(size_t index_of_front = 1; index_of_front < r-1; ++index_of_front){
        // Si f_max y f_min son iguales la distancia se aumenta en nada
        if(std::abs(population[F[r-1]].objectives[index_objective] - population[F[0]].objectives[index_objective]) < EPS){
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
std::vector<std::vector<size_t>> fast_non_dominated_sort(std::vector<Individual> &population){
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
std::vector<Individual> merge(const std::vector<Individual> &parents, const std::vector<Individual> &offspring){

  // R en el artículo
  std::vector<Individual> mixed_population = parents;

  // Mezclar
  mixed_population.insert(mixed_population.end(), offspring.begin(), offspring.end());

  if((mixed_population.size()) != 2*global_problem.population_size)
    throw std::runtime_error("La mezcla de padres e hijos definitivamente no es 2N!!!!");

  return mixed_population;
}

/// Llenamos los nuevos padres como pide NSGA-II aquí está el elitismo
/// \param mixed_population
void fill_new_parents(std::vector<Individual> &mixed_population){
  auto F = fast_non_dominated_sort(mixed_population);

  // Ciclo para llenar P_t+1 en el artículo
  std::vector<Individual> parent_population;
  size_t index_of_F= 0;
  while(parent_population.size() + F[index_of_F].size() <= global_problem.population_size){
    crowding_distance_assignment(F[index_of_F], mixed_population);
    for(auto &p: F[index_of_F]){
      parent_population.push_back(mixed_population[p]);
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
  while(parent_population.size() != global_problem.population_size){
    parent_population.push_back(mixed_population[*itF_l]);
    ++itF_l;
  }

  // Última verificación
  if(parent_population.size() != global_problem.population_size)
    throw std::runtime_error("Algo ocurrió y no se llenó al tamaño deseado");

  mixed_population = parent_population;
}

#endif //NSGAII__SORT_H_
