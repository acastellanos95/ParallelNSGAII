//
// Created by andre on 7/28/22.
//

#ifndef NSGAII__INDIVIDUAL_H_
#define NSGAII__INDIVIDUAL_H_

#include <vector>

/* Estructura de individuo */
struct Individual {
  std::vector<size_t> solutions_dominated; /* Soluciones dominadas S_p */
  int np; /* Contador de dominación */
  int rank; /* Rango del individuo */
  std::vector<double> x; /* Cromosoma es el vector x */
  std::vector<double> objectives; /* Valor de cada función objetivo */
  std::vector<double> constraints; /* Valor de restricción */
  double constraint_violation; /* Valor de violación de restricciones */
  double crowding_distance; /* Valor de crowding-distance */
};

#endif //NSGAII__INDIVIDUAL_H_
