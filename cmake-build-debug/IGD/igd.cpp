#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "omp.h"

int main(int argc, char *argv[]){
  std::filesystem::path path_pareto, path_algorithm;
  if( argc == 1){
    std::cout << "Usage: pareto_file_data algorithm_file_data";
  } else if(argc == 3){
    path_pareto = argv[1];
    path_algorithm = argv[2];
  } else{
    std::cout << "Usage: pareto_file_data algorithm_file_data";
    throw std::runtime_error("Mal uso del ejecutable");
  }

  if(!is_regular_file(path_pareto) || !is_regular_file(path_algorithm))
    throw std::runtime_error("No son archivos!!");

  // Debug
//  std::cout << path_pareto << " es un " << is_regular_file(path_pareto) << "\n";
//  std::cout << path_algorithm << " es un " << is_regular_file(path_algorithm) << "\n";

  std::ifstream pareto_front_file(path_pareto.c_str());
  std::ifstream algorithm_front_file(path_algorithm.c_str());

  std::vector<std::vector<double>> pareto_front, algorithm_front;

  std::string line;
  while (std::getline(pareto_front_file, line)) {
    std::stringstream ss(line);
    std::vector<double> point;
    while (!ss.eof()){
      double component;
      std::string tmp;

      ss >> tmp;
      if (std::stringstream(tmp) >> component)
        point.push_back(component);
    }
    pareto_front.push_back(point);
  }

  while (std::getline(algorithm_front_file, line)) {
    std::stringstream ss(line);
    std::vector<double> point;
    while (!ss.eof()){
      double component;
      std::string tmp;

      ss >> tmp;
      if (std::stringstream(tmp) >> component)
        point.push_back(component);
    }
    algorithm_front.push_back(point);
  }

  pareto_front_file.close();
  algorithm_front_file.close();

  // Debug
//  for(auto &point: pareto_front){
//    for(auto &component: point){
//      std::cout << component << " ";
//    }
//    std::cout << "\n";
//  }
//
//  for(auto &point: algorithm_front){
//    for(auto &component: point){
//      std::cout << component << " ";
//    }
//    std::cout << "\n";
//  }

  // Calculation

  double igd = 0.0;

  #pragma omp parallel for reduction(+:igd)
  for(size_t index_pareto = 0; index_pareto < pareto_front.size(); ++index_pareto){
    auto point = pareto_front[index_pareto];
    // Obtenemos el punto en el frente del algoritmo con menor distancia al punto del frente
    auto min_distance = *std::min_element(algorithm_front.begin(), algorithm_front.end(), [&point](std::vector<double> point_1, std::vector<double> point_2){
      double distance_1 = 0.0;
      for(size_t index_component = 0; index_component < point.size(); ++index_component){
        distance_1 += (point[index_component] - point_1[index_component]) * (point[index_component] - point_1[index_component]);
      }
      distance_1 = sqrt(distance_1);

      double distance_2 = 0.0;
      for(size_t index_component = 0; index_component < point.size(); ++index_component){
        distance_2 += (point[index_component] - point_2[index_component]) * (point[index_component] - point_2[index_component]);
      }
      distance_2 = sqrt(distance_2);

      return distance_1 < distance_2;
    });

    double distance_min = 0.0;
    for(size_t index_component = 0; index_component < point.size(); ++index_component){
      distance_min += (point[index_component] - min_distance[index_component]) * (point[index_component] - min_distance[index_component]);
    }

    igd += distance_min;
  }

  igd = sqrt(igd);

  igd = (1.0/pareto_front.size())*igd;

  std::cout << igd << "\n";

  return 0;
}