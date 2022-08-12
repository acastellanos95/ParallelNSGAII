#include <iostream>
#include <iomanip>
#include <queue>
#include <fstream>
#include <random>
#include "Individual.h"
#include "Problem.h"
#include "NSGA.h"
#include "omp.h"


int main(int argc, char *argv[]) {
//  Pedimos los valores
  size_t problemIndex, number_demes, size_demes, number_of_generations, migration_interval; /* Usaremos esto y no una probabilidad porque así corre más tiempo sin bloquear hilos además de que son recíprocas: http://staffwww.dcs.shef.ac.uk/people/d.sudholt/parallel-eas.pdf */
  std::string report_filename;
  std::vector<NSGA> demes;
  if (argc == 1) {
    std::cout << "Introduzca índice de problema (1) DTLZ1, (2) DTLZ2, (3) DTLZ4, (4) WFG1, (5) WFG2, (6) Un problema real: ";
    std::cin >> problemIndex;
    std::cout << "Introduzca número de islas (máximo " << omp_get_max_threads() << "): ";
    std::cin >> number_demes;
    omp_set_num_threads((int) number_demes);
    std::cout << "Tamaño de población de las islas: ";
    std::cin >> size_demes;
    std::cout << "Introduzca el número máximo de generaciones: ";
    std::cin >> number_of_generations;
    std::cout << "Intervalo de migración (debe ser positivo): ";
    std::cin >> migration_interval;

    // Inicializar islas con número de islas, tamaño de población e intervalo de migración
    demes = std::vector<NSGA>(number_demes);
    for(auto &deme: demes){
      deme.problem.population_size = size_demes;
      deme.problem.problem_index = problemIndex;
      deme.problem.number_of_generations = number_of_generations;
      deme.problem.migration_interval = migration_interval;
    }

    // Múltiples probabilidades
    u_int8_t is_multiple_probabilities;
    std::cout << "¿Quiere diferentes probabilidades y densidades por isla? (1) Si, (0) No: ";
    std::cin >> is_multiple_probabilities;

    if(is_multiple_probabilities == 1){
      uint8_t deme_number = 1;
      for(auto &deme: demes){
        std::cout << "Probabilidad de cruza para la isla " << deme_number << " : ";
        std::cin >> deme.problem.pcross_real;
        std::cout << "Probabilidad de mutación para la isla " << deme_number << " : ";
        std::cin >> deme.problem.pmut_real;
        std::cout << "Densidad de cruza para la isla " << deme_number << " : ";
        std::cin >> deme.problem.eta_c;
        std::cout << "Densidad de mutación para la isla " << deme_number << " : ";
        std::cin >> deme.problem.eta_m;
        deme_number++;
      }
    } else {
      double pcross_real;
      double pmut_real;
      double eta_c;
      double eta_m;
      std::cout << "Probabilidad de cruza: ";
      std::cin >> pcross_real;
      std::cout << "Probabilidad de mutación: ";
      std::cin >> pmut_real;
      std::cout << "Densidad de cruza: ";
      std::cin >> eta_c;
      std::cout << "Densidad de mutación: ";
      std::cin >> eta_m;

      for(auto &deme: demes){
        deme.problem.pcross_real = pcross_real;
        deme.problem.pmut_real = pmut_real;
        deme.problem.eta_c = eta_c;
        deme.problem.eta_m = eta_m;
      }
    }

    // Introducir semillas para cada generador
    u_int8_t is_seed_user_generated;
    std::cout << "Quiere asignar usted las semillas de cada isla (Si selecciona no se usa mersenne twister)? (1) Si, (0) No: ";
    std::cin >> is_seed_user_generated;
    if(is_seed_user_generated == 1) {
      std::cout << "Introduzca las semillas de cada isla: \n";
      uint8_t deme_number = 1;
      for(auto &deme: demes){
        std::cout << "Semilla para la isla " << (int) deme_number << " : ";
        std::cin >> deme.seed;
        deme_number++;
      }
    } else if(is_seed_user_generated == 0){
      double fracc;
      std::cout << "Ingresa un número entre 0.0 y 1.0: ";
      std::cin >> fracc;
      auto Gseed = (unsigned) round(4294967295.0 * fracc);
      std::mt19937 g(Gseed);
      std::uniform_real_distribution<float> distribution(0.0,1.0);
      uint8_t deme_number = 1;
      for(auto &deme: demes){
        deme.seed = distribution(g);
        std::cout << "Semilla para la isla " << (int) deme_number << " : " << deme.seed << "\n";
        deme_number++;
      }
    } else
      throw std::runtime_error("Esa opción no existe");

    std::cout << "Nombrar archivo con reporte: ";
    std::cin >> report_filename;
  } else if(argc == 12) {

    double pcross_real;
    double pmut_real;
    double eta_c;
    double eta_m;
    double fracc;
    problemIndex = std::strtoul(argv[1], nullptr, 0);
    number_demes = std::strtoul(argv[2], nullptr, 0);
    size_demes = std::strtoul(argv[3], nullptr, 0);
    number_of_generations = std::strtoul(argv[4], nullptr, 0);
    migration_interval = std::strtoul(argv[5], nullptr, 0);
    pcross_real = std::strtod(argv[6], nullptr);
    pmut_real = std::strtod(argv[7], nullptr);
    eta_c = std::strtod(argv[8], nullptr);
    eta_m = std::strtod(argv[9], nullptr);
    fracc = std::strtod(argv[10], nullptr);
    report_filename = std::string(argv[11]);

    // Inicializar islas con número de islas, tamaño de población e intervalo de migración
    demes = std::vector<NSGA>(number_demes);
    for(auto &deme: demes){
      deme.problem.population_size = size_demes;
      deme.problem.problem_index = problemIndex;
      deme.problem.number_of_generations = number_of_generations;
      deme.problem.migration_interval = migration_interval;
    }

    for(auto &deme: demes){
      deme.problem.pcross_real = pcross_real;
      deme.problem.pmut_real = pmut_real;
      deme.problem.eta_c = eta_c;
      deme.problem.eta_m = eta_m;
    }

    auto Gseed = (unsigned) round(4294967295.0 * fracc);
    std::mt19937 g(Gseed);
    std::uniform_real_distribution<float> distribution(0.0,1.0);
    uint8_t deme_number = 1;
    for(auto &deme: demes){
      deme.seed = distribution(g);
      std::cout << "Semilla para la isla " << (int) deme_number << " : " << deme.seed << "\n";
      deme_number++;
    }
  } else {
    throw std::runtime_error("programa no admite entradas en línea de comandos");
  }

  // Población final de todas las islas
  std::vector<Individual> all_population;
  // Fila para migrar
  std::vector<std::queue<Individual>> migration_queue(number_demes);

  // Ejecución muli-hilo
  #pragma omp parallel
  {
    // Inicializar
    size_t deme_index = omp_get_thread_num();
    demes[deme_index].random_struct.randomize();
    demes[deme_index].parent_population = std::vector<Individual>(demes[deme_index].problem.population_size);
    demes[deme_index].child_population = std::vector<Individual>(demes[deme_index].problem.population_size);

    if(problemIndex == 1 || problemIndex == 2 || problemIndex == 3){
      demes[deme_index].problem.x_ranges = std::vector<std::pair<double, double>>(12, {0.0,1.0});
      demes[deme_index].problem.number_of_objectives = 3;
      demes[deme_index].problem.number_of_real_variables = 12;
      demes[deme_index].problem.number_of_constraints = 0;
    } else if(problemIndex == 4 || problemIndex == 5){
      demes[deme_index].problem.x_ranges = std::vector<std::pair<double, double>>();
      double i_d = 1.0;
      for (size_t index_ranges = 0; index_ranges < 24; ++index_ranges) {
        demes[deme_index].problem.x_ranges.emplace_back(0.0, 2.0 * i_d);
        i_d += 1.0;
      }
      demes[deme_index].problem.number_of_objectives = 3;
      demes[deme_index].problem.number_of_real_variables = 24;
      demes[deme_index].problem.number_of_constraints = 0;
      // Solo por diversión tomé un problema del mundo real de: https://www.sciencedirect.com/science/article/abs/pii/S1568494620300181
    } else if(problemIndex == 6){
      demes[deme_index].problem.x_ranges = std::vector<std::pair<double, double>>(4, {0.0, 1.0});
      demes[deme_index].problem.number_of_objectives = 3;
      demes[deme_index].problem.number_of_real_variables = 4;
      demes[deme_index].problem.number_of_constraints = 0;
    } else {
      throw std::runtime_error("No elegiste un problema!!");
    }

    demes[deme_index].initialize_population(demes[deme_index].parent_population);
    demes[deme_index].initialize_population(demes[deme_index].child_population);
    demes[deme_index].evaluate_population(demes[deme_index].parent_population);
    auto F = demes[deme_index].fast_non_dominated_sort(demes[deme_index].parent_population);

    std::vector<Individual> tmp;
    for(auto &F_i: F){
      for(auto &p: F_i){
        tmp.push_back(demes[deme_index].parent_population[p]);
      }
    }

    demes[deme_index].parent_population = tmp;

    // NSGA-II main loop
    for(size_t index_generation = 2; index_generation <= demes[deme_index].problem.number_of_generations; ++index_generation){
      demes[deme_index].selection(demes[deme_index].parent_population, demes[deme_index].child_population);
      demes[deme_index].mutation(demes[deme_index].child_population);
      demes[deme_index].evaluate_population(demes[deme_index].child_population);
      auto mixed_population = demes[deme_index].merge(demes[deme_index].parent_population, demes[deme_index].child_population);
      if(index_generation % demes[deme_index].problem.migration_interval == 0){
        if(omp_get_thread_num() == 0)
          std::cout << index_generation << "\n";
        // Esto no es paralelo pero es necesario para tener consistencia de datos
        #pragma omp critical
        {
          // Mandamos el mejor de la iteración anterior
          migration_queue[(deme_index + 1) % number_demes].push(demes[deme_index].parent_population[0]);
        }
        #pragma omp barrier
        #pragma omp critical
        {
          // Recibimos al mejor de otra isla y lo ponemos en la mixta antes de ser rankeada
          mixed_population.push_back(migration_queue[deme_index].front());
          migration_queue[deme_index].pop();
        }
        #pragma omp barrier
      }
      demes[deme_index].fill_new_parents(mixed_population);
      demes[deme_index].parent_population = mixed_population;
    }

    // Juntamos las poblaciones
    #pragma omp critical
    {
      // Insertamos
      all_population.insert(all_population.end(), demes[deme_index].parent_population.begin(), demes[deme_index].parent_population.end());
    }
  }

  demes[0].fast_non_dominated_sort(all_population);

  std::vector<Individual> final_dominated_individuals;

  // Solo nos quedamos con los no dominados
  for(auto &individuo: all_population){
    if(individuo.rank == 1){
      final_dominated_individuals.push_back(individuo);
    }
  }

  std::ofstream report(report_filename);

  for(auto &individuo: final_dominated_individuals){
    report << individuo.objectives[0] << "\t" << individuo.objectives[1] << "\t" << individuo.objectives[2] << "\n";
  }

  report.close();

  return 0;
}
