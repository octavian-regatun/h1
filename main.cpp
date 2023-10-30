#include <bitset>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <omp.h>
#include <random>
#include <string>
#include <vector>
#include <chrono>

#define M_PI 3.14159265358979323846
#define DIMENSIONS 30

//variable paramaters

#define ITERATIONS 20000
#define PRECISION 5
#define MAX_IMPROVEMENT_ATTEMPTS 1000
#define MUTATION_RATE 0.001
#define INITIAL_TEMPERATURE 100
#define COOLING_RATE 0.3

struct bounds {
  double min, max;
};

int length;

double rastrigin(const std::vector<double> &x, int dimensions) {
  double sum = 10 * dimensions;
  for (int i = 0; i < dimensions; i++)
    sum += x[i] * x[i] - 10 * cos(2 * M_PI * x[i]);
  return sum;
}

#define M 10

double michalewicz(const std::vector<double> &x, int dimensions) {
  double sum = 0;
  for (int i = 0; i < dimensions; i++)
    sum -= sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / M_PI), 2 * M);
  return sum;
}

double dejong1(const std::vector<double> &x, int dimensions) {
  double sum = 0.0;
  for (int i = 0; i < dimensions; i++)
    sum += x[i] * x[i];
  return sum;
}

double schwefel(const std::vector<double> &x, int dimensions) {
  double sum = 0.0;
  for (int i = 0; i < dimensions; i++) {
    sum += x[i] * sin(sqrt(abs(x[i])));
  }
  return 418.9829 * dimensions - sum;
}

int calculate_length(const bounds &bounds, int dimensions) {
  return static_cast<int>(std::ceil(
      dimensions * std::log2(pow(10, PRECISION) * (bounds.max - bounds.min))));
}

std::random_device random_device;
std::mt19937 random_generator(omp_get_thread_num() + random_device());

std::vector<char> generate_random_bits(int length) {
  std::uniform_int_distribution<int> distribution(0, 1);

  std::vector<char> bits;

  for (int i = 0; i < length; i++) {
    int bit = distribution(random_generator);

    if (bit == 0)
      bits.push_back('0');
    else
      bits.push_back('1');
  }

  return bits;
}

void flipBit(std::vector<char> &bits, int index) {
  auto &bit = bits.at(index);

  if (bit == '0')
    bit = '1';
  else
    bit = '0';
}

std::vector<double> bits_to_number(const std::vector<char> &bits,
                                   const bounds &bounds, int length,
                                   int dimensions) {
  std::vector<double> result;
  int bits_per_dimension = length / dimensions;

  for (int d = 0; d < dimensions; d++) {
    double sum = 0.0;
    for (int i = 0; i < bits_per_dimension; i++) {
      sum = sum * 2 + (bits[d * bits_per_dimension + i] - '0');
    }

    double real_val = bounds.min + sum * (bounds.max - bounds.min) /
                                       (pow(2, bits_per_dimension) - 1);
    result.push_back(real_val);
  }

  return result;
}

double calculate_function(const std::vector<char> &bits,
                          double (*func)(const std::vector<double> &, int),
                          const bounds &bounds, int dimensions) {
  std::vector<double> values = bits_to_number(bits, bounds, length, dimensions);
  return func(values, dimensions);
}

std::vector<char> best_improvement(const std::vector<char> &vc,
                                   double (*func)(const std::vector<double> &,
                                                  int),
                                   const bounds &bounds, int dimensions) {
  std::vector<char> best_neighbour = vc;
  double best_score =
      calculate_function(best_neighbour, func, bounds, dimensions);
  // int bits_per_dimension = length / dimensions;

  /*for (int i = 0; i < bits_per_dimension; i++) {
      std::vector<char> neighbour = vc;
      for (int d = 0; d < dimensions; d++) {
          flipBit(neighbour, d * bits_per_dimension + i);
      }

      double score = calculate_function(neighbour, func, bounds, dimensions);
      if (score < best_score) {
          best_neighbour = neighbour;
          best_score = score;
      }
  }*/

  for (int i = 0; i < length; i++) {
    std::vector<char> neighbour = vc;
    flipBit(neighbour, i);

    double score = calculate_function(neighbour, func, bounds, dimensions);
    if (score < best_score) {
      best_neighbour = neighbour;
      best_score = score;
    }
  }

  return best_neighbour;
}

std::vector<char> first_improvement(const std::vector<char> &vc,
                                    double (*func)(const std::vector<double> &,
                                                   int),
                                    const bounds &bounds, int dimensions) {
  std::vector<char> current = vc;
  double current_score = calculate_function(current, func, bounds, dimensions);
  // int bits_per_dimension = length / dimensions;

  // for (int i = 0; i < bits_per_dimension; i++) {
  //     std::vector<char> neighbour = vc;
  //     for (int d = 0; d < dimensions; d++) {
  //         flipBit(neighbour, d * bits_per_dimension + i);
  //     }

  //    double score = calculate_function(neighbour, func, bounds, dimensions);
  //    if (score < current_score) {
  //        return neighbour; // First improvement found, return it immediately.
  //    }
  //}

  for (int i = 0; i < length; i++) {
    std::vector<char> neighbour = vc;
    flipBit(neighbour, i);

    double score = calculate_function(neighbour, func, bounds, dimensions);
    if (score < current_score) {
      return neighbour;
    }
  }

  return current; // If no improvement is found, return the original solution.
}

std::vector<char> worst_improvement(const std::vector<char> &vc,
                                    double (*func)(const std::vector<double> &,
                                                   int),
                                    const bounds &bounds, int dimensions) {
  std::vector<char> worst_neighbour = vc;
  double worst_score =
      calculate_function(worst_neighbour, func, bounds, dimensions);
  // int bits_per_dimension = length / dimensions;

  /*for (int i = 0; i < bits_per_dimension; i++) {
      std::vector<char> neighbour = vc;
      for (int d = 0; d < dimensions; d++) {
          flipBit(neighbour, d * bits_per_dimension + i);
      }

      double score = calculate_function(neighbour, func, bounds, dimensions);
      if (score > worst_score) {
          worst_neighbour = neighbour;
          worst_score = score;
      }
  }*/
  for (int i = 0; i < length; i++) {
    std::vector<char> neighbour = vc;
    flipBit(neighbour, i);

    double score = calculate_function(neighbour, func, bounds, dimensions);
    if (score > worst_score) {
      worst_neighbour = neighbour;
      worst_score = score;
    }
  }
  return worst_neighbour;
}

enum { BEST = 1, FIRST = 2, WORST = 3 };

std::vector<char> hill_climbing(const std::vector<char> &x,
                                double (*func)(const std::vector<double> &,
                                               int),
                                const bounds &bounds, int dimensions, int type,
                                double mutation_rate) {
  std::vector<char> best = x;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

#pragma omp parallel for
  for (int t = 0; t < ITERATIONS; t++) {
    std::vector<char> current = generate_random_bits(length);

    bool local = false;
    int improvement_attempts = 0;
    std::vector<char> next;

    while (!local && improvement_attempts < MAX_IMPROVEMENT_ATTEMPTS) {
      if (type == BEST)
        next = best_improvement(current, func, bounds, dimensions);
      else if (type == FIRST)
        next = first_improvement(current, func, bounds, dimensions);
      else
        next = worst_improvement(current, func, bounds, dimensions);

      // Apply mutation with probability mutation_rate
      for (int i = 0; i < length; i++) {
        if (distribution(gen) < mutation_rate) {
          flipBit(next, i);
        }
      }

      if (calculate_function(next, func, bounds, dimensions) <=
          calculate_function(current, func, bounds, dimensions)) {
        current = next;
      } else {
        local = true;
      }
      improvement_attempts++;
    }

    if (calculate_function(current, func, bounds, dimensions) <
        calculate_function(best, func, bounds, dimensions)) {
      best = current;
    }
  }
  return best;
}

std::vector<char> generate_random_neighbour(const std::vector<char> &bits,
                                            bounds &bounds, int dimensions) {
  int bits_per_dimension = length / dimensions;
  std::uniform_int_distribution<int> distribution1(0, bits_per_dimension - 1);
  int i = distribution1(random_generator);
  std::vector<char> neighbour(bits.begin(), bits.end());

  for (int d = 0; d < dimensions; d++)
    flipBit(neighbour, d * bits_per_dimension + i);

  return neighbour;
}


std::vector<char>
hybrid_simulated_annealing(const std::vector<char> &x,
                           double (*func)(const std::vector<double> &, int),
                           const bounds &bounds, int dimensions) {
  double initial_temperature = INITIAL_TEMPERATURE;
  double cooling_rate = COOLING_RATE;
  const int MAX_INNER_ITERATIONS = 20;
  double min_temperature = 1e-4;
  double mutation_rate = MUTATION_RATE;
  std::uniform_real_distribution<double> distribution(0.0, 1.0);

  double global_value = calculate_function(x, func, bounds, dimensions);
  std::vector<char> global = x;

  omp_set_num_threads(omp_get_max_threads());

  #pragma omp parallel for shared(global, global_value)
  for (int it = 0; it < 50; it++) {
    int t = 0;
    double T = initial_temperature;
    std::vector<char> vc = x;
    double vc_value = calculate_function(vc, func, bounds, dimensions);

    while (T > min_temperature) {

      for(int inner_iterations = 0;inner_iterations < MAX_INNER_ITERATIONS; inner_iterations++) {
        std::vector<char> vn = best_improvement(vc, func, bounds, dimensions);
        double vn_value = calculate_function(vn, func, bounds, dimensions);

        if (vn_value < vc_value)
         {
          vc = vn;
          vc_value = vn_value;
          break;
        }
        else if(distribution(random_generator) < std::exp(-std::abs(vn_value - vc_value) / T))
        {
            vc = vn;
            vc_value = vn_value;
            break;
        }
      }

      T *= cooling_rate;
      t++;
    }

    std::vector<char> refined_solution =
        hill_climbing(vc, func, bounds, dimensions, BEST, mutation_rate);

    double refined_solution_value =
        calculate_function(refined_solution, func, bounds, dimensions);
    std::vector<char> local =
        (refined_solution_value < vc_value) ? refined_solution : vc;
    double local_value = calculate_function(local, func, bounds, dimensions);

    #pragma omp critical  // This ensures that this section is executed by one thread at a time
    {
      if (local_value < global_value) {
        global_value = local_value;
        global = local;
      }
    }
  }

  return global;
}


int main() {
  bounds dejong1_bounds = {-5.12, 5.12};
  bounds schwefel_bounds = {-500, 500};
  bounds rastrigin_bounds = {-5.12, 5.12};
  bounds michalewicz_bounds = {0, M_PI};

  // Define the dimensions for testing.
  std::vector<int> dimensions = {30};
  std::chrono::time_point<std::chrono::system_clock> start, finish;
  std::chrono::duration<double> elapsed_seconds;

  for (int dim : dimensions) {
    printf("Results for %d-dimensional functions:\n", dim);

    // // Run Hill Climbing with Best Improvement for De Jong 1
    // length = calculate_length(dejong1_bounds, dim);
    // auto dejong1_bits = generate_random_bits(length);
    // auto dejong1_result_bits = hill_climbing(dejong1_bits, dejong1,
    // dejong1_bounds, dim, BEST, MUTATION_RATE); // here double dejong1_result
    // = calculate_function(dejong1_result_bits, dejong1, dejong1_bounds, dim);
    // printf("De Jong 1 (Hill Climbing - Best Improvement): %f\n",
    // dejong1_result);

    // // Run Hill Climbing with Best Improvement for Schwefel's
    // length = calculate_length(schwefel_bounds, dim);
    // auto schwefel_bits = generate_random_bits(length);
    // auto schwefel_result_bits = hill_climbing(schwefel_bits, schwefel,
    // schwefel_bounds, dim, BEST, MUTATION_RATE); //here double schwefel_result
    // = calculate_function(schwefel_result_bits, schwefel, schwefel_bounds,
    // dim); printf("Schwefel's (Hill Climbing - Best Improvement): %f\n",
    // schwefel_result);

    // Run Hill Climbing with Best Improvement for Rastrigin's
    // length = calculate_length(rastrigin_bounds, dim);
    // auto rastrigin_bits = generate_random_bits(length);
    //start = std::chrono::system_clock::now();
    // auto rastrigin_result_bits =
    //     hill_climbing(rastrigin_bits, rastrigin, rastrigin_bounds, dim, BEST,
    //                   MUTATION_RATE); // here
    // double rastrigin_result = calculate_function(
    //     rastrigin_result_bits, rastrigin, rastrigin_bounds, dim);
    // printf("Rastrigin's (Hill Climbing - Best Improvement): %f\n",
    //        rastrigin_result);

    // // Run Hill Climbing with Best Improvement for Michalewicz's
    // length = calculate_length(michalewicz_bounds, dim);
    // auto michalewicz_bits = generate_random_bits(length);
    //start = std::chrono::system_clock::now();
    // auto michalewicz_result_bits = hill_climbing(michalewicz_bits,
    // michalewicz, michalewicz_bounds, dim, BEST, MUTATION_RATE);
    //finish = std::chrono::system_clock::now();
    //elapsed_seconds = finish - start;
    // double michalewicz_result = calculate_function(michalewicz_result_bits,
    // michalewicz, michalewicz_bounds, dim); printf("Michalewicz's (Hill
    // Climbing - Best Improvement): %f\n", michalewicz_result);

    length = calculate_length(rastrigin_bounds, dim);
    auto rastrigin_bits = generate_random_bits(length);
    start = std::chrono::system_clock::now();
    auto rastrigin_result_bits =
        hybrid_simulated_annealing(rastrigin_bits, rastrigin, rastrigin_bounds, dim);
    finish = std::chrono::system_clock::now();
    elapsed_seconds = finish - start;
    double rastrigin_result = calculate_function(
        rastrigin_result_bits, rastrigin, rastrigin_bounds, dim);
    printf("Rastrigin's (Simulated Annealing): %f\nTime: %f",
           rastrigin_result, elapsed_seconds);
    printf("\n");
  }

  return 0;
}
