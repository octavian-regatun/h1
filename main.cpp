#include <bitset>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#define M_PI 3.14159265358979323846
#define DIMENSIONS 30
#define ITERATIONS 10000
#define EPSILON 0.001
#define PRECISION 5

struct bounds {
  double min, max;
};

long long int length;

double rastrigin(const std::vector<double>& x, int dimensions) {
  double sum = 10 * dimensions;
  for (auto i = 0; i < dimensions; i++)
    sum += x[i] * x[i] - 10 * cos(2 * M_PI * x[i]);

  return sum;
}

#define M 10

double michalewicz(const std::vector<double>& x, int dimensions) {
    double sum = 0;
    for (auto i = 0; i < dimensions; i++)
        sum -= sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / M_PI), 2 * M);

    return sum;
}

double dejong1(const std::vector<double>& x, int dimensions)
{
    double sum = 0.0;
    for (auto i = 0; i < dimensions; i++)
        sum += x[i] * x[i];
    return sum;
}

double schwefel(const std::vector<double>& x, int dimensions) {
    double sum = 0.0;
    for (int i = 0; i < dimensions; i++) {
        sum += x[i] * sin(sqrt(abs(x[i])));
    }
    return 418.9829 * dimensions - sum;
}

long long int calculate_length(bounds bounds,int dimensions)
{
  return std::ceil(dimension * std::log2(pow(10, 5) * (rastrigin_bounds.max - rastrigin_bounds.min)));
}

std::random_device random_device;
std::mt19937 random_generator(random_device());

std::vector<char> generate_random_bits() {
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

std::vector<std::vector<char>> generate_neighbours(const std::vector<char>& bits) {
  std::vector<std::vector<char>> neighbours;

  for (int i = 0; i < bits.size(); i++) {
    std::vector<char> neighbour(bits.begin(), bits.end());
    flipBit(neighbour, i);
    neighbours.push_back(neighbour);
  }

  return neighbours;
}

double bits_to_number(const std::vector<char>& bits, bounds bounds) {
  double result = 0.0;

  for (char bit : bits) {
    result = result * 2 + (bit - '0');
  }

  return bounds.min +
         result * (bounds.max - bounds.min) / (pow(2, length) - 1);
}

double calculate_function(const std::vector<char>& bits, double (*func)(const vector<double)>&, int) , bounds bounds) {
  double number = bits_to_number(bits, bounds, length);
  auto function_value = func(number);
  return function_value;
}

std::vector<char> improve(const std::vector<char>& vc, double (*func)(const vector<double)>&, int), bounds bounds)
{
  std::vector<char> best_neighbour = vc;
  double best_score = calculate_function(best_neighbour, func, bounds);
  for(int i = 0; i < length, i++)
  {
    std::vector<char> neighbour = vc;
    neighbour.flip(i);

    double score = func(neighbour); 
    if(score < best_score)
    {
      best_neighbour = neighbour;
      best_score = score;
    }
  }
  return best_neighbour;
}

void hill_climbing(std::vector<double> x, double (*func)(const std::vector<double>&, int), std::vector<double>
    (*grad_func)(const std::vector<double>&, int), bounds bounds, int dimensions) {
    int best = INFINITY;
    bool local;
    for(int t = 0; t < ITERATIONS; t++) {
      local = false;
      std::vector<char> = generate_random_bits(length);
      while(!local)
      {
        std::vector<char> next = improve(current);
        if(calculate_function(next, func, bounds) >= calculate_function(current,func, bounds))
          next = current;
        else 
          local = true;
      }

      if(func(current) < func(best)) best = current;
    }

}

int main() {
  bounds rastrigin_bounds = {-5.12, 5.12};
  bounds michalewicz_bounds = {0, M_PI};
  bounds dejong1_bounds = {-5.12, 5.12};
  bounds schwefel_bounds = {-500, 500};

  auto bits = generate_random_bits(number_of_bits);

  auto neighbours = generate_neighbours(bits);

  auto initial_function_value =
      calculate_function(bits, rastrigin_bounds, number_of_bits);
  printf("Initial Function Value: %f\n\n", initial_function_value);

  for (auto neighbour : neighbours) {
    auto gradient =
        calculate_gradient(neighbour, rastrigin_bounds, number_of_bits);
    auto function_value =
        calculate_function(neighbour, rastrigin_bounds, number_of_bits);
    printf("Gradient: %f\n", gradient.at(0));
    printf("Function Value: %f\n\n", function_value);
  }

  return 0;
}