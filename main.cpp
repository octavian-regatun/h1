#include <bitset>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <math.h>
#include <random>
#include <string>
#include <vector>

struct bounds {
  double min, max;
};

double rastrigin(std::vector<double> x) {
  int dimensions = x.size();

  double sum = 10 * dimensions;
  for (auto i = 0; i < dimensions; i++)
    sum += x[i] * x[i] - 10 * cos(2 * M_PI * x[i]);

  return sum;
}

std::vector<double> gradient_rastrigin(std::vector<double> x) {
  std::vector<double> gradients;

  int dimensions = x.size();
  for (auto i = 0; i < dimensions; i++)
    gradients.push_back(2 * x.at(i) + 20 * M_PI * sin(2 * M_PI * x.at(i)));

  return gradients;
}

std::random_device random_device;
std::mt19937 random_generator(random_device());

std::vector<char> generate_random_bits(int number_of_bytes) {
  std::uniform_int_distribution<int> distribution(0, 1);

  std::vector<char> bits;

  for (int i = 0; i < number_of_bytes; i++) {
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

std::vector<std::vector<char>> generate_neighbours(std::vector<char> bits) {
  std::vector<std::vector<char>> neighbours;

  for (int i = 0; i < bits.size(); i++) {
    std::vector<char> neighbour(bits.begin(), bits.end());
    flipBit(neighbour, i);
    neighbours.push_back(neighbour);
  }

  return neighbours;
}

double bits_to_number(std::vector<char> bits, bounds bounds,
                      int number_of_bits) {
  double result = 0.0;

  for (char bit : bits) {
    result = result * 2 + (bit - '0');
  }

  return bounds.min +
         result * (bounds.max - bounds.min) / (pow(2, number_of_bits) - 1);
}

std::vector<double> calculate_gradient(std::vector<char> bits, bounds bounds,
                                       int number_of_bits) {
  auto number = bits_to_number(bits, bounds, number_of_bits);

  return gradient_rastrigin({number});
}

double calculate_function(std::vector<char> bits, bounds bounds,
                          int number_of_bits) {
  auto number = bits_to_number(bits, bounds, number_of_bits);

  auto function_value = rastrigin({number});

  return function_value;
}

int main() {
  bounds rastrigin_bounds = {-5.12, 5.12};

  int number_of_bits =
      std::log2(pow(10, 10) * (rastrigin_bounds.max - rastrigin_bounds.min));

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