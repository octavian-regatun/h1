#include <bitset>
#include <cstdio>
#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#define M_PI 3.14159265358979323846
#define DIMENSIONS 30
#define ITERATIONS 10000
#define PRECISION 5
#define MAX_IMPROVEMENT_ATTEMPTS 100  

struct bounds {
    double min, max;
};

int length;

double rastrigin(const std::vector<double>& x, int dimensions) {
    double sum = 10 * dimensions;
    for (int i = 0; i < dimensions; i++)
        sum += x[i] * x[i] - 10 * cos(2 * M_PI * x[i]);

    return sum;
}

#define M 10

double michalewicz(const std::vector<double>& x, int dimensions) {
    double sum = 0;
    for (int i = 0; i < dimensions; i++)
        sum -= sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / M_PI), 2 * M);

    return sum;
}

double dejong1(const std::vector<double>& x, int dimensions)
{
    double sum = 0.0;
    for (int i = 0; i < dimensions; i++)
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

int calculate_length(const bounds& bounds, int dimensions)
{
    return static_cast<int>(std::ceil(dimensions * std::log2(pow(10, PRECISION) * (bounds.max - bounds.min))));
}

std::random_device random_device;
std::mt19937 random_generator(random_device());

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

void flipBit(std::vector<char>& bits, int index) {
    auto& bit = bits.at(index);

    if (bit == '0')
        bit = '1';
    else
        bit = '0';
}

std::vector<std::vector<char>> generate_neighbours(const std::vector<char>& bits) {
    std::vector<std::vector<char>> neighbours;

    for (int i = 0; i < length; i++) {
        std::vector<char> neighbour(bits.begin(), bits.end());
        flipBit(neighbour, i);
        neighbours.push_back(neighbour);
    }

    return neighbours;
}
std::vector<double> bits_to_number(const std::vector<char>& bits, const bounds& bounds, int length, int dimensions) {
    std::vector<double> result;
    int bits_per_dimension = length / dimensions;

    for (int d = 0; d < dimensions; d++) {
        double sum = 0.0;
        for (int i = 0; i < bits_per_dimension; i++) {
            sum = sum * 2 + (bits[d * bits_per_dimension + i] - '0');
        }

        double real_val = bounds.min + sum * (bounds.max - bounds.min) / (pow(2, bits_per_dimension) - 1);
        result.push_back(real_val);
    }

    return result;
}


double calculate_function(const std::vector<char>& bits, double (*func)(const std::vector<double>&, int), const bounds& bounds, int dimensions) {
    std::vector<double> values = bits_to_number(bits, bounds, length, dimensions);
    return func(values, dimensions);
}


std::vector<char> improve(const std::vector<char>& vc, double (*func)(const std::vector<double>&, int), const bounds& bounds, int dimensions) {
    std::vector<char> best_neighbour = vc;
    double best_score = calculate_function(best_neighbour, func, bounds, dimensions);
    int bits_per_dimension = length / dimensions;

    for (int i = 0; i < bits_per_dimension; i++) { // For each bit position in a dimension
        std::vector<char> neighbour = vc;
        for (int d = 0; d < dimensions; d++) { // Flip the same bit for each dimension
            flipBit(neighbour, d * bits_per_dimension + i);
        }

        double score = calculate_function(neighbour, func, bounds, dimensions);
        if (score < best_score) {
            best_neighbour = neighbour;
            best_score = score;
        }
    }

    return best_neighbour;
}




std::vector<char> hill_climbing(const std::vector<char>& x, double (*func)(const std::vector<double>&, int), const bounds& bounds, int dimensions) {
    std::vector<char> best = x;
    for (int t = 0; t < ITERATIONS; t++) {
        std::vector<char> current = generate_random_bits(length);

        bool local = false;
        int improvement_attempts = 0;

        while (!local && improvement_attempts < MAX_IMPROVEMENT_ATTEMPTS) {
            std::vector<char> next = improve(current, func, bounds, dimensions);
            if (calculate_function(next, func, bounds, dimensions) <= calculate_function(current, func, bounds, dimensions)) {
                current = next;
            }
            else {
                local = true;
            }
            improvement_attempts++;
        }

        if (calculate_function(current, func, bounds, dimensions) < calculate_function(best, func, bounds, dimensions)) {
            best = current;
        }
    }

    return best;
}


int main() {
    bounds rastrigin_bounds = { -5.12, 5.12 };
    length = calculate_length(rastrigin_bounds, DIMENSIONS);

    auto bits = generate_random_bits(length);

    auto result_bits = hill_climbing(bits, rastrigin, rastrigin_bounds, 5);
    double result = calculate_function(result_bits, rastrigin, rastrigin_bounds,5);

    printf("Result: %f\n", result);

    return 0;
}

