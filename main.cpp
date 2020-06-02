#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <omp.h>
#include <algorithm>

#define POPULATION_SIZE 2000
#define SIZE 100

/*
struct Timer {
public:
    std::chrono::high_resolution_clock::time_point start_time;

    void start_timer() {
        start_time = std::chrono::steady_clock::now();
    }

    double get_time_elapsed() {
        std::chrono::high_resolution_clock::time_point now = std::chrono::steady_clock::now();;
        std::chrono::duration<double> time_span =  now - start_time;
        return time_span.count();
    }

};
*/
/* --- Calculation functions --- */

double target_function(double x) {
    return std::sin(x);
}

double get_random(double start, double end, std::mt19937 &rand) {
    std::uniform_real_distribution<> dist(start, end);

    return dist(rand);
}

unsigned int get_random_int(double end, std::mt19937 &rand) {
    std::uniform_int_distribution<> dist(0, end);

    return dist(rand);
}

struct XSort {
    bool operator() (const double &l, const double &r) const {
        return target_function(r) < target_function(l);
    }
};

/* --- Genetic algorithms functions --- */

void mutate(std::vector<double> &population, std::mt19937 &rand) {
    unsigned int population_size = population.size();

    for (unsigned int i = 0; i < (population_size * 2 / 10); i++) {
        unsigned int person = get_random_int(population_size, rand);

        population[person] = population[person] + get_random(-0.125, 0.125, rand);
    }
}

void spawn(std::vector<double> &population, std::mt19937 &rand) {
    unsigned int population_size = population.size();

    for (unsigned int i = 0; i < population_size; i++) {
        double father = population[i];
        double mother = population[get_random_int(population_size, rand)];

        population.push_back((father + mother) / 2);
    }

    while (population.size() < 500) {
        double father = population[get_random_int(population_size, rand)];
        double mother = population[get_random_int(population_size, rand)];

        population.push_back((father + mother) / 2);
    }
}

void kill(std::vector<double> &population) {
    population.erase(population.begin() + population.size() / 2, population.end());
}

double get_effectiveness(std::vector<double> &population) {
    double sum = 0;

    for (double i : population) {
        sum = target_function(i);
    }

    return sum / population.size();
}

/* --- Region calculation functions --- */

double calculate_region(double start, double end, std::mt19937 &rand) {
    std::vector<double> population(POPULATION_SIZE);

    for (int i = 0; i < POPULATION_SIZE; i++) {
        population.push_back(get_random(start, end, rand));
    }

    double prev_effectiveness = get_effectiveness(population);

    while (true) {
        double effectiveness = get_effectiveness(population);

        if (std::abs(effectiveness - prev_effectiveness) < 0.0025) {
            // return here best of population
            std::sort(population.begin(), population.end(), XSort());

            return population[0];
        }
        prev_effectiveness = effectiveness;

        // sort here
        std::sort(population.begin(), population.end(), XSort());
        kill(population);

        mutate(population, rand);
        spawn(population, rand);
    }
 }

int main() {

    //std::cout << "Max threads: " << omp_get_max_threads() << "\n";
    std::mt19937 rand { std::random_device{}() };
   // Timer timer = Timer();
    std::vector<double> picks;
   // timer.start_timer();
    
    {
	#pragma omp parallel for   
    	for (int i = 0; i < SIZE; i++) {
        	double pick = calculate_region(i * 2 * 3.14, i * 2 * 3.14 + 2 * 3.14, rand);
                #pragma omp critical (insert_picks)
		{
        		picks.push_back(pick);
		}
    	}
    }
   // double elapsed_time = timer.get_time_elapsed();

    for (int i = 0; i < SIZE; i++) {
        std::cout << "Pick #" << i << " x: " << picks[i] << "\n";
    }
   // std::cout << "Elapsed time: " << elapsed_time << "s\n";

    return 0;
}
