#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <omp.h>
#include <climits>
#include <algorithm>
#include "queue.hpp"
#include "algorithms.hpp"

int main(int argc, char *argv[]) {
    double exec_time;

    if (argc != 3) {
        std::cout << "Usage: tsp <cities file> <max-value>\n";
        return 1;
    }

    char *cities_file = argv[1];
    int max_value = atoi(argv[2]);

    FILE *fp = fopen(cities_file, "r");
    if (fp == NULL) {
        std::cout << "Error: Unable to open file " << cities_file << std::endl;
        return 1;
    }

    int num_cities, num_roads;
    if (!fscanf(fp, "%d %d", &num_cities, &num_roads)){
        std::cout << "Error";
        return 1;
    }

    std::vector<std::vector<double>> distances (num_cities, std::vector<double>(num_cities));
    std::vector<std::vector<int>> neighbors(num_cities);

    for (int row = 0; row < num_cities; row++) {
        for (int column = 0; column < num_cities; column++) {
            distances[row][column] = INT_MAX;
        }
    }

    for (int row = 0; row < num_roads; row++) {
        int city1, city2;
        double distance;
        if (!fscanf(fp, "%d %d %lf", &city1, &city2, &distance)){
            std::cout << "Error";
            return 1;
        }
        distances[city1][city2] = distance;
        distances[city2][city1] = distance;
        neighbors[city1].push_back(city2);
        neighbors[city2].push_back(city1);
    }

    fclose(fp);

    exec_time = -omp_get_wtime();

    Tour best_tour = Serial_tsp_bb(distances, num_cities, max_value, neighbors);
    
    exec_time += omp_get_wtime();

    fprintf(stderr, "%lfs\n", exec_time);

    //No solution that has a better value than the max admited
    if (best_tour.cost > max_value){
        std::cout << "NO SOLUTION";
    }
    //This means the graph is disconnected
    else if (best_tour.tour.size() != num_cities +1){
       std::cout << "NO SOLUTION";
    }
    //Valid solution
    else{
        std::cout << best_tour.cost << std::endl;

        for (int i : best_tour.tour) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
