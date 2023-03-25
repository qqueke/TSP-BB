#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <limits.h>
#include <string.h>
#include <vector>
#include <climits>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <mpi.h>
#include "queue.hpp"
#include "algorithms.hpp"


int main(int argc, char *argv[]) {
    int node_id, num_nodes;
    int num_threads;
    int layer_cap;
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
    int slices = (1+ num_roads/num_cities)*2;


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

    char* num_threads_str = std::getenv("OMP_NUM_THREADS");
    if (num_threads_str != nullptr) {
        num_threads = std::atoi(num_threads_str);
    } else {
        num_threads = omp_get_max_threads();
    }

    if (num_threads <= 5){
        layer_cap = 1;
    }
    else if (num_threads < 14){
        layer_cap = 2;
    }
    else{
        layer_cap = 3;
    }

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &node_id);
    MPI_Comm_size (MPI_COMM_WORLD, &num_nodes);

    double exec_time = -omp_get_wtime();

    Tour best_tour = Parallel_MPI_tsp_bb(MPI_COMM_WORLD, num_nodes, node_id, distances, num_cities, max_value, neighbors, layer_cap);
    
    exec_time += omp_get_wtime();
    MPI_Finalize ()
    fprintf(stderr, "%lfs\n", exec_time);

    best_tour.tour.shrink_to_fit();
    //No solution that has a better value than the max admited
    if (best_tour.cost > max_value){
        std::cout << "NO SOLUTION";
    }
    //This means the graph is disconnected
    else if (best_tour.tour.size() != num_cities +1){
        printf("NO SOLUTION");
    }
    //Valid solution
    else{
        std::cout << best_tour.cost << std::endl;

        for (int i = 0; i<best_tour.tour.size(); i++) {
            std::cout << best_tour.tour[i] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}

