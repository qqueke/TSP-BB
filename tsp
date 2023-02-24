#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <vector>
#include <climits>
#include <algorithm>
#include <cfloat>

#include "/home/qqueke/PDC/PDC---Projeto/queue/queue.hpp"

typedef struct Tour{
    std::vector<int> tour;
    double cost;
    double bound;
    int length;
    int current_node;

    bool operator>(const Tour& other){
        return bound > other.bound;
    }
}Tour;


double Initial_LB(const std::vector<std::vector<double>>& Distances, std::vector<double>& min1, std::vector<double>& min2) {
    double lowerbound = 0;


    //Prov paralelizável com uma soma final de todas as partições
    for (int row = 0; row < Distances.size(); row++) {
        for (int column = 0; column < Distances.size(); column++) {
            if (Distances[row][column] < min1[row]) {
                min2[row] = min1[row];
                min1[row] = Distances[row][column];
            }
            else if (Distances[row][column] >= min1[row] && Distances[row][column] < min2[row]) {
                min2[row] = Distances[row][column];
            }
        }
        lowerbound += min1[row] + min2[row];
    }
    return lowerbound / 2;
}

double Compute_LB(const std::vector<std::vector<double>>& Distances, const std::vector<double>& min1, const std::vector<double>& min2, int f, int t, double LB) {
    double cf = 0;
    double ct = 0;

    if (Distances[f][t] >= min2[f]) {
        cf = min2[f];
    }
    else {
        cf = min1[f];
    }

    if (Distances[t][f] >= min2[t]){
        ct = min2[t];
    }
    else {
        ct = min1[t];
    }

    double lower_bound = LB + cf + ct;

    return lower_bound;
}

Tour TSPBB(const std::vector<std::vector<double>>& Distances, int N, double BestTourCost){
    
    //Estas variaveis ainda nem tem tipos certos,  simplesmente nao quis que estivessem
    //sublinhadas c os erros
    int row;
    std::vector<double> min1(Distances.size(), INT_MAX);
    std::vector<double> min2(Distances.size(), INT_MAX);
    PriorityQueue<Tour> queue;
    Tour tour, newTour, BestTour;

    newTour.tour.push_back(N+1);
    BestTour.tour.push_back(N+1);
    BestTour.cost = BestTourCost;

    tour.tour.push_back(N+1);
    tour.bound = Initial_LB(Distances, min1, min2);
    tour.cost = 0;
    tour.length = 1;
    tour.current_node = 0;


    while (!queue.empty()){
        tour = queue.pop();
        
        if (tour.bound >= BestTour.cost){
            return BestTour;
        }

        if (tour.length == N){
            if (tour.cost + Distances[tour.current_node][0] < BestTour.cost){
                BestTour.tour = tour.tour;
                BestTour.tour[tour.tour.size()] = 0;
                BestTour.cost = tour.cost + Distances[tour.current_node][0];
            }
        }
        else{
            //for each neighbor v of Node and v doesnt belong to Tour
            for (int v = 0; v < N; v++){
                //If v is neighbor and doesnt belong to tour then:
                //If it belongs to tour quit

                //Trocar por loop paralelizável
                if (std::find(tour.tour.begin(), tour.tour.end(), v) != tour.tour.end()) {
                    continue;
                }
                //If it is not a neighbot quit
                if (Distances[tour.current_node][v] == DBL_MAX){
                    continue;
                }

                newTour.bound = Compute_LB(Distances, min1, min2, tour.current_node, v, tour.bound);
                
                if (newTour.bound > BestTour.cost){
                    continue;
                }

                newTour.tour = tour.tour;
                newTour.tour[tour.tour.size()] = v;
                newTour.cost = tour.cost + Distances[tour.current_node][v];
                newTour.length = tour.length + 1;
                newTour.current_node = v;
                queue.push(newTour);
            }
        }
    }
    return BestTour;
}


int main(int argc, char *argv[]) {
    double exec_time;

    if (argc != 3) {
        printf("Usage: tsp <cities file> <max-value>\n");
        return 1;
    }

    char *cities_file = argv[1];
    int max_value = atoi(argv[2]);

    FILE *fp = fopen(cities_file, "r");
    if (fp == NULL) {
        printf("Error: Unable to open file %s\n", cities_file);
        return 1;
    }

    int num_cities, num_roads;
    fscanf(fp, "%d %d", &num_cities, &num_roads);

    std::vector<std::vector<double>> Distances (num_cities, std::vector<double>(num_cities));
    


    for (int row = 0; row < num_cities; row++) {
        for (int column = 0; column < num_cities; column++) {
            Distances[row][column] = DBL_MAX;
        }
    }

    for (int row = 0; row < num_roads; row++) {
        int city1, city2;
        double distance;
        fscanf(fp, "%d %d %lf", &city1, &city2, &distance);
        Distances[city1][city2] = distance;
        Distances[city2][city1] = distance;
    }

    fclose(fp);

    //exec_time = -omp_get_wtime();

    Tour BestTour = TSPBB(Distances, num_cities, max_value);
    
    //exec_time += omp_get_wtime();

    
    fprintf(stderr, "%.1fs\n", exec_time);

    //No solution that has a better value than the max admited
    if (BestTour.cost > max_value){
        printf("NO SOLUTION");
    }
    //This means the graph is disconnected
    else if (BestTour.tour.size() != num_cities +1){
        printf("NO SOLUTION");
    }
    //Valid solution
    else{
        printf("%.1lf\n", BestTour.cost);
        for (int row = 0; row < num_cities; row ++){
            printf("%d", BestTour.tour[row]);
        }
    }

    return 0;
}
