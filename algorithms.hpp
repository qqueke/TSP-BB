#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_set>
#include <vector>

typedef struct Tour{
    std::vector<int> tour;
    double cost;
    double bound;
    

    friend bool operator>(const Tour& self, const Tour& other){
        if (self.bound != other.bound) {
            return self.bound > other.bound;
        }
        else {
            return self.tour.back() > other.tour.back();
        }
    }
}Tour;

void print_tour(const Tour& tour);

double Serial_compute_lbound(const std::vector<std::vector<double>> &distances, const std::vector<std::vector<double>> &min, int f, int t, double LB);

double Serial_first_lbound(const std::vector<std::vector<double>> &distances, std::vector<std::vector<double>> &min);

double Parallel_first_lbound(const std::vector<std::vector<double>> &distances, std::vector<std::vector<double>> &min);

Tour Serial_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double BestTourCost);

Tour Parallel_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double BestTourCost);

Tour Parallel2_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double BestTourCost);