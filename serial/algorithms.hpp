#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
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

double Serial_compute_lbound(const double distance, const std::vector<std::vector<double>> &min, int f, int t, double LB);

double Serial_first_lbound(const std::vector<std::vector<double>> &distances, std::vector<std::vector<double>> &min);

Tour Serial_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors);

