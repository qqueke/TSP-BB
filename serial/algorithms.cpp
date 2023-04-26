#include <iostream>
#include <vector>
#include <omp.h>
#include <climits>
#include <algorithm>
#include "queue.hpp"
#include "algorithms.hpp"

bool cmp(const PriorityQueue<Tour>& a, const PriorityQueue<Tour>& b) {
    
    double a_bound =  a.top().bound;
    double b_bound =  b.top().bound;
    if (a_bound != b_bound) {
        return a_bound < b_bound;
    }
    else {

        return a.top().tour.back() < b.top().tour.back();
    }
}

double Serial_compute_lbound(const double distance, std::vector<std::vector<double>> &min, int f, int t, double LB) {
    double cf = 0;
    double ct = 0;
    double dist = distance;

    if (dist >= min[f][1]) {
        cf = min[f][1];
    }
    else {
        cf = min[f][0];
    }

    if (dist >= min[t][1]){
        ct = min[t][1];
    }
    else {
        ct = min[t][0];
    }

    double lower_bound = LB + dist - (cf + ct)/2;

    return lower_bound;
}

double Serial_first_lbound(const std::vector<std::vector<double>> &distances, std::vector<std::vector<double>> &min)
{
    double distance;
    double lowerbound = 0;

    for (int row = 0; row < distances.size(); row++) {
            double min1 = INT_MAX;
            double min2 = INT_MAX;
        for (int column = 0; column < distances[row].size(); column++) {
            distance = distances[row][column];
            if (distance < min1) {
                min2 = min1;
                min1 = distance;
            }

            else if (distance >= min1 && distance < min2) {
                min2 = distance;
            }
        }

        min[row][0] = min1;
        min[row][1] = min2;

        lowerbound += min1 + min2;
    }
    return lowerbound / 2;
}

Tour Serial_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors){
    int neighbor;
    double dist;

    std::vector<std::vector<double>> min (N, std::vector<double>(2, INT_MAX));

    PriorityQueue<Tour> queue;
    
    Tour tour, best_tour, new_tour;

    tour.bound = Serial_first_lbound(distances, min);
    tour.tour.push_back(0); 
    tour.cost = 0; 
    queue.push(tour); 

    best_tour.tour.push_back(0);
    best_tour.cost = max_value;

    while (!queue.empty()){ 
        tour = queue.pop(); 

        if (tour.bound >= best_tour.cost){
            return best_tour; 
        }

        if (tour.tour.size() == N){ 
            dist = distances[tour.tour.back()][0];
            if (tour.cost + dist < best_tour.cost){ 
                best_tour.tour = tour.tour; 
                best_tour.tour.push_back(0); 
                best_tour.cost = tour.cost + dist; 
            } 
        }
        else{
            for (int v = 0; v < neighbors[tour.tour.back()].size(); v++){
                neighbor = neighbors[tour.tour.back()][v];
                dist = distances[tour.tour.back()][neighbor];
                
                if (std::find(tour.tour.begin(), tour.tour.end(), neighbor) != tour.tour.end()) {
                    continue;
                }

                new_tour.bound = Serial_compute_lbound(dist, min, tour.tour.back(), neighbor, tour.bound);
                
                if (new_tour.bound > best_tour.cost){ 
                    continue;
                }

                new_tour.tour = tour.tour;
                new_tour.tour.push_back(neighbor); 
                new_tour.cost = tour.cost + dist; 
                queue.push(new_tour); 
            }
        }
    }
    return best_tour;
}

