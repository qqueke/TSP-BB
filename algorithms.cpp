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

Tour Parallel_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors, const int layer_cap){
    //---------------------------------Private variables -----------------------------------
    int neighbor;
    double distance;
    
    //---------------------------------Shared variables -----------------------------------

    std::vector<std::vector<double>> min (N, std::vector<double>(2));
    
    std::vector<PriorityQueue<Tour>> queues;

    Tour tour, best_tour, new_tour;

    std::vector<std::vector<Tour>> tour_matrix(layer_cap + 1);

    double lowerbound;

    //--------------------------------------------------------------------------------------

    tour.tour.push_back(0); 
    tour.cost = 0; 
    best_tour.tour.push_back(0);
    best_tour.cost = max_value;

    #pragma omp parallel shared(lowerbound)
    {
        double private_lb = 0;   
        double min1;
        double min2;

        #pragma omp for schedule(static) nowait
        for (int row = 0; row < distances.size(); row++) {
            min1 = INT_MAX;
            min2 = INT_MAX;

            for (int column = 0; column < distances[row].size(); column++) {
                double distance = distances[row][column];

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

            private_lb += min1 + min2;
        }

        #pragma omp atomic
        lowerbound += private_lb;
        
        #pragma omp barrier
        #pragma omp single
        {
            tour.bound = lowerbound/2;
        
            if (tour.bound >= best_tour.cost){
                omp_set_num_threads(0); 
            }
        }

        #pragma omp for private(new_tour, distance, neighbor) schedule(static)
        for (int v = 0; v < neighbors[0].size(); v++){
            neighbor = neighbors[0][v];

            distance = distances[0][neighbor];
            new_tour.bound = Serial_compute_lbound(distance, min, 0, neighbor, tour.bound);
            
            if (new_tour.bound > max_value){ 
                continue;
            }

            new_tour.tour = tour.tour;
            new_tour.tour.push_back(neighbor); 
            new_tour.cost = tour.cost + distance;
            
            #pragma omp critical
            tour_matrix[0].push_back(new_tour);
            
        }


        for (int layer = 0; layer < layer_cap; layer ++){
            #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
            for (int i= 0; i < tour_matrix[layer].size(); i++){
                for (int v = 0; v < neighbors[tour_matrix[layer][i].tour.back()].size(); v++){
                    
                    neighbor = neighbors[tour_matrix[layer][i].tour.back()][v];

                    if (std::find(tour_matrix[layer][i].tour.begin(), tour_matrix[layer][i].tour.end(), neighbor) != tour_matrix[layer][i].tour.end()) {
                        continue;
                    }

                    distance = distances[tour_matrix[layer][i].tour.back()][neighbor];
                    new_tour.bound = Serial_compute_lbound(distance, min, tour_matrix[layer][i].tour.back(), neighbor, tour_matrix[layer][i].bound);
                    
                    if (new_tour.bound > max_value){ 
                        continue;
                    }

                    new_tour.tour = tour_matrix[layer][i].tour;
                    new_tour.tour.push_back(neighbor); 
                    new_tour.cost = tour_matrix[layer][i].cost + distance;
                    
                    #pragma omp critical
                    tour_matrix[layer+1].push_back(new_tour);

                }
            }
        }
        

        #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
        for (int i = 0; i <  tour_matrix[tour_matrix.size()-1].size(); i++){
            for (int v = 0; v < neighbors[tour_matrix[tour_matrix.size()-1][i].tour.back()].size(); v++){
                
                neighbor = neighbors[tour_matrix[tour_matrix.size()-1][i].tour.back()][v];

                if (std::find(tour_matrix[tour_matrix.size()-1][i].tour.begin(), tour_matrix[tour_matrix.size()-1][i].tour.end(), neighbor) != tour_matrix[tour_matrix.size()-1][i].tour.end()) {
                    continue;
                }

                distance = distances[tour_matrix[tour_matrix.size()-1][i].tour.back()][neighbor];
                new_tour.bound = Serial_compute_lbound(distance, min, tour_matrix[tour_matrix.size()-1][i].tour.back(), neighbor, tour_matrix[tour_matrix.size()-1][i].bound);
                
                if (new_tour.bound > max_value){ 
                    continue;
                }

                new_tour.tour = tour_matrix[tour_matrix.size()-1][i].tour;

                new_tour.tour.push_back(neighbor); 
                new_tour.cost = tour_matrix[tour_matrix.size()-1][i].cost + distance;
                
                PriorityQueue<Tour> queue;
                queue.push(new_tour);

                #pragma omp critical
                queues.push_back(queue);
            }
        }

        #pragma omp single
        std::sort(queues.begin(), queues.end(), cmp);

        #pragma omp for private(tour, new_tour, distance, neighbor) schedule(dynamic) nowait
        for (int i = 0; i < queues.size(); i++){
            while (!queues[i].empty()){ 
                tour = queues[i].pop(); 
                
                if (tour.bound >= best_tour.cost){
                    break;
                }

                if (tour.tour.size() == N){
                    distance = distances[0][tour.tour.back()];
                    #pragma omp critical
                    { 
                        if (tour.cost + distance < best_tour.cost){
                                best_tour.cost = tour.cost + distance;
                                best_tour.tour = tour.tour; 
                                best_tour.tour.push_back(0); 
                        }
                    }
                }
                else{
                    for (int v = 0; v < neighbors[tour.tour.back()].size(); v++){
                        neighbor = neighbors[tour.tour.back()][v];

                        if (std::find(tour.tour.begin(), tour.tour.end(), neighbor) != tour.tour.end()) {
                            continue;
                        }

                        distance = distances[tour.tour.back()][neighbor];
                        new_tour.bound = Serial_compute_lbound(distance, min, tour.tour.back(), neighbor, tour.bound);
                        
                        if (new_tour.bound > best_tour.cost){ 
                            continue;
                        }

                        new_tour.tour = tour.tour;
                        new_tour.tour.push_back(neighbor); 
                        new_tour.cost = tour.cost + distances[tour.tour.back()][neighbor];

                        queues[i].push(new_tour);
                    }  
                }
            }   
        } 
    }
    return best_tour;
}
