#include <iostream>
#include <vector>
#include <omp.h>
#include <climits>
#include <unordered_set>
#include <algorithm>
#include "queue.hpp"
#include "algorithms.hpp"



//This is a print function to help with the queue.print
void print_tour(const Tour& tour) {
    std::cout << "Tour: ";
    for (int i : tour.tour) {
        std::cout << i << " ";
    }
}

double Serial_compute_lbound(const std::vector<std::vector<double>>& distances, std::vector<std::vector<double>> &min, int f, int t, double LB) {
    double cf = 0;
    double ct = 0;
    double dist = distances[f][t];

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

double Parallel_first_lbound(const std::vector<std::vector<double>> &distances, std::vector<std::vector<double>> &min)
{
    double lowerbound = 0;

    #pragma omp parallel
    {
        double private_lb = 0;

        #pragma omp for
        for (int row = 0; row < distances.size(); row++) {
            double min1 = INT_MAX;
            double min2 = INT_MAX;

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
    }
    return lowerbound / 2;
}

Tour Serial_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors){
    int neighbor;
    double dist;
    //Lowerbound vectors for min1 and min2
    std::vector<std::vector<double>> min (N, std::vector<double>(2, INT_MAX));

    //Priority queue
    PriorityQueue<Tour> queue;
    
    //Tours include the path, cost and bound
    Tour tour, best_tour, new_tour;

    tour.bound = Serial_first_lbound(distances, min);//Tour lowerbound <- Initial lowerbound guess
    tour.tour.push_back(0); // Tour ← {0}
    tour.cost = 0; // Tour Cost <- 0
    queue.push(tour); //Queue ← (Tour, 0, LB, 1, 0) (Tour, Cost, Bound, Length, Current city)

    best_tour.tour.push_back(0);
    best_tour.cost = max_value;

    while (!queue.empty()){ //while Queue ̸= {} do
        tour = queue.pop(); // (Tour, Cost, Bound, Length, N ode) ← Queue.pop()

        if (tour.bound >= best_tour.cost){// if Bound ≥ BestT ourCost then
            return best_tour; // return BestT our, BestT ourCost
        }

        if (tour.tour.size() == N){ // if Length = N then
            dist = distances[tour.tour.back()][0];
            if (tour.cost + dist < best_tour.cost){ // if Cost + distances(Node, 0) < best_tourCost then
                best_tour.tour = tour.tour; // best_tour ← Tour ∪ {0}
                best_tour.tour.push_back(0); // best_tour ← Tour ∪ {0}
                best_tour.cost = tour.cost + dist; // BestT ourCost ← Cost + distances(N ode, 0)
            } 
        }
        else{
            //for each neighbor v of Node and v ̸∈ Tour do
            for (int v = 0; v < neighbors[tour.tour.back()].size(); v++){
                neighbor = neighbors[tour.tour.back()][v];

                //If it already belongs to tour quit
                if (std::find(tour.tour.begin(), tour.tour.end(), neighbor) != tour.tour.end()) {
                    continue;
                }

                //newBound ← updated lower bound on tour cost
                new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), neighbor, tour.bound);
                
                if (new_tour.bound > best_tour.cost){ //if newBound > BestT ourCost then
                    continue;
                }

                new_tour.tour = tour.tour;
                new_tour.tour.push_back(neighbor); // new_tour ← Tour ∪ {v}
                new_tour.cost = tour.cost + distances[tour.tour.back()][neighbor]; //newCost ← cost + distances(N ode, v)
                queue.push(new_tour); //Queue.add((new_tour, newCost, newBound, Length + 1, v)), v is the new current node
            }//end for
        }//end if
    }//end while
    return best_tour;
}//end procedure

Tour Parallel_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors){
    
    int index = 0;
    double min_cost;
    //---------------------------------Shared variables -----------------------------------

    std::vector<std::vector<double>> min (N, std::vector<double>(2));
    
    std::vector<PriorityQueue<Tour>> queues;

    Tour tour, best_tour, new_tour;
    std::vector<Tour> best_tours;

    //--------------------Setting tours information -----------------------------------------
    tour.tour.push_back(0); 
    tour.cost = 0; 
    tour.bound = Parallel_first_lbound(distances, min); 

    best_tour.tour.push_back(0);
    best_tour.cost = max_value;
    //--------------------------------------------------------------------------------------

    if (tour.bound >= best_tour.cost){
        return best_tour; 
    }

    #pragma omp parallel
    {
        int total_td = omp_get_num_threads();
        std::cout << "Total number of threads: " << total_td << std::endl;
        int neighbor;
        int total_tdm = omp_get_max_threads();
        std::cout << "Max number of threads: " << total_tdm << std::endl;

        #pragma omp for private(new_tour) 
        for (int v = 0; v < neighbors[tour.tour.back()].size(); v++){
            neighbor = neighbors[tour.tour.back()][v];
            
            PriorityQueue<Tour> queue;

            if (std::find(tour.tour.begin(), tour.tour.end(), neighbor) != tour.tour.end()) {
                continue;
            }

            new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), neighbor, tour.bound);
            
            if (new_tour.bound > max_value){ 
                continue;
            }

            new_tour.tour = tour.tour;
            new_tour.tour.push_back(neighbor); 
            new_tour.cost = tour.cost + distances[tour.tour.back()][neighbor];
            queue.push(new_tour);

            //#pragma omp atomic
            queues.push_back(queue);
        }

        #pragma omp for private(tour, new_tour) 
        for (int i = 0; i < queues.size(); i++){
            while (!queues[i].empty()){ 
                tour = queues[i].pop(); 
                double distance = distances[0][tour.tour.back()];

                if (tour.bound >= best_tour.cost){
                    break;
                }

                if (tour.tour.size() == N){
                    #pragma omp critical
                    {
                        if (tour.cost + distance < best_tour.cost){ 
                                best_tour.tour = tour.tour; 
                                best_tour.tour.push_back(0); 
                                best_tour.cost = tour.cost + distance;
                            }
                    } 
                }
                else{
                    for (int v = 0; v < neighbors[tour.tour.back()].size(); v++){
                        neighbor = neighbors[tour.tour.back()][v];

                        if (std::find(tour.tour.begin(), tour.tour.end(), neighbor) != tour.tour.end()) {
                            continue;
                        }

                        new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), neighbor, tour.bound);
                        
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

Tour Parallel2_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors){
    
       int index = 0;
    double min_cost;
    //---------------------------------Shared variables -----------------------------------

    std::vector<std::vector<double>> min (N, std::vector<double>(2));
    
    std::vector<PriorityQueue<Tour>> queues;

    Tour tour, best_tour, new_tour;
    std::vector<Tour> best_tours;

    //--------------------Setting tours information -----------------------------------------
    tour.tour.push_back(0); 
    tour.cost = 0; 
    tour.bound = Parallel_first_lbound(distances, min); 

    best_tour.tour.push_back(0);
    best_tour.cost = max_value;
    //--------------------------------------------------------------------------------------

    if (tour.bound >= best_tour.cost){
        return best_tour; 
    }


        int total_td = omp_get_num_threads();
        std::cout << "Total number of threads: " << total_td << std::endl;
        int neighbor;
        int total_tdm = omp_get_max_threads();
        std::cout << "Max number of threads: " << total_tdm << std::endl;

        #pragma omp parallel for private(new_tour) 
        for (int v = 0; v < neighbors[tour.tour.back()].size(); v++){
            neighbor = neighbors[tour.tour.back()][v];
            
            PriorityQueue<Tour> queue;

            if (std::find(tour.tour.begin(), tour.tour.end(), neighbor) != tour.tour.end()) {
                continue;
            }

            new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), neighbor, tour.bound);
            
            if (new_tour.bound > max_value){ 
                continue;
            }

            new_tour.tour = tour.tour;
            new_tour.tour.push_back(neighbor); 
            new_tour.cost = tour.cost + distances[tour.tour.back()][neighbor];
            queue.push(new_tour);

            //#pragma omp atomic
            queues.push_back(queue);
        }

        #pragma omp parallel for private(tour, new_tour) 
        for (int i = 0; i < queues.size(); i++){
            while (!queues[i].empty()){ 
                tour = queues[i].pop(); 
                double distance = distances[0][tour.tour.back()];

                if (tour.bound >= best_tour.cost){
                    break;
                }

                if (tour.tour.size() == N){
                    #pragma omp critical
                    {
                        if (tour.cost + distance < best_tour.cost){ 
                                best_tour.tour = tour.tour; 
                                best_tour.tour.push_back(0); 
                                best_tour.cost = tour.cost + distance;
                            }
                    } 
                }
                else{
                    for (int v = 0; v < neighbors[tour.tour.back()].size(); v++){
                        neighbor = neighbors[tour.tour.back()][v];

                        if (std::find(tour.tour.begin(), tour.tour.end(), neighbor) != tour.tour.end()) {
                            continue;
                        }

                        new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), neighbor, tour.bound);
                        
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
        
        return best_tour;
    
    }
