#include <iostream>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include <climits>
#include <algorithm>
#include "queue.hpp"
#include "algorithms.hpp"

#define TAG 100000

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

bool v_cmp(const Tour& a, const Tour& b) {
    return a > b;
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


Tour Serial_MPI_tsp_bb(const MPI_Comm comm, const int num_nodes, const int node_id, const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors, const int layer_cap){
    //---------------------------------MPI variables -----------------------------------
    
    double exec_time = -MPI_Wtime();

    MPI_Request request;
 
    double best_cost = max_value;
    int best_cost_flag = 1;

    int iterations = 0;    
    int claimed;
    int aux;
    //---------------------------------Private variables -----------------------------------

    int neighbor;
    double distance;
    
    //---------------------------------Shared variables -----------------------------------

    std::vector<std::vector<double>> min (N, std::vector<double>(2));
    
    std::vector<PriorityQueue<Tour>> queues;

    Tour tour, best_tour, new_tour;
    Tour aux_tour, aux2_tour;

    std::vector<std::vector<Tour>> tour_matrix(layer_cap + 1);

    //--------------------------------------------------------------------------------------

    tour.tour.push_back(0); 
    tour.cost = 0; 
    best_tour.tour.push_back(0);
    best_tour.cost = max_value;

    tour.bound = Serial_first_lbound(distances,min);

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
        
        tour_matrix[0].push_back(new_tour);
        
    }


    for (int layer = 0; layer < layer_cap; layer ++){
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
                
                tour_matrix[layer+1].push_back(new_tour);

            }
        }
    }
    
    std::sort(tour_matrix[tour_matrix.size()-1].begin(), tour_matrix[tour_matrix.size()-1].end(), v_cmp);

    for (int i = 0; i <  tour_matrix[tour_matrix.size()-1].size(); i++){

        if (i < 2){
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
                
                for (int vo = 0; vo < neighbors[neighbor].size(); vo++){
                
			int neighbor2 = neighbors[neighbor][vo];

			if (std::find(tour_matrix[tour_matrix.size()-1][i].tour.begin(), tour_matrix[tour_matrix.size()-1][i].tour.end(), neighbor2) != tour_matrix[tour_matrix.size()-1][i].tour.end()) 					 		{
			    continue;
			}
			
			int distance2 = distances[neighbor][neighbor2];
			aux_tour.bound = Serial_compute_lbound(distance2, min, neighbor, neighbor2, new_tour.bound);
			
			if (aux_tour.bound > max_value){ 
			    continue;
			}
                
                	aux_tour.tour = tour_matrix[tour_matrix.size()-1][i].tour;

			aux_tour.tour.push_back(neighbor);
			aux_tour.tour.push_back(neighbor2); 
			
			aux_tour.cost = tour_matrix[tour_matrix.size()-1][i].cost + distance + distance2;
			
			PriorityQueue<Tour> queue;
			queue.push(aux_tour);

			queues.push_back(queue);
        
                }
                
            }
        }
        else{
            new_tour = tour_matrix[tour_matrix.size()-1][i];
            PriorityQueue<Tour> queue;
            queue.push(new_tour);

            queues.push_back(queue);

        }

    }

    std::sort(queues.begin(), queues.end(), cmp);

    for (int i = node_id; i < queues.size(); i++){
    
        claimed = 0;
        
        for (int j = 0; j < num_nodes; j++) {
            if (j != node_id) {
                int probe_flag = 0;
                MPI_Iprobe(j, i, comm, &probe_flag, MPI_STATUS_IGNORE);
                if (probe_flag) {
                    MPI_Request rr;
                    MPI_Irecv(&aux, 1, MPI_INT, j, i, comm, &rr);
                    claimed = 1;
                }
                else{
                    MPI_Request r;
                    MPI_Isend(&i, 1, MPI_INT, j, i, comm, &r);
                }
            }
        }

        if (claimed) {
            continue;
        }
        
        
        for (int j = 0; j < num_nodes; j++) {
            if (j != node_id) {
                int probe_flag = 0;
                MPI_Iprobe(j, i, comm, &probe_flag, MPI_STATUS_IGNORE);
                if (probe_flag) {
                    MPI_Request rr;
                    MPI_Irecv(&aux, 1, MPI_INT, j, i, comm, &rr);
                    if (node_id > j){
                    	claimed = 1;
                    }
                }
            }
        }
        
        if (claimed) {
            continue;
        }
        

        iterations = 0;

        while (!queues[i].empty()){

            if (iterations == 0){
                if (best_cost_flag){
                    MPI_Iallreduce(&best_tour.cost, &best_cost, 1, MPI_DOUBLE, MPI_MIN, comm, &request);
                    best_cost_flag = 0;
                }

                MPI_Test(&request, &best_cost_flag, MPI_STATUS_IGNORE);

            }

            iterations++;
            
            if (iterations >= 100){
                iterations = 0;
            }
            
            tour = queues[i].pop(); 

            if (tour.bound >= best_cost){
                break;
            }

            if (tour.tour.size() == N){
                distance = distances[0][tour.tour.back()];
                if (tour.cost + distance < best_cost){
                    //Update local and global best cost
                    best_cost = tour.cost + distance;
                    best_tour.cost = best_cost;
                    best_tour.tour = tour.tour; 
                    best_tour.tour.push_back(0);
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
                    
                    if (new_tour.bound > best_cost){ 
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



