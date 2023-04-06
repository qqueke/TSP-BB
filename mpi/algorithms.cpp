#include <iostream>
#include <vector>
#include <omp.h>
#include <mpi.h>
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

Tour Parallel_MPI_tsp_bb(const MPI_Comm comm, const int num_nodes, const int node_id, const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors, const int layer_cap){
    //---------------------------------MPI variables -----------------------------------
    
    double exec_time = -MPI_Wtime();

    MPI_Request recv_req[2];
    MPI_Request request;
    MPI_Request send_req[4];
    MPI_Status status;

    double best_cost = max_value;
    double mpi_cost=max_value;
    std::vector<int> mpi_tour(N+1, 0);

    std::vector<double> costs(num_nodes);    
    std::vector<int> tours(num_nodes*(N+1));

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

    #pragma omp parallel private(status, mpi_cost, mpi_tour) shared(comm)
    {
        double private_lb = 0;   
        double min1;
        double min2;
        int num_threads = omp_get_num_threads();

        #pragma omp single
        std::cout << "Node " << node_id << ", actual number of threads: " << num_threads << std::endl;

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

        int flag = 1;
        int iterations = 0;

        #pragma omp for private(tour, new_tour, distance, neighbor) schedule(dynamic)
        for (int i = node_id*num_threads; i < queues.size(); i+=num_nodes*num_threads){
            //std::cout << "Iteration: " << i << std::endl;
            while (!queues[i].empty()){

                #pragma omp critical
                //if (omp_get_thread_num() == 0)
                {
                    iterations++;
                    if (iterations > 5){
                        if (flag){
                            //std::cout << "Node " << node_id << ": is calling for a new reduction. "<< std::endl;
                            MPI_Iallreduce(&best_tour.cost, &best_cost, 1, MPI_DOUBLE, MPI_MIN, comm, &request);
                            flag = 0;
                        }

                        MPI_Test(&request, &flag, &status);
                        //if (flag){
                          //  std::cout << "Node " << node_id << ": Reduction is complete and the final value is: "<< best_cost << std::endl;
                        //}
                        iterations = 0;
                    }
                }
                
                tour = queues[i].pop(); 
                
                if (tour.bound >= best_cost){
                    break;
                }

                if (tour.tour.size() == N){
                    distance = distances[0][tour.tour.back()];
                    #pragma omp critical
                    { 
                        if (tour.cost + distance < best_cost){
                                //Update local and global best cost
                                best_cost = tour.cost + distance;
                                best_tour.cost = best_cost;
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
    }

    exec_time += MPI_Wtime();
    
    std::cout << "Node "<< node_id << " finished in " << exec_time << "s and cost: "<< best_tour.cost <<std::endl;

    return best_tour;
}

Tour Serial_MPI_tsp_bb(const MPI_Comm comm, const int num_nodes, const int node_id, const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors, const int layer_cap){
    //---------------------------------MPI variables -----------------------------------
    
    double exec_time = -MPI_Wtime();

    MPI_Request request[2];
    MPI_Request req;

    double best_cost = max_value;

    double my_iteration = node_id;
    int global_last_iter;

    int best_cost_flag = 1;
    int iterations_flag = 1;

    std::vector<double> costs(num_nodes);    
    std::vector<int> tours(num_nodes*(N+1));

    //---------------------------------Private variables -----------------------------------

    int neighbor;
    double distance;
    
    //---------------------------------Shared variables -----------------------------------

    std::vector<std::vector<double>> min (N, std::vector<double>(2));
    
    std::vector<PriorityQueue<Tour>> queues;

    Tour tour, best_tour, new_tour;

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

            queues.push_back(queue);
        }

    }
    int aux;

    std::sort(queues.begin(), queues.end(), cmp);



    int iterations = 0;
    
    
    
    int claimed;
    //std::vector<int> its(queues.size(), 0);

    for (int i = node_id; i < queues.size(); i++){

        claimed = 0;
        
        for (int j = 0; j < num_nodes; j++) {
            if (j != node_id) {
                int probe_flag = 0;
                
                //MPI_Request rr;
                //MPI_Irecv(&aux, 1, MPI_INT, j, i, comm, &rr);
                //MPI_Test(&rr, &probe_flag, MPI_STATUS_IGNORE);
                MPI_Iprobe(j, i, comm, &probe_flag, MPI_STATUS_IGNORE);
                if (probe_flag) {
                    MPI_Request rr;
                    MPI_Irecv(&aux, 1, MPI_INT, j, i, comm, &rr);
                    claimed = 1;
                    break;
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

        //its[i] = 1;
        iterations=0;
        //std::cout << "Node "<< node_id << " claimed " << i <<std::endl;
        while (!queues[i].empty()){

            /*if (iterations == 0){
                if (iterations_flag){
                    MPI_Iallreduce(MPI_IN_PLACE, &its, its.size(), MPI_INT, MPI_MAX, comm, &request[1]);
                    iterations_flag = 0;
                }

                MPI_Test(&request[1], &iterations_flag, MPI_STATUS_IGNORE);
            }*/



            if (iterations == 0){
                if (best_cost_flag){
                    MPI_Iallreduce(&best_tour.cost, &best_cost, 1, MPI_DOUBLE, MPI_MIN, comm, &request[0]);
                    best_cost_flag = 0;
                }



                MPI_Test(&request[0], &best_cost_flag, MPI_STATUS_IGNORE);
                //if (best_cost_flag ){
                    //std::cout << "Best cost updated to: " << best_cost << std::endl;
                //}
                if (iterations_flag){
                    int my_iteration = i;
                    MPI_Iallreduce(&my_iteration, &global_last_iter, 1, MPI_INT, MPI_MAX, comm, &request[1]);
                    iterations_flag = 0;
                }

                MPI_Test(&request[1], &iterations_flag, MPI_STATUS_IGNORE);
                if (iterations_flag ){
                    std::cout << "Last iter updated to: " << global_last_iter << std::endl;
                }

                //MPI_Iallreduce(&my_iteration, &last_iteration, 1, MPI_INT, MPI_MAX, comm, &request[1]);
                //MPI_Test(&request[1], &iterations_flag, MPI_STATUS_IGNORE);
                iterations=0;
            }

            iterations++;
            
            if (iterations >= 5){
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
        /*my_iteration = i;

        if (iterations_flag){
            my_iteration = -i;
            std::cout << "My iter: " << my_iteration << std::endl;
            MPI_Iallreduce(&my_iteration, &global_last_iter, 1, MPI_INT, MPI_MAX, comm, &req);
            iterations_flag = 0;
        }

        MPI_Test(&req, &iterations_flag, MPI_STATUS_IGNORE);
        if (iterations_flag){
            std::cout << "last iteration updated to: " << -global_last_iter << std::endl;

        }*/

    //std::cout << "Node "<< node_id << " jumping to " << last_iteration <<std::endl;
        //i = last_iteration;

    }

    exec_time += MPI_Wtime();
    std::cout << "Node "<< node_id << " finished in " << exec_time << "s and path cost:"<< best_tour.cost <<std::endl;

    //MPI_Allreduce(MPI_IN_PLACE, &its, queues.size(), MPI_INT, MPI_MAX, comm);

    /*if (node_id == 0){
        std::cout << "Visited queues: " << std::endl;
        for (int i = 0; i < its.size(); i ++)
        {
            std::cout << "iterations[" << i << "] = " << its[i] << std::endl;
        }
    }*/
    
    return best_tour;
}


