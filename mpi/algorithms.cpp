#include <iostream>
#include <vector>
#include <omp.h>
#include <mpi.h>
#include <climits>
#include <algorithm>
#include "queue.hpp"
#include "algorithms.hpp"

#define TAG 123
#define TAG 123
#define TAG 123
#define TAG 123

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

        #pragma omp barrier
        #pragma omp single
        std::cout << "Node " << node_id << " -> number of queues: " << queues.size() <<std::endl;

        #pragma omp single
        std::sort(queues.begin(), queues.end(), cmp);

        int flag = 1;
        int iterations = 0;

        #pragma omp for private(tour, new_tour, distance, neighbor) schedule(dynamic) nowait
        for (int i = node_id; i < queues.size(); i++){
        
            int claimed = 0;
            
            #pragma omp critical
            {
                // Notify other processes that the current iteration has been claimed

                // Check if any other process has claimed the current iteration
                for (int j = 0; j < num_nodes; j++) {
                    if (j != node_id) {
                        int probe_flag = 0;
                        MPI_Iprobe(j, i, comm, &probe_flag, MPI_STATUS_IGNORE);
                        if (probe_flag) {
                            claimed = 1;
                            break;
                        }
                    }
                }

                if (claimed != 1){
                    for (int j = 0; j < num_nodes; j++) {
                        if (j != node_id) {
                            MPI_Request r;
                            MPI_Isend(&i, 1, MPI_INT, j, i, comm, &r);
                        }
                    }
                }
            }

            if (claimed) {
                continue;
            }

            while (!queues[i].empty()){

                #pragma omp critical
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
                            //std::cout << "Node " << node_id << ": Reduction is complete and the final value is: "<< best_cost << std::endl;
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
    
    /*std::cout << "Node "<< node_id << " finished in " << exec_time << "s and path:"<< std::endl;
    std::cout << best_tour.cost << std::endl;

    for (int i = 0; i<N +1; i++) {
        std::cout << best_tour.tour[i] << " ";
    }
    std::cout << std::endl;*/
    
    MPI_Gather(&best_tour.cost, 1, MPI_DOUBLE, &costs[0], 1, MPI_DOUBLE, 0, comm);
    MPI_Gather(&best_tour.tour[0], N+1 , MPI_INT, &tours[0], N+1 , MPI_INT, 0, comm);


    if (node_id == 0){
        int minimum = INT_MAX;
        int index;
        for (int i = 0; i < costs.size(); i++){
            if (costs[i] < minimum){
                minimum = costs[i];
                index = i;
            }
        }

        /*int counter = 0;
        for (int i = 0; i < num_nodes*(N+1); i++){
            if (counter == 0){
                std::cout << "New path: " << std::endl;
                std::cout << tours[i] << " ";
                counter++;
            }
            else if (counter < N +1 ){
                std::cout << tours[i] << " ";
                counter++;
            }
            else{
                std::cout << "\nNew path: " << std::endl;
                std::cout << tours[i] << " ";
                counter = 1;
            }

        }

        std::cout << "\nIndex for the best path: " << index << std::endl;*/
        
        best_tour.cost = costs[index];
        best_tour.tour.clear();
        best_tour.tour.resize(N+1);
 
        int ind = 0;
        for (int i = index * (N+1); i < (index +1) * (N+1); i++){
            best_tour.tour[ind] = tours[i];
            ind++;
        }
    }

    return best_tour;
}

Tour Serial_MPI_tsp_bb(const MPI_Comm comm, const int num_nodes, const int node_id, const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors, const int layer_cap){
    //---------------------------------MPI variables -----------------------------------
    
    double exec_time = -MPI_Wtime();

    MPI_Request request;
    MPI_Status status;

    double best_cost = max_value;

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

    std::sort(queues.begin(), queues.end(), cmp);

    int flag = 1;
    int iterations = 0;

    for (int i = node_id; i < queues.size(); i++){
        
        int claimed = 0;
        
        for (int j = 0; j < num_nodes; j++) {
            if (j != node_id) {
                int probe_flag = 0;
                MPI_Iprobe(j, i, comm, &probe_flag, MPI_STATUS_IGNORE);
                if (probe_flag) {
                    claimed = 1;
                    break;
                }
            }
        }

        if (claimed != 1){
            for (int j = 0; j < num_nodes; j++) {
                if (j != node_id) {
                    MPI_Request r;
                    MPI_Isend(&i, 1, MPI_INT, j, i, comm, &r);
                }
            }
        }
        
        if (claimed) {
            continue;
        }

        while (!queues[i].empty()){

            iterations++;
            if (iterations > 5){
                if (flag){
                    MPI_Iallreduce(&best_tour.cost, &best_cost, 1, MPI_DOUBLE, MPI_MIN, comm, &request);
                    flag = 0;
                }

                MPI_Test(&request, &flag, &status);
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
    

    exec_time += MPI_Wtime();
    /*std::cout << "Node "<< node_id << " finished in " << exec_time << "s and path:"<< std::endl;
    std::cout << best_tour.cost << std::endl;

    for (int i = 0; i<N +1; i++) {
        std::cout << best_tour.tour[i] << " ";
    }
    std::cout << std::endl;*/
    
    MPI_Gather(&best_tour.cost, 1, MPI_DOUBLE, &costs[0], 1, MPI_DOUBLE, 0, comm);
    MPI_Gather(&best_tour.tour[0], N+1 , MPI_INT, &tours[0], N+1 , MPI_INT, 0, comm);


    if (node_id == 0){
        int minimum = INT_MAX;
        int index;
        for (int i = 0; i < costs.size(); i++){
            if (costs[i] < minimum){
                minimum = costs[i];
                index = i;
            }
        }

        /*int counter = 0;
        for (int i = 0; i < num_nodes*(N+1); i++){
            if (counter == 0){
                std::cout << "New path: " << std::endl;
                std::cout << tours[i] << " ";
                counter++;
            }
            else if (counter < N +1 ){
                std::cout << tours[i] << " ";
                counter++;
            }
            else{
                std::cout << "\nNew path: " << std::endl;
                std::cout << tours[i] << " ";
                counter = 1;
            }

        }

        std::cout << "\nIndex for the best path: " << index << std::endl;*/
        
        best_tour.cost = costs[index];
        best_tour.tour.clear();
        best_tour.tour.resize(N+1);
 
        int ind = 0;
        for (int i = index * (N+1); i < (index +1) * (N+1); i++){
            best_tour.tour[ind] = tours[i];
            ind++;
        }
    }

    return best_tour;
}