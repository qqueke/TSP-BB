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

    if (distances[f][t] >= min[f][1]) {
        cf = min[f][1];
    }
    else {
        cf = min[f][0];
    }

    if (distances[t][f] >= min[t][1]){
        ct = min[t][1];
    }
    else {
        ct = min[t][0];
    }

    double lower_bound = LB + distances[f][t] - (cf + ct)/2;

    return lower_bound;
}

double Serial_first_lbound(const std::vector<std::vector<double>> &distances, std::vector<std::vector<double>> &min)
{
    double lowerbound = 0;

    for (int row = 0; row < distances.size(); row++) {
        for (int column = 0; column < distances[row].size(); column++) {
            if (distances[row][column] < min[row][0]) {
                min[row][1] = min[row][0];
                min[row][0] = distances[row][column];
            }
            else if (distances[row][column] >= min[row][0] && distances[row][column] < min[row][1]) {
                min[row][1] = distances[row][column];
            }
        }
        lowerbound += min[row][0] + min[row][1];
    }
    return lowerbound / 2;
}

double Parallel_first_lbound(const std::vector<std::vector<double>> &distances, std::vector<std::vector<double>> &min)
{
    double lowerbound = 0;

    //Using distance, min1 and min2 we won't be accessing the matrix too many times
    #pragma omp parallel for collapse(2) reduction(+: lowerbound)
    for (int row = 0; row < distances.size(); row++) {
        double min1 = INT_MAX;
        double min2 = INT_MAX;
        //#pragma omp parallel for
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

        lowerbound += min[row][0] + min[row][1];
    }
    return lowerbound / 2;
}

Tour Serial_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value){

    //Lowerbound vectors for min1 and min2
    std::vector<std::vector<double>> min (N, std::vector<double>(2));

    for (int row = 0; row < N; row++) {
        for (int column = 0; column < 2; column++) {
            min[row][column] = INT_MAX;
        }
    }

    //Priority queue
    PriorityQueue<Tour> queue;
    

    //Tours include the path, cost and bound
    Tour tour, best_tour, new_tour;

    tour.tour.push_back(0); // Tour ← {0}
    tour.cost = 0; // Tour Cost <- 0
    tour.bound = Serial_first_lbound(distances, min); //Tour lowerbound <- Initial lowerbound guess
    queue.push(tour); //Queue ← (Tour, 0, LB, 1, 0) (Tour, Cost, Bound, Length, Current city)

    best_tour.tour.push_back(0);
    best_tour.cost = max_value;

    while (!queue.empty()){ //while Queue ̸= {} do
        tour = queue.pop(); // (Tour, Cost, Bound, Length, N ode) ← Queue.pop()

        if (tour.bound >= best_tour.cost){// if Bound ≥ BestT ourCost then
            return best_tour; // return BestT our, BestT ourCost
        }

        if (tour.tour.size() == N){ // if Length = N then
            if (tour.cost + distances[tour.tour.back()][0] < best_tour.cost){ // if Cost + distances(Node, 0) < best_tourCost then
                best_tour.tour = tour.tour; // best_tour ← Tour ∪ {0}
                best_tour.tour.push_back(0); // best_tour ← Tour ∪ {0}
                best_tour.cost = tour.cost + distances[tour.tour.back()][0];//BestT ourCost ← Cost + distances(N ode, 0)
            } 
        }
        else{
            //for each neighbor v of Node and v ̸∈ Tour do
            for (int v = 0; v < N; v++){
                //If v is neighbor and doesnt belong to tour then:
                //If it is not a neighbor quit
                if (distances[tour.tour.back()][v] == INT_MAX){
                    continue;
                }

                //This find method will be replaced
                //If it already belongs to tour quit
                if (std::find(tour.tour.begin(), tour.tour.end(), v) != tour.tour.end()) {
                    continue;
                }

                //newBound ← updated lower bound on tour cost
                new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), v, tour.bound);
                
                if (new_tour.bound > best_tour.cost){ //if newBound > BestT ourCost then
                    continue;
                }

                new_tour.tour = tour.tour;
                new_tour.tour.push_back(v); // new_tour ← Tour ∪ {v}
                new_tour.cost = tour.cost + distances[tour.tour.back()][v]; //newCost ← cost + distances(N ode, v)
                queue.push(new_tour); //Queue.add((new_tour, newCost, newBound, Length + 1, v)), v is the new current node
            }//end for
        }//end if
    }//end while
    return best_tour;
}//end procedure

Tour Parallel_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value){
    
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
        #pragma omp for private(new_tour) no wait
        for (int v = 0; v < N; v++){
            PriorityQueue<Tour> queue;

            //double distance = distances[tour.tour.back()][v];

            if (distances[tour.tour.back()][v] == INT_MAX){
                continue;
            }

            if (std::find(tour.tour.begin(), tour.tour.end(), v) != tour.tour.end()) {
                continue;
            }

            new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), v, tour.bound);
            
            if (new_tour.bound > max_value){ 
                continue;
            }

            new_tour.tour = tour.tour;
            new_tour.tour.push_back(v); 
            new_tour.cost = tour.cost + distances[tour.tour.back()][v];
            queue.push(new_tour);

            #pragma omp atomic
            queues.push_back(queue);
        }

            #pragma omp for private(tour, best_tour, new_tour) nowait
            for (int i = 0; i < queues.size(); i++){
                while (!queues[i].empty()){ 
                    tour = queues[i].pop(); 
                    //double distance = distances[tour.tour.back()][0];

                    if (tour.bound >= best_tour.cost){
                        #pragma omp atomic
                        best_tours.push_back(best_tour);
                        break;
                    }

                    if (tour.tour.size() == N){
                        if (tour.cost + distances[tour.tour.back()][0] < best_tour.cost){ 
                            best_tour.tour = tour.tour; 
                            best_tour.tour.push_back(0); 
                            best_tour.cost = tour.cost + distances[tour.tour.back()][0];
                        } 
                    }
                    else{
                        //#pragma omp for private(new_tour)
                        for (int v = 0; v < N; v++){

                            //double dist = distances[tour.tour.back()][v];

                            if (distances[tour.tour.back()][v] == INT_MAX){
                                continue;
                            }
                            //O(n)
                            if (std::find(tour.tour.begin(), tour.tour.end(), v) != tour.tour.end()) {
                                continue;
                            }

                            new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), v, tour.bound);
                            
                            if (new_tour.bound > best_tour.cost){ 
                                continue;
                            }

                            new_tour.tour = tour.tour;
                            new_tour.tour.push_back(v); 
                            new_tour.cost = tour.cost + distances[tour.tour.back()][v]; 
                            
                            queues[i].push(new_tour); 
                        }
                    }
                }
                #pragma omp atomic
                best_tours.push_back(best_tour);
            }
        

        //Here all threads should have their result computed to best_tours
        //So we're going to find the minimum value among all the results
        
        #pragma omp single copyprivate(min_cost)
        {
            min_cost = best_tours[index].cost;
        }

        #pragma omp for
        for (int i = 0; i < best_tours.size(); i++){
            #pragma omp critical
            {
                if (best_tours[i].cost < min_cost){
                    min_cost = best_tours[i].cost;
                    index = i;
                }
            }
        }
        return best_tours[index];
    }
}

Tour Parallel2_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value){

    /*---------------------------------Min 2D array-------------------------------------------*/
    std::vector<std::vector<double>> min (N, std::vector<double>(2));

    for (int row = 0; row < N; row++) {
        for (int column = 0; column < 2; column++) {
            min[row][column] = INT_MAX;
        }
    }
    /*----------------------------------Priority Queues-------------------------------------------*/
    PriorityQueue<Tour> queue;
    PriorityQueue<Tour> result_queue;
    /*-----------------------------------Tours------------------------------------------*/
    Tour tour, best_tour, new_tour;

    tour.tour.push_back(0); 
    tour.cost = 0; 
    tour.bound = Serial_first_lbound(distances, min); 
    queue.push(tour); 

    best_tour.tour.push_back(0);
    best_tour.cost = max_value;
    /*-----------------------------------------------------------------------------*/


    if (tour.bound >= best_tour.cost){
        return best_tour; 
    }
        
    #pragma omp parallel private(tour, best_tour, new_tour)
    {

        #pragma omp for nowait
        for (int v = 0; v < N; v++){

            //double distance = distances[tour.tour.back()][v];

            if (distances[tour.tour.back()][v] == INT_MAX){
                continue;
            }

            if (std::find(tour.tour.begin(), tour.tour.end(), v) != tour.tour.end()) {
                continue;
            }

            new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), v, tour.bound);
            
            if (new_tour.bound > max_value){ 
                continue;
            }

            new_tour.tour = tour.tour;
            new_tour.tour.push_back(v); 
            new_tour.cost = tour.cost + distances[tour.tour.back()][v];
            
            #pragma omp atomic
            queue.push(new_tour);
        }


        while (!queue.empty()){ 
            #pragma omp atomic
            tour = queue.pop(); 

            if (tour.bound >= best_tour.cost){
                #pragma omp atomic 
                result_queue.push(best_tour);
                break;
            }

            if (tour.tour.size() == N){ 
                if (tour.cost + distances[tour.tour.back()][0] < best_tour.cost){ 
                    best_tour.tour = tour.tour; 
                    best_tour.tour.push_back(0);
                    best_tour.cost = tour.cost + distances[tour.tour.back()][0];
                }
            }
            else{
                #pragma omp for nowait
                for (int v = 0; v < N; v++){

                    if (distances[tour.tour.back()][v] == INT_MAX){
                        continue;
                    }

                    if (std::find(tour.tour.begin(), tour.tour.end(), v) != tour.tour.end()) {
                        continue;
                    }
                    

                    new_tour.bound = Serial_compute_lbound(distances, min, tour.tour.back(), v, tour.bound);
                    
                    if (new_tour.bound > best_tour.cost){ 
                        continue;
                    }

                    new_tour.tour = tour.tour;
                    new_tour.tour.push_back(v); 
                    new_tour.cost = tour.cost + distances[tour.tour.back()][v];

                    #pragma omp atomic
                    queue.push(new_tour); 
                }
            }
        }
    }
    return result_queue.pop();
}