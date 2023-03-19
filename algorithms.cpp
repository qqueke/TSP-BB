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

bool cmp(const PriorityQueue<Tour>& a, const PriorityQueue<Tour>& b) {
    
    double a_bound =  a.top().bound;
    double b_bound =  b.top().bound;

    return a_bound < b_bound;
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

double Parallel_first_lbound(const std::vector<std::vector<double>> &distances, std::vector<std::vector<double>> &min)
{
    double lowerbound = 0;

    #pragma omp parallel
    {
        double private_lb = 0;
        double min1;
        double min2;


        #pragma omp for  private(min1, min2) schedule(dynamic)
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
                dist = distances[tour.tour.back()][neighbor];
                //If it already belongs to tour quit
                if (std::find(tour.tour.begin(), tour.tour.end(), neighbor) != tour.tour.end()) {
                    continue;
                }

                //newBound ← updated lower bound on tour cost
                new_tour.bound = Serial_compute_lbound(dist, min, tour.tour.back(), neighbor, tour.bound);
                
                if (new_tour.bound > best_tour.cost){ //if newBound > BestT ourCost then
                    continue;
                }

                new_tour.tour = tour.tour;
                new_tour.tour.push_back(neighbor); // new_tour ← Tour ∪ {v}
                new_tour.cost = tour.cost + dist; //newCost ← cost + distances(N ode, v)
                queue.push(new_tour); //Queue.add((new_tour, newCost, newBound, Length + 1, v)), v is the new current node
            }//end for
        }//end if
    }//end while
    return best_tour;
}//end procedure

Tour Parallel_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors){
        //Dynamic schedule alone improved alot 90s to 70s
    //Chunk size did not help
    int neighbor;
    double distance;
    int td_av = 0;
    //---------------------------------Shared variables -----------------------------------

    std::vector<std::vector<double>> min (N, std::vector<double>(2));
    
    std::vector<PriorityQueue<Tour>> queues;

    Tour tour, best_tour, new_tour;

    //--------------------------------------------------------------------------------------

    double lowerbound;

    #pragma omp parallel
    {

        tour.bound = Parallel_first_lbound(distances, min);

        double private_lb = 0;
        double min1;
        double min2;

        #pragma omp for  private(min1, min2) schedule(dynamic) nowait
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
        
        #pragma omp single
        {
            //std::cout << "Lowerbound: " << lowerbound << std::endl;
            tour.bound = lowerbound/2;
            tour.tour.push_back(0); 
            tour.cost = 0; 
            best_tour.tour.push_back(0);
            best_tour.cost = max_value;

            if (tour.bound >= best_tour.cost){
                omp_set_num_threads(0); 
            }
        }

        
        #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
        for (int v = 0; v < neighbors[tour.tour.back()].size(); v++){
            neighbor = neighbors[tour.tour.back()][v];

            if (std::find(tour.tour.begin(), tour.tour.end(), neighbor) != tour.tour.end()) {
                continue;
            }

            distance = distances[tour.tour.back()][neighbor];
            new_tour.bound = Serial_compute_lbound(distance, min, tour.tour.back(), neighbor, tour.bound);
            
            if (new_tour.bound > max_value){ 
                continue;
            }

            PriorityQueue<Tour> queue;
            new_tour.tour = tour.tour;
            new_tour.tour.push_back(neighbor); 
            new_tour.cost = tour.cost + distance;
            queue.push(new_tour);

            #pragma omp critical
            queues.push_back(queue);
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

Tour Parallel2_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors, int slices){
        //Dynamic schedule alone improved alot 90s to 70s
    //Chunk size did not help
    int neighbor;
    double distance;
    int td_av = 0;
    //---------------------------------Shared variables -----------------------------------

    std::vector<std::vector<double>> min (N, std::vector<double>(2));
    
    std::vector<PriorityQueue<Tour>> queues;

    Tour tour, best_tour, new_tour;

//neighbors[0].size()
    //--------------------------------------------------------------------------------------
    std::vector<Tour> tour_vector_2;
    std::vector<Tour> tour_vector_3;
    std::vector<Tour> tour_vector_4;
    std::vector<std::vector<Tour>> tour_matrix;

    int counter=0;
    int counter2=0;

    double lowerbound;
    //slices = 13;
    std::cout << "Number of slices: " << slices << std::endl;


    #pragma omp parallel
    {
        //tour.bound = Parallel_first_lbound(distances, min);

        double private_lb = 0;
        double min1;
        double min2;

        #pragma omp for  private(min1, min2) schedule(dynamic) nowait
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
        
        #pragma omp single
        {
            std::cout << "Lb: " << lowerbound << std::endl;
            tour.bound = lowerbound/2;
            tour.tour.push_back(0); 
            tour.cost = 0; 
            best_tour.tour.push_back(0);
            best_tour.cost = max_value;

            if (tour.bound >= best_tour.cost){
                omp_set_num_threads(0); 
            }
        }

        #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
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
            {
                tour_matrix[0].push_back(new_tour);
            }
        }


        for (int layer = 0; layer < 2; layer ++){
            tour_matrix.push_back(std::vector<Tour>());
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
                    {
                        tour_matrix[layer+1].push_back(new_tour);
                    }
                }
            }
        }

        #pragma omp single
        std::cout << "Number of queues: " << tour_vector_2.size() << std::endl;

        #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
        for (int i = 0; i < tour_matrix[tour_matrix.size()-1].size(); i++){
            std::vector<PriorityQueue<Tour>> queue_slices(slices);

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
                


                for (int partition = 0; partition < slices; partition++){
                    int part = (partition+1) * neighbors[tour_matrix[tour_matrix.size()-1][i].tour.back()].size()/slices;
                    if (v <= part){
                        queue_slices[partition].push(new_tour);
                        break;
                    }
                }
            }

            for (int partition = 0; partition < slices; partition++){
                if (!queue_slices[partition].empty()){
                    #pragma omp critical
                    queues.push_back(queue_slices[partition]);
                }
            }
        }
        
        #pragma omp single
        std::cout << "New number of queues: " << queues.size() << std::endl;

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

Tour Parallel3_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double max_value, const std::vector<std::vector<int>> &neighbors, int slices){    
        //Dynamic schedule alone improved alot 90s to 70s
    //Chunk size did not help
    int neighbor;
    double distance;
    int td_av = 0;
    //---------------------------------Shared variables -----------------------------------

    std::vector<std::vector<double>> min (N, std::vector<double>(2));
    
    std::vector<PriorityQueue<Tour>> queues;

    Tour tour, best_tour, new_tour;

//neighbors[0].size()
    //--------------------------------------------------------------------------------------
    std::vector<Tour> tour_vector_2;
    std::vector<Tour> tour_vector_3;
    std::vector<Tour> tour_vector_4;
    std::vector<Tour> tour_vector_5;

    int counter=0;
    int counter2=0;

    double lowerbound;
    //slices = 13;
    std::cout << "Number of slices: " << slices << std::endl;
    #pragma omp parallel
    {
        tour.bound = Parallel_first_lbound(distances, min);

        double private_lb = 0;
        double min1;
        double min2;

        #pragma omp for  private(min1, min2) schedule(dynamic) nowait
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
        
        #pragma omp single
        {
            tour.bound = lowerbound/2;
            tour.tour.push_back(0); 
            tour.cost = 0; 
            best_tour.tour.push_back(0);
            best_tour.cost = max_value;

            if (tour.bound >= best_tour.cost){
                omp_set_num_threads(0); 
            }
        }

        #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
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
            {
                tour_vector_2.push_back(new_tour);
            }
        }


        #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
        for (int i= 0; i < tour_vector_2.size(); i++){
            for (int v = 0; v < neighbors[tour_vector_2[i].tour.back()].size(); v++){
                
                neighbor = neighbors[tour_vector_2[i].tour.back()][v];

                if (std::find(tour_vector_2[i].tour.begin(), tour_vector_2[i].tour.end(), neighbor) != tour_vector_2[i].tour.end()) {
                    continue;
                }

                distance = distances[tour_vector_2[i].tour.back()][neighbor];
                new_tour.bound = Serial_compute_lbound(distance, min, tour_vector_2[i].tour.back(), neighbor, tour_vector_2[i].bound);
                
                if (new_tour.bound > max_value){ 
                    continue;
                }

                new_tour.tour = tour_vector_2[i].tour;
                new_tour.tour.push_back(neighbor); 
                new_tour.cost = tour_vector_2[i].cost + distance;
                

                #pragma omp critical
                {
                    tour_vector_3.push_back(new_tour);
                }
            }
        }

        #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
        for (int i= 0; i < tour_vector_3.size(); i++){
            for (int v = 0; v < neighbors[tour_vector_3[i].tour.back()].size(); v++){
                
                neighbor = neighbors[tour_vector_3[i].tour.back()][v];

                if (std::find(tour_vector_3[i].tour.begin(), tour_vector_3[i].tour.end(), neighbor) != tour_vector_3[i].tour.end()) {
                    continue;
                }

                distance = distances[tour_vector_3[i].tour.back()][neighbor];
                new_tour.bound = Serial_compute_lbound(distance, min, tour_vector_3[i].tour.back(), neighbor, tour_vector_3[i].bound);
                
                if (new_tour.bound > max_value){ 
                    continue;
                }

                new_tour.tour = tour_vector_3[i].tour;
                new_tour.tour.push_back(neighbor); 
                new_tour.cost = tour_vector_3[i].cost + distance;
                

                #pragma omp critical
                {
                    tour_vector_4.push_back(new_tour);
                }
            }
        }

        #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
        for (int i= 0; i < tour_vector_4.size(); i++){
            for (int v = 0; v < neighbors[tour_vector_4[i].tour.back()].size(); v++){
                
                neighbor = neighbors[tour_vector_4[i].tour.back()][v];

                if (std::find(tour_vector_4[i].tour.begin(), tour_vector_4[i].tour.end(), neighbor) != tour_vector_4[i].tour.end()) {
                    continue;
                }

                distance = distances[tour_vector_4[i].tour.back()][neighbor];
                new_tour.bound = Serial_compute_lbound(distance, min, tour_vector_4[i].tour.back(), neighbor, tour_vector_4[i].bound);
                
                if (new_tour.bound > max_value){ 
                    continue;
                }

                new_tour.tour = tour_vector_4[i].tour;
                new_tour.tour.push_back(neighbor); 
                new_tour.cost = tour_vector_4[i].cost + distance;
                

                #pragma omp critical
                {
                    tour_vector_5.push_back(new_tour);
                }
            }
        }

        #pragma omp single
        std::cout << "Number of queues: " << tour_vector_2.size() << std::endl;

        #pragma omp for private(new_tour, distance, neighbor) schedule(dynamic)
        for (int i = 0; i < tour_vector_5.size(); i++){
            std::vector<PriorityQueue<Tour>> queue_slices(slices);

            for (int v = 0; v < neighbors[tour_vector_5[i].tour.back()].size(); v++){
                
                neighbor = neighbors[tour_vector_5[i].tour.back()][v];

                if (std::find(tour_vector_5[i].tour.begin(), tour_vector_5[i].tour.end(), neighbor) != tour_vector_5[i].tour.end()) {
                    continue;
                }

                distance = distances[tour_vector_5[i].tour.back()][neighbor];
                new_tour.bound = Serial_compute_lbound(distance, min, tour_vector_5[i].tour.back(), neighbor, tour_vector_5[i].bound);
                
                if (new_tour.bound > max_value){ 
                    continue;
                }

                new_tour.tour = tour_vector_5[i].tour;

                new_tour.tour.push_back(neighbor); 
                new_tour.cost = tour_vector_5[i].cost + distance;
                


                for (int partition = 0; partition < slices; partition++){
                    int part = (partition+1) * neighbors[tour_vector_5[i].tour.back()].size()/slices;
                    if (v <= part){
                        queue_slices[partition].push(new_tour);
                        break;
                    }
                }
            }

            for (int partition = 0; partition < slices; partition++){
                if (!queue_slices[partition].empty()){
                    #pragma omp critical
                    queues.push_back(queue_slices[partition]);
                }
            }
        }
        
        #pragma omp single
        std::cout << "New number of queues: " << queues.size() << std::endl;

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