#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <limits.h>
#include <string.h>
#include <vector>
#include <climits>
#include <algorithm>
#include <cfloat>
#include "queue.hpp"


// #pragma omp for reduction(+: NodesNoOutEdges) schedule(dynamic,64)
void print_tour(const Tour& tour) {
    std::cout << "Tour: {";
    for (int i : tour.tour) {
        std::cout << i << " ";
    }
    std::cout << "}, Current Node: " << tour.current_node << ", Tour cost: " << tour.cost << ", Tour Lowerbound: " << tour.bound << std::endl;
}

double Initial_LB(const std::vector<std::vector<double>> &Distances, std::vector<double> &min1, std::vector<double> &min2)
{
    double lowerbound = 0;
    //Este threshold é a partir do qual vale a pena o overhead de criar threads 
    int threshold = 100;
    
    //Parallel
    if (Distances.size() > threshold){
        #pragma omp parallel for reduction(+: lowerbound) schedule(dynamic,64)
            for (int row = 0; row < Distances.size(); row++) {
                #pragma omp parallel for schedule(dynamic,64)
                for (int column = 0; column < Distances.size(); column++) {
                    if (Distances[row][column] < min1[row]) {
                        min2[row] = min1[row];
                        min1[row] = Distances[row][column];
                    }
                    else if (Distances[row][column] >= min1[row] && Distances[row][column] < min2[row]) {
                        min2[row] = Distances[row][column];
                    }
                }
            lowerbound += min1[row] + min2[row];
        }
    }
    //Serial
    else{
        for (int row = 0; row < Distances.size(); row++) {
            for (int column = 0; column < Distances.size(); column++) {
                if (Distances[row][column] < min1[row]) {
                    min2[row] = min1[row];
                    min1[row] = Distances[row][column];
                }
                else if (Distances[row][column] >= min1[row] && Distances[row][column] < min2[row]) {
                    min2[row] = Distances[row][column];
                }
            }
            lowerbound += min1[row] + min2[row];
        }
    }

    return lowerbound / 2;
}

double Compute_LB(const std::vector<std::vector<double>>& Distances, const std::vector<double>& min1, const std::vector<double>& min2, int f, int t, double LB) {
    double cf = 0;
    double ct = 0;

    if (Distances[f][t] >= min2[f]) {
        cf = min2[f];
    }
    else {
        cf = min1[f];
    }

    if (Distances[t][f] >= min2[t]){
        ct = min2[t];
    }
    else {
        ct = min1[t];
    }

    double lower_bound = LB + Distances[f][t] - (cf + ct)/2;

    return lower_bound;
}

Tour TSPBB(const std::vector<std::vector<double>>& Distances, int N, double BestTourCost){

    //Lowerbound vectors for min1 and min2
    std::vector<double> min1(N, INT_MAX);
    std::vector<double> min2(N, INT_MAX);
    
    //Priority queue
    PriorityQueue<Tour> queue;
    //Tours include the path, cost, bound and the current node
    Tour initialTour, BestTour, newTour;

    initialTour.tour.push_back(0); // Tour ← {0}
    initialTour.cost = 0; // Tour Cost <- 0
    initialTour.bound = Initial_LB(Distances, min1, min2); //Tour lowerbound <- Initial lowerbound guess
    initialTour.current_node = 0; //Tour current node  <- 0
    //printf("Current node: %d\n", initialTour.current_node);   
    //printf("Lowerbound: %.1lf\n", initialTour.bound);
    queue.push(initialTour); //Queue ← (Tour, 0, LB, 1, 0) (Tour, Cost, Bound, Length, Current city)
    //std::cout << "Push:\n";
    //queue.print(&print_tour);
    BestTour.tour.push_back(0);
    BestTour.cost = BestTourCost;

    while (!queue.empty()){ //while Queue ̸= {} do
        initialTour = queue.pop(); // (Tour, Cost, Bound, Length, N ode) ← Queue.pop()
        //std::cout << "Pop:\n";
        //queue.print(&print_tour);

        if (initialTour.bound >= BestTour.cost){// if Bound ≥ BestT ourCost then
            return BestTour; // return BestT our, BestT ourCost
        }// end if

        if (initialTour.tour.size() == N){ // if Length = N then
            if (initialTour.cost + Distances[initialTour.current_node][0] < BestTour.cost){ // if Cost + Distances(Node, 0) < BestTourCost then
                BestTour.tour = initialTour.tour; // BestTour ← Tour ∪ {0}
                BestTour.tour.push_back(0); // BestTour ← Tour ∪ {0}
                BestTour.cost = initialTour.cost + Distances[initialTour.current_node][0];//BestT ourCost ← Cost + Distances(N ode, 0)
            } // end if
        }
        else{
            //Testar aqui v= 0
            //for each neighbor v of Node and v ̸∈ Tour do
            for (int v = 0; v < N; v++){
                //If v is neighbor and doesnt belong to tour then:
                //If it is not a neighbor quit
                if (Distances[initialTour.current_node][v] == DBL_MAX){
                    continue;
                }

                //Trocar por loop paralelizável
                //If it already belongs to tour quit
                if (std::find(initialTour.tour.begin(), initialTour.tour.end(), v) != initialTour.tour.end()) {
                    continue;
                }

                //newBound ← updated lower bound on tour cost
                newTour.bound = Compute_LB(Distances, min1, min2, initialTour.current_node, v, initialTour.bound);
                //initialTour.bound = Compute_LB(Distances, min1, min2, initialTour.current_node, v, initialTour.bound);

                if (newTour.bound > BestTour.cost){ //if newBound > BestT ourCost then
                    continue;
                }// end if

                //std::cout << "Current node: " << initialTour.current_node << std::endl;
                //std::cout << "Visited node: " << v << std::endl;
                //std::cout << "New lowerbound: " << initialTour.bound << std::endl;

                newTour.tour = initialTour.tour;
                newTour.tour.push_back(v);                                                // newTour ← Tour ∪ {v}
                newTour.cost = initialTour.cost + Distances[initialTour.current_node][v]; //newCost ← cost + Distances(N ode, v)
                newTour.current_node = v;
                queue.push(newTour); //Queue.add((newTour, newCost, newBound, Length + 1, v)), v is the new current node
                //std::cout << "Push:\n";
                //queue.print(&print_tour);
            }//end for
        }//end if
    }//end while
    return BestTour;
}//end procedure


int main(int argc, char *argv[]) {
    double exec_time;

    if (argc != 3) {
        std::cout << "Usage: tsp <cities file> <max-value>\n";
        return 1;
    }

    char *cities_file = argv[1];
    int max_value = atoi(argv[2]);

    FILE *fp = fopen(cities_file, "r");
    if (fp == NULL) {
        std::cout << "Error: Unable to open file " << cities_file << std::endl;
        return 1;
    }

    int num_cities, num_roads;
    if (!fscanf(fp, "%d %d", &num_cities, &num_roads)){
        std::cout << "Error";
        return 1;
    }

    std::vector<std::vector<double>> Distances (num_cities, std::vector<double>(num_cities));
    


    for (int row = 0; row < num_cities; row++) {
        for (int column = 0; column < num_cities; column++) {
            Distances[row][column] = DBL_MAX;
        }
    }

    for (int row = 0; row < num_roads; row++) {
        int city1, city2;
        double distance;
        if (!fscanf(fp, "%d %d %lf", &city1, &city2, &distance)){
            std::cout << "Error";
            return 1;
        }
        Distances[city1][city2] = distance;
        Distances[city2][city1] = distance;
    }

    fclose(fp);

    exec_time = -omp_get_wtime();

    Tour BestTour = TSPBB(Distances, num_cities, max_value);
    
    exec_time += omp_get_wtime();

    fprintf(stderr, "%lfs\n", exec_time);

    //No solution that has a better value than the max admited
    if (BestTour.cost > max_value){
        std::cout << "NO SOLUTION";
    }
    //This means the graph is disconnected
    else if (BestTour.tour.size() != num_cities +1){
        printf("NO SOLUTION");
    }
    //Valid solution
    else{
        std::cout << BestTour.cost << std::endl;
        std::cout << "Tour: {";

        for (int i : BestTour.tour) {
            std::cout << i << " ";
        }
        std::cout << "}" << std::endl;
    }

    return 0;
}
