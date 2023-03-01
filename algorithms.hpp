#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <omp.h>
#include <climits>
#include <algorithm>
#include <cfloat>
#include "queue.hpp"

typedef struct Tour{
    std::vector<int> tour;
    double cost;
    double bound;
    int current_node;

    friend bool operator>(const Tour& self, const Tour& other){
        if (self.bound != other.bound) {
            return self.bound > other.bound;
        }
        else {
            return self.current_node > other.current_node;
        }
    }
}Tour;

//This is a print function to help with the queue.print
void print_tour(const Tour& tour) {
    std::cout << "Tour: {";
    for (int i : tour.tour) {
        std::cout << i << " ";
    }
    std::cout << "}, Current Node: " << tour.current_node << ", Tour cost: " << tour.cost << ", Tour Lowerbound: " << tour.bound << std::endl;
}

double compute_lbound(const std::vector<std::vector<double>>& Distances, const std::vector<double>& min1, const std::vector<double>& min2, int f, int t, double LB) {
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

double Serial_first_lbound(const std::vector<std::vector<double>> &Distances, std::vector<double> &min1, std::vector<double> &min2)
{
    double lowerbound = 0;

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
    return lowerbound / 2;
}

double Parallel_first_lbound(const std::vector<std::vector<double>> &Distances, std::vector<double> &min1, std::vector<double> &min2)
{
    double lowerbound = 0;

    //Alternitively we could have a lowerbound array for each thread that would then be added up to the total lowerbound
    //Cannot remember what the +: meant
    #pragma omp parallel for reduction(+: lowerbound) schedule(dynamic,64)
    for (int row = 0; row < Distances.size(); row++) {
        #pragma omp parallel for schedule(dynamic, 64)
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
    return lowerbound / 2;
}

Tour Serial_tsp_bb(const std::vector<std::vector<double>>& Distances, int N, double BestTourCost){

    //Lowerbound vectors for min1 and min2
    std::vector<double> min1(Distances.size(), INT_MAX);
    std::vector<double> min2(Distances.size(), INT_MAX);
    
    //Priority queue
    PriorityQueue<Tour> queue;
    //Tours include the path, cost, bound and the current node
    Tour initialTour, BestTour, newTour;

    initialTour.tour.push_back(0); // Tour ← {0}
    initialTour.cost = 0; // Tour Cost <- 0
    initialTour.bound = Serial_first_lbound(Distances, min1, min2); //Tour lowerbound <- Initial lowerbound guess
    initialTour.current_node = 0; //Tour current node  <- 0
    queue.push(initialTour); //Queue ← (Tour, 0, LB, 1, 0) (Tour, Cost, Bound, Length, Current city)

    BestTour.tour.push_back(0);
    BestTour.cost = BestTourCost;

    while (!queue.empty()){ //while Queue ̸= {} do
        initialTour = queue.pop(); // (Tour, Cost, Bound, Length, N ode) ← Queue.pop()

        if (initialTour.bound >= BestTour.cost){// if Bound ≥ BestT ourCost then
            return BestTour; // return BestT our, BestT ourCost
        }

        if (initialTour.tour.size() == N){ // if Length = N then
            if (initialTour.cost + Distances[initialTour.current_node][0] < BestTour.cost){ // if Cost + Distances(Node, 0) < BestTourCost then
                BestTour.tour = initialTour.tour; // BestTour ← Tour ∪ {0}
                BestTour.tour.push_back(0); // BestTour ← Tour ∪ {0}
                BestTour.cost = initialTour.cost + Distances[initialTour.current_node][0];//BestT ourCost ← Cost + Distances(N ode, 0)
            } 
        }
        else{
            //for each neighbor v of Node and v ̸∈ Tour do
            for (int v = 0; v < N; v++){
                //If v is neighbor and doesnt belong to tour then:
                //If it is not a neighbor quit
                if (Distances[initialTour.current_node][v] == DBL_MAX){
                    continue;
                }

                //This find method will be replaced
                //If it already belongs to tour quit
                if (std::find(initialTour.tour.begin(), initialTour.tour.end(), v) != initialTour.tour.end()) {
                    continue;
                }

                //newBound ← updated lower bound on tour cost
                newTour.bound = compute_lbound(Distances, min1, min2, initialTour.current_node, v, initialTour.bound);
                //initialTour.bound = compute_lbound(Distances, min1, min2, initialTour.current_node, v, initialTour.bound);

                if (newTour.bound > BestTour.cost){ //if newBound > BestT ourCost then
                    continue;
                }

                newTour.tour = initialTour.tour;
                newTour.tour.push_back(v); // newTour ← Tour ∪ {v}
                newTour.cost = initialTour.cost + Distances[initialTour.current_node][v]; //newCost ← cost + Distances(N ode, v)
                newTour.current_node = v;
                queue.push(newTour); //Queue.add((newTour, newCost, newBound, Length + 1, v)), v is the new current node
            }//end for
        }//end if
    }//end while
    return BestTour;
}//end procedure
