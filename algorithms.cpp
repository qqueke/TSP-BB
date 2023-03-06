#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <omp.h>
#include <climits>
#include <algorithm>
#include <cfloat>
#include "queue.hpp"
#include "algorithms.hpp"


//This is a print function to help with the queue.print
void print_tour(const Tour& tour) {
    std::cout << "Tour: {";
    for (int i : tour.tour) {
        std::cout << i << " ";
    }
    std::cout << "}, Current Node: " << tour.current_node << ", Tour cost: " << tour.cost << ", Tour Lowerbound: " << tour.bound << std::endl;
}

double compute_lbound(const std::vector<std::vector<double>>& distances, std::vector<std::vector<double>> &min, int f, int t, double LB) {
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
        for (int column = 0; column < distances.size(); column++) {
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

    //Alternitively we could have a lowerbound array for each thread that would then be added up to the total lowerbound
    //Cannot remember what the +: meant
    #pragma omp parallel for reduction(+: lowerbound) schedule(dynamic,64)
    for (int row = 0; row < distances.size(); row++) {
        #pragma omp parallel for schedule(dynamic, 64)
        for (int column = 0; column < distances.size(); column++) {
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

Tour Serial_tsp_bb(const std::vector<std::vector<double>>& distances, int N, double BestTourCost){

    //Lowerbound vectors for min1 and min2
    std::vector<std::vector<double>> min (N, std::vector<double>(2));

    for (int row = 0; row < N; row++) {
        for (int column = 0; column < 2; column++) {
            min[row][column] = INT_MAX;
        }
    }


    //Priority queue
    PriorityQueue<Tour> queue;
    //Tours include the path, cost, bound and the current node
    Tour initialTour, BestTour, newTour;

    initialTour.tour.push_back(0); // Tour ← {0}
    initialTour.cost = 0; // Tour Cost <- 0
    initialTour.bound = Serial_first_lbound(distances, min); //Tour lowerbound <- Initial lowerbound guess
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
            if (initialTour.cost + distances[initialTour.current_node][0] < BestTour.cost){ // if Cost + distances(Node, 0) < BestTourCost then
                BestTour.tour = initialTour.tour; // BestTour ← Tour ∪ {0}
                BestTour.tour.push_back(0); // BestTour ← Tour ∪ {0}
                BestTour.cost = initialTour.cost + distances[initialTour.current_node][0];//BestT ourCost ← Cost + distances(N ode, 0)
            } 
        }
        else{
            //for each neighbor v of Node and v ̸∈ Tour do
            for (int v = 0; v < N; v++){
                //If v is neighbor and doesnt belong to tour then:
                //If it is not a neighbor quit
                if (distances[initialTour.current_node][v] == INT_MAX){
                    continue;
                }

                //This find method will be replaced
                //If it already belongs to tour quit
                if (std::find(initialTour.tour.begin(), initialTour.tour.end(), v) != initialTour.tour.end()) {
                    continue;
                }

                //newBound ← updated lower bound on tour cost
                newTour.bound = compute_lbound(distances, min, initialTour.current_node, v, initialTour.bound);
                //initialTour.bound = compute_lbound(distances, min1, min2, initialTour.current_node, v, initialTour.bound);

                if (newTour.bound > BestTour.cost){ //if newBound > BestT ourCost then
                    continue;
                }

                newTour.tour = initialTour.tour;
                newTour.tour.push_back(v); // newTour ← Tour ∪ {v}
                newTour.cost = initialTour.cost + distances[initialTour.current_node][v]; //newCost ← cost + distances(N ode, v)
                newTour.current_node = v;
                queue.push(newTour); //Queue.add((newTour, newCost, newBound, Length + 1, v)), v is the new current node
            }//end for
        }//end if
    }//end while
    return BestTour;
}//end procedure

