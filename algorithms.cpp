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
    std::cout << "}, Current Node: " << tour.tour.back() << ", Tour cost: " << tour.cost << ", Tour Lowerbound: " << tour.bound << std::endl;
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

    //Tours include the path, cost, bound and the current node
    Tour tour, best_tour, new_tour;


    tour.tour.push_back(0); // Tour ← {0}
    tour.cost = 0; // Tour Cost <- 0
    tour.bound = Serial_first_lbound(distances, min); //Tour lowerbound <- Initial lowerbound guess
    //tour.current_node = 0; //Tour current node  <- 0
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
                //new_tour.current_node = v;
                queue.push(new_tour); //Queue.add((new_tour, newCost, newBound, Length + 1, v)), v is the new current node
            }//end for
        }//end if
    }//end while
    return best_tour;
}//end procedure

