#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <climits>

typedef struct Tour{
    std::vector<int> tour;
    double cost;
    double bound;
    int length;
    int current_node;

    friend bool operator>(const Tour& other) const {
        return bound > other.bound;
    }
}Tour;


double Initial_LB(const std::vector<std::vector<double>>& distances, std::vector<double>& min1, std::vector<double>& min2) {
    double lowerbound = 0;

    for (int row = 0; row < distances.size(); row++) {
        for (int column = 0; column < distances.size(); column++) {
            if (distances[row][column] < min1[row]) {
                min2[row] = min1[row];
                min1[row] = distances[row][column];
            }
            else if (distances[row][column] >= min1[row] && distances[row][column] < min2[row]) {
                min2[row] = distances[row][column];
            }
        }
        lowerbound += min1[row] + min2[row];
    }
    return lowerbound / 2;
}

double Compute_LB(const std::vector<std::vector<double>>& distances, const std::vector<double>& min1, const std::vector<double>& min2, int f, int t, double LB) {
    double cf = 0;
    double ct = 0;

    if (distances[f][t] >= min2[f]) {
        cf = min2[f];
    }
    else {
        cf = min1[f];
    }

    if (distances[t][f] >= min2[t]){
        ct = min2[t];
    }
    else {
        ct = min1[t];
    }

    double lower_bound = LB + cf + ct;

    return lower_bound;
}

Tour TSPBB(const std::vector<std::vector<double>>& Distances, int N, double BestTourCost){
    
    //Estas variaveis ainda nem tem tipos certos,  simplesmente nao quis que estivessem
    //sublinhadas c os erros
    int row;
    std::vector<double> min1(Distances.size(), INT_MAX);
    std::vector<double> min2(Distances.size(), INT_MAX);
    PriorityQueue<Tour> queue;
    Tour tour, newTour, BestTour;

    newTour.tour.push_back(N+1);
    BestTour.tour.push_back(N+1);
    BestTour.cost = BestTourCost;

    tour.tour.push_back(N+1);
    tour.bound = Initial_LB(Distances, min1, min2);
    tour.cost = 0;
    tour.length = 1;
    tour.current_node = 0;


    while (!queue.empty()){
        tour = queue.pop();
        
        if (tour.bound >= BestTour.cost){
            return BestTour;
        }

        if (tour.length == N){
            if (tour.cost + Distances[aux.current_node][0] < BestTour.cost){
                BestTour.tour = tour.tour;
                BestTour.tour[tour.tour.size()] = 0;
                BestTour.cost = tour.cost + Distances[tour.current_node][0];
            }
        }
        else{
            //for each neighbor v of Node and v doesnt belong to Tour
            for (int v = 0; v < N; v++){
                //If v is neighbor and doesnt belong to tour then:
                //Aqui damos update considerando v e Node né?
                newTour.bound = Compute_LB(Distances, min1, min2, tour.current_node, v, tour.bound);
                
                if (newTour.bound > BestTour.cost){
                    continue;
                }

                newTour.tour = tour.tour;
                newTour.tour[tour.tour.size()] = v;
                newTour.cost = tour.cost + Distances[tour.current_node][v];
                newTour.length = tour.length + 1;
                newTour.current_node = v;
                queue.push(newTour);
            }
        }
    }
    return BestTour;
}


int main(int argc, char *argv[]) {
    double exec_time;

    if (argc != 3) {
        printf("Usage: tsp <cities file> <max-value>\n");
        return 1;
    }

    char *cities_file = argv[1];
    int max_value = atoi(argv[2]);

    FILE *fp = fopen(cities_file, "r");
    if (fp == NULL) {
        printf("Error: Unable to open file %s\n", cities_file);
        return 1;
    }

    int num_cities, num_roads;
    fscanf(fp, "%d %d", &num_cities, &num_roads);

    std::vector<std::vector<double>> Distances (num_cities, std::vector<double>(num_cities));
    
    for (int row = 0; row < num_roads; row++) {
        int city1, city2;
        double distance;
        fscanf(fp, "%d %d %lf", &city1, &city2, &distance);
        Distances[city1][city2] = distance;
        Distances[city2][city1] = distance;
    }

    fclose(fp);

    exec_time = -omp_get_wtime();

    Tour BestTour = TSPBB(Distances, num_cities, max_value);
    
    exec_time += omp_get_wtime();

    
    fprintf(stderr, "%.1fs\n", exec_time);

    //No solution that has a better value than the max admited
    if (BestTour.cost > max_value){
        printf("NO SOLUTION");
    }
    //This means the graph is disconnected
    else if (BestTour.tour.size() != num_cities +1){
        printf("NO SOLUTION");
    }
    //Valid solution
    else{
        printf("%.1lf\n", BestTour.cost);
        for (int row = 0; row < num_cities; row ++){
            printf("%d", BestTour.tour[row]);
        }
    }

    return 0;
}
/*

double TSPBB(const std::vector<std::vector<double>>& distances) {
    int n = distances.size();
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    double best_dist = std::numeric_limits<double>::infinity();
    std::vector<int> best_perm;

    // Compute lower bound using the initial lower bound function.
    double LB = Initial_LB(distances);

    // Create a priority queue to store nodes in the search tree.
    std::priority_queue<Node, std::vector<Node>, NodeComparator> pq;

    // Initialize the search tree with the root node.
    pq.push(Node(0, perm, LB));

    while (!pq.empty()) {
        // Get the node with the highest priority.
        Node curr_node = pq.top();
        pq.pop();

        // If the lower bound of this node is worse than the current best distance, prune this branch.
        if (curr_node.LB >= best_dist) {
            continue;
        }

        // If the current node is a leaf node, update the best distance and permutation.
        if (curr_node.level == n) {
            double curr_dist = Compute_Distance(distances, curr_node.perm);
            if (curr_dist < best_dist) {
                best_dist = curr_dist;
                best_perm = curr_node.perm;
            }
            continue;
        }

        // Branch on all possible next cities.
        for (int i = curr_node.level; i < n; i++) {
            std::vector<int> child_perm = curr_node.perm;
            std::swap(child_perm[curr_node.level], child_perm[i]);

            // Compute the lower bound for the child node using the Compute_LB function.
            double child_LB = Compute_LB(distances, curr_node.min1, curr_node.min2, curr_node.level, i, curr_node.LB);

            // Add the child node to the priority queue.
            pq.push(Node(curr_node.level + 1, child_perm, child_LB, curr_node.min1, curr_node.min2));
        }
    }

    return best_dist;
}



*/