#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
//#include "/home/qqueke/PDC/PDC---Projeto/queue"
#include <string.h>

double Initial_LB(double **Distances, int N){
    double lowerbound = 0;
    double min1 = INT_MAX;
    double min2 = INT_MAX;

    //Paralelizável se utilizarmos vetores min1 e min2 com indíces
    //para guardar resultados de cada thread e finalmente somar tudo.
    
    for (int row = 0; row < N; row++){
        for (int column = 0; column < N; column++){
            if (Distances[row][column] < min1) {
                min1 = Distances[row][column];
            }
            else if (Distances[row][column] >= min1 && Distances[row][column] < min2){
                min2 = Distances[row][column];
            }
        }
        lowerbound = min1 + min2;
    }
    return lowerbound/2;
}

double Compute_LB(double **Distances, int N, int f, int t, double LB){
    double cf;
    double ct;
    double min1 = INT_MAX;
    double min2 = INT_MAX;

    //Melhor abordagem parece ser separar os dois for

    int row = f;
    for (int column = 0; column < N; column++){
        if (Distances[row][column] < min1) {
            min1 = Distances[row][column];
        }
        else if (Distances[row][column] >= min1 && Distances[row][column] < min2){
            min2 = Distances[row][column];
        }
    }

    if (Distances[f][t] >= min2) {
        cf = min2;
    }
    else{
        cf = min1;
    }


    row = t;
    for (int column = 0; column < N; column++){
        if (Distances[row][column] < min1) {
            min1 = Distances[row][column];
        }
        else if (Distances[row][column] >= min1 && Distances[row][column] < min2){
            min2 = Distances[row][column];
        }
    }

    if (Distances[f][t] >= min2) {
        ct = min2;
    }
    else{
        ct = min1;
    }

    return LB + Distances[f][t] - (cf - ct)/2;

}

void TSPBB(double **Distances, int N, double BestTourCost){
    
    //Estas variaveis ainda nem tem tipos certos,  simplesmente nao quis que estivessem
    //sublinhadas c os erros
    int row;
    int Tour[N+1];
    int newTour[N+1];
    int BestTour[N+1];
    for (int row = 0; row < N+1; row++)
    {
        if (row == 0){
            Tour[row] = 0;
            newTour[row] = 0;
            BestTour[row] = 0;
        }
        else{
            Tour[row] = -1;
            newTour[row] = -1;
            BestTour[row] = -1;
        }
    }

    double LB = Initial_LB(Distances, N);
    double newBound;
    double Cost, newCost, Bound, Length;
    int Node;

    priority_queue_t * Queue = queue_create(Tour, 0, LB, 1, 0);

    while (queue_pop(Queue) != NULL){
        (Tour, Cost, Bound, Length, Node) = queue_pop(Queue);
        
        if (Bound >= BestTourCost){
            return BestTour, BestTourCost;
        }

        if (Length == N){
            if (Cost + Distances[Node][0] < BestTourCost){
                for (row = 0; row < N; row++){
                    if (Tour[row] == -1){
                        break;
                    }
                    BestTour[row] = Tour[row];
                }
                BestTour[row] = 0;
                BestTourCost = Cost + Distances[Node][0];
            }
        }
        else{
            //for each neighbor v of Node and v doesnt belong to Tour
            for (int v = 0; v < N; v++){
                //If v is neighbor and doesnt belong to tour then:
                //Aqui damos update considerando v e Node né?
                newBound = Compute_LB(Distances, N, Node, v, LB);
                
                if (newBound > BestTourCost){
                    continue;
                }
                
                for (row = 0; row < N; row++){
                    if (Tour[row] == -1){
                        break;
                    }
                    newTour[row] = Tour[row];
                }

                newTour[row] = v;
                newCost = Cost + Distances[Node][v];
                queue_push(Queue, (newTour, newCost, newBound, Length + 1, v));
            }
        }
    }
    return BestTour, BestTourCost;
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

    double **Distances = (double **)malloc(num_cities * sizeof(double *));
    for (int row = 0; row < num_cities; row++) {
        Distances[row] = (double *)malloc(num_cities * sizeof(double));
    }

    for (int row = 0; row < num_cities; row++) {
        for (int column= 0; column< num_cities; column++) {
            Distances[row][column] = 0;
        }
    }
    
    for (int row = 0; row < num_roads; row++) {
        int city1, city2;
        double distance;
        fscanf(fp, "%d %d %lf", &city1, &city2, &distance);
        Distances[city1][city2] = distance;
        Distances[city2][city1] = distance;
    }

    fclose(fp);

    int BestTour[num_cities];
    for (int row = 0; row < num_cities; row ++){
        BestTour[row] = -1;
    }
    int BestTourCost;

    exec_time = -omp_get_wtime();

    (BestTour, BestTourCost) = TSPBB(Distances, num_cities, max_value);
    
    exec_time += omp_get_wtime();

    
    fprintf(stderr, "%.1fs\n", exec_time);

    //This means the graph is disconnected
    for (int row = num_cities - 1; row >= 0; row --){
        if (BestTour[row] == -1){
            printf("NO SOLUTION");
            break;
        }
    }
    //No solution that has a better value than the max admited
    if (BestTourCost > BestTourCost){
        printf("NO SOLUTION");
    }
    else{
        printf("%.1lf\n", BestTourCost);
        for (int row = 0; row < num_cities; row ++){
            printf("%d", BestTour[row]);
        }
    }

    return 0;
}
