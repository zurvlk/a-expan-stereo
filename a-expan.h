#include "graph.h"
#include "bmp.h"

double energy(Graph *G ,int *label, int *I, double T, int lambda);
double pairwise(double i, double j, double T, int lambda);
double data(int *I, int a, double i);
int update_labels(Graph *G, int alpha, int *label, int *t, int *I, double T, int lambda);
void set_capacity(Graph *G, int alpha, int *label, int *I, double T, int lambda);
void set_all_edge(Graph *G, int height, int width);
void capacity(Graph *G,int x[], int y[],int a);
