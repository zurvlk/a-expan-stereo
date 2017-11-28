#include "graph.h"
#include "bmp.h"

double energy(Graph *G, int *label,  double T, int lambda, int width, int *I_left, int *I_right);
double pairwise(double i, double j, double T, int lambda);
double data(int i, int label, int width, int *I_left, int *I_right);
int update_labels(Graph *G, int *label, int alpha, int width, double T, int lambda, int *t, int *I_left, int *I_right);
void set_capacity(Graph *G, int *label, int width, int alpha, double T, int lambda, int *I_left, int *I_right);
void set_all_edge(Graph *G, int height, int width);
double err_rate(img output, img truth, int scale);

