#include "graph.h"
#include "bmp.h"
typedef struct __INPUT_BITMAPS__ {
    img raw_left;
    img raw_right;
    img truth;
    img output;
    int *left;
    int *right;
    int width;
    int height;
}Image;

void label2Bmp(int *label, Image image, int scale, char filename[]);
double readStrBmp(Image *image, char filename[], int scale);
double energy(Graph *G, int *label,  double T, int lambda, Image image);
double pairwise(double i, double j, double T, int lambda);
double data(int i, int label, Image image);
// int update_labels(Graph *G, int *label, int alpha, double T, int lambda, int *t, int *I);
void set_capacity(Graph *G, int *label, int alpha, double T, int lambda, Image image);
void set_all_edge(Graph *G, int height, int width);
double err_rate(img output, img truth, int scale);

void label2Bmp(int *label, Image image, int scale, char output_file[]);