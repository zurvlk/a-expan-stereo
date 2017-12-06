#include "bmp.h"
#include "a_estr.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <string.h>
#define INF DBL_MAX

double readStrBmp(Image *image, char filename[], int scale) {
    int i, j, grids_node;
    char imgleft[100];
    char imgright[100];
    char imgtruth[100];

    strcpy(imgleft, filename);
    strcpy(imgright, filename);
    strcpy(imgtruth, filename);

    strcat(imgleft, "left.bmp");
    strcat(imgright, "right.bmp");
    strcat(imgtruth, "truth.bmp");

    ReadBmp(imgleft, &(image->raw_left));
    ReadBmp(imgright, &(image->raw_right));
    ReadBmp(imgtruth, &(image->truth));
    ReadBmp(imgtruth, &(image->output));

    image->width = image->raw_left.width;
    image->height = image->raw_left.height;
    grids_node = image->height * image->width;

    Gray(&(image->raw_left), &(image->raw_left));
    Gray(&(image->raw_right), &(image->raw_right));
    Gray(&(image->truth), &(image->truth));


    if ((image->left = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->image.left]\n");
        exit(EXIT_FAILURE);
    }
    if ((image->right = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->image.right]\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i <  image->height; i++) {
        for (j = 0; j < image->width; j++) {
            image->left[i * image->width + j + 1] = image->raw_left.data[i][j].r / scale;
        }
    }
    for (i = 0; i <  image->height; i++) {
        for (j = 0; j < image->width; j++) {
            image->right[i * image->width + j + 1] = image->raw_right.data[i][j].r / scale;
        }
    }
    return grids_node;
}

double dabs(double a, double b) {
    return a - b > 0 ? a - b : b - a;
}
double fmin(double i, double j) {
    return i < j  ? i : j;
}

double fmin3(double i, double j, double k) {
    return fmin(fmin(i, j), k);
}

double between(double a, double b, double c){
    if((a <= b && b <= c) || (c <= b && b <= a))    return 0;
    else return 1;
}

double pairwise(double i, double j, double T, int lambda) {
    return lambda * fmin((i - j) * (i - j), T);
}

double data(int i, int label, int width, int *I_left, int *I_right) {
    double data = 0;
    // return 1.0 * dabs(label, I_left[i]);
    //leftの中のものがどこにあるか
    if ((i - 1) / width == (i - label - 1) / width) {
        data = (I_left[i] - I_right[i - label]) * (I_left[i] - I_right[i - label]);
    }else data = INF;

    return sqrt(data);
}

double Dt(int x, Image *image, int i, int j) {
    double d_m, d_l, d_r, d_m2, d_l2, d_r2;
    double  I, I_1 = 0, I_2 = 0;

	if(x < 0 || x > image->label_max ) return INF;
    if (j - x < 0) return 10000;

    int flag = 1;
    while(flag){
    
	    d_m = image->raw_left.data[i][j].r;
	    d_m2 =  image->raw_right.data[i][j - x].r;

	    if((j - x - 1 >= 0) && (j - x + 1 <= image->raw_right.width - 1)){

	        d_l2 = (image->raw_right.data[i][j - x].r + image->raw_right.data[i][j - x - 1].r) / 2.0;
            d_r2 = (image->raw_right.data[i][j - x].r + image->raw_right.data[i][j - x + 1].r) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        if(between(d_r2, d_m, d_m2) == 0)    break;
	        I_1 = fmin3(fabs(d_m - d_l2), fabs(d_m - d_m2), fabs(d_m - d_r2));              /*特に制約なし*/

	    }else if(j - x - 1 < 0){

	        d_r2 = (image->raw_right.data[i][j - x].r + image->raw_right.data[i][j - x + 1].r) / 2.0;
	        if(between(d_m2, d_m, d_r2) == 0)    break;
	        I_1 = fmin(fabs(d_m - d_m2), fabs(d_m - d_r2));                                /*Rの左端が出る*/

	    }else{

	        d_l2 = (image->raw_right.data[i][j - x].r + image->raw_right.data[i][j - x - 1].r) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        I_1 = fmin(fabs(d_m - d_l2), fabs(d_m - d_m2));                                /*Rの右端が出る*/
	    }

	    if((j != 0) && (j != image->raw_left.width - 1)){

	        d_l = (image->raw_left.data[i][j].r + image->raw_left.data[i][j - 1].r) / 2.0;
	        d_r = (image->raw_left.data[i][j].r + image->raw_left.data[i][j + 1].r) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        if(between(d_r, d_m2, d_m) == 0)    break;
	        I_2 = fmin3(fabs(d_l - d_m2), fabs(d_m - d_m2), fabs(d_r - d_m2));              /*特に制約なし*/

	    }else if(j == 0){ 

	        d_r = (image->raw_left.data[i][j].r + image->raw_left.data[i][j + 1].r) / 2.0;
	        if(between(d_m, d_m2, d_r) == 0)    break;
	        I_2 = fmin(fabs(d_m - d_m2), fabs(d_r - d_m2));                                /*Lの左端が出る*/

	    }else{

	        d_l = (image->raw_left.data[i][j].r + image->raw_left.data[i][j - 1].r) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        I_2 = fmin(fabs(d_l - d_m2), fabs(d_m - d_m2));                                /*Lの右端が出る*/
	    }
	    flag = 0;
	}
    I = fmin(I_1, I_2);
    I = fmin(I, 20);
    return I;
    //I *= I; 
    //return I * I;                                   /*Dを2乗する*/
    //return C_D * (I);

}



double energy(Graph *G, int *label,  double T, int lambda, Image image) {

    double energy = 0;
    //* Dterm
    // for (i = 1; i <= G->n - 2; i++) {
    //     energy += data(i, label[i], image.width, image.left, image.right);
    //     // energy += data(I_left[i], label[i]);
    // }
    // */

    for(int i = 1; i < G->n - 1; i++){
        energy += Dt(label[i], &image, i / image.width, i % image.width);
    }
    // Vterm
    for (int i = 1; i <= G->m - 2 * (G->n - 2); i++) {
        energy += pairwise(label[G->tail[i]], label[G->head[i]], T, lambda);
    }
    return energy;
}

void set_capacity(Graph *G, int *label, int alpha, double T, int lambda, Image image) {
    int i, s2i_begin, i2t_begin, grids_node;
    double A, B, C, D, temp;
    //　source->各ノード　を示す枝の開始位置
    s2i_begin = G->m - 2 * (G->n - 2) + 1;
    //　各ノード->sink　を示す枝の開始位置
    i2t_begin = G->m - (G->n - 2) + 1;
    // G中のノードのうちsource, sinkを除く画像部分のノードの総数
    grids_node = G->n - 2;

    for (i = 0; i <= G->m + 1; i++) G->capa[i] = 0;
    for (i = s2i_begin; i < i2t_begin; i++) {
        //labelがalphaでない->capacity設定s2i
        //labelがalpha->枝を削除
        if (label[G->head[i]] == alpha) {
            deleteAdjEdge(G, i);
            deleteAdjEdge(G, i + grids_node);
        }
    }

    // set Dterm
    for(i = 1; i < G->n - 1; i++){
        temp = Dt(label[i], &image, i / image.width, i % image.width) - Dt(alpha, &image, i / image.width, i % image.width);
        if(temp > 0){
            G->capa[s2i_begin - 1 + i] +=  temp;
        }else{
            G->capa[i2t_begin - 1  + i] -= temp;
        }
    }
    // set Vterm
    for(i = 1; i < s2i_begin; i++) {
        //head,tailが共にalpha->枝を削除
        if (label[G->tail[i]] != alpha && label[G->head[i]] != alpha) {
            A = pairwise(label[G->tail[i]], label[G->head[i]], T, lambda);
            B = pairwise(label[G->tail[i]], alpha, T, lambda);
            C = pairwise(alpha, label[G->head[i]], T, lambda);
            D = pairwise(alpha, alpha, T, lambda);

            A = fmin(A, B + C - D);
            G->capa[i] = B + C - A - D;
            // source->i
            G->capa[s2i_begin - 1 + G->tail[i]] += B - D;
            // i->sink
            temp = B - A;
            if(temp > 0) {
                G->capa[i2t_begin - 1 + G->head[i]] += temp;
            } else {
                G->capa[i2t_begin - 1 + G->head[i]] += temp;
            }
        } else {
            deleteAdjEdge(G, i);
        }
    }
    // 各 s->i->tに対し、capacityが小さいEdgeのcapacityを、両方のcapacityから引く(高速化処理)
    for (i = s2i_begin; i < i2t_begin; i++) {
        temp = fmin(G->capa[i], G->capa[i + grids_node]);
        G->capa[i] -= temp;
        G->capa[i + grids_node] -= temp;
    }
    return;
}

void set_all_edge(Graph *G, int height, int width) {
    int i, j, tail, head, grids_node, source, sink, edge_count;

    grids_node = height * width;
    edge_count = 1;

    // horizontal
    for (i = 1; i < height + 1; i++){
        for (j = 1; j < width; j++){
            tail = (i - 1) * width + j;
            head = (i - 1) * width + j + 1;
            setEdge(G, edge_count, tail, head, 0);
            edge_count++;
        }
    }

    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++){
            tail = (i - 1) * width + j;
            head = i * width + j;
            setEdge(G, edge_count, tail, head, 0);
            edge_count++;
        }
    }

    source = grids_node + 1;
    sink = source + 1;

    setSource(G, source);
    setSink(G, sink);

    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, source, i, 0);
        setEdge(G, edge_count + G->n - 2, i, sink, 0);
        edge_count++;
    }
    return;
}

double err_rate(img output, img truth, int scale) {
    int i, error_count = 0;
    double err;
    if (truth.data[0][0].r) {
        for(i = 1; i <= (output.height) * (output.width); i++) {   
            if (abs(output.data[(i - 1) / output.width][(i - 1) % output.width].r - truth.data[(i - 1) / truth.width][(i - 1) % truth.width].r ) 
                >= scale + 1) {
                error_count++;
            }                            
        }
    } else {
        for(i = 1; i <= (output.height) * (output.width); i++) {
            if ((i - 1) / output.width >= scale && (i - 1) % output.width >= scale &&
                (i - 1) / output.width <= output.height - scale && (i - 1) % output.width <= output.width - scale) {
                if (abs(output.data[(i - 1) / output.width][(i - 1) % output.width].r - truth.data[(i - 1) / truth.width][(i - 1) % truth.width].r ) 
                    >= scale + 1) {
                    error_count++;
                }
            }                            
        }
    }

    err = 100 * error_count / (double)(truth.height * truth.width);
    return err;
}