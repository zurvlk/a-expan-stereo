#include "bmp.h"
#include "a-expan.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#define INF DBL_MAX

double dabs(double a, double b) {
    return a - b > 0 ? a - b : b - a;
}
double fmin(double i, double j) {
    return i < j  ? i : j;
}

double pairwise(double i, double j, double T, int lambda) {
    return fmin(1.0 * dabs(i, j) , T);
}

double data(int *I, int a, double i) {
    return 1.0 * dabs(I[a], i) ;
}
double Dt(int x, img *tmp_l2, img *tmp_r2, int i, int j, int max_label, int image)
{

    double d_m, d_l, d_r, d_m2, d_l2, d_r2;
    double  I, I_1 = 0, I_2 = 0;
    double IR, IG, IB;
    int RED = 1, GREEN = 1, BLUE = 1;

	if(x == -1 || x == max_label) return INF;

    if(j - x < 0)              return INF;

    while(RED){

	    d_m = tmp_l2->data[i][j].r;
	    d_m2 =  tmp_r2->data[i][j - x].r;

	    if((j - x - 1 >= 0) && (j - x + 1 <= tmp_r2->width - 1)){

	        d_l2 = (tmp_r2->data[i][j - x].r + tmp_r2->data[i][j - x - 1].r) / 2.0;
	        d_r2 = (tmp_r2->data[i][j - x].r + tmp_r2->data[i][j - x + 1].r) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        if(between(d_r2, d_m, d_m2) == 0)    break;
	        I_1 = min3(fabs(d_m - d_l2), fabs(d_m - d_m2), fabs(d_m - d_r2));              /*特に制約なし*/

	    }else if(j - x - 1 < 0){

	        d_r2 = (tmp_r2->data[i][j - x].r + tmp_r2->data[i][j - x + 1].r) / 2.0;
	        if(between(d_m2, d_m, d_r2) == 0)    break;
	        I_1 = min(fabs(d_m - d_m2), fabs(d_m - d_r2));                                /*Rの左端が出る*/

	    }else{

	        d_l2 = (tmp_r2->data[i][j - x].r + tmp_r2->data[i][j - x - 1].r) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        I_1 = min(fabs(d_m - d_l2), fabs(d_m - d_m2));                                /*Rの右端が出る*/
	    }

	    if((j != 0) && (j != tmp_l2->width - 1)){

	        d_l = (tmp_l2->data[i][j].r + tmp_l2->data[i][j - 1].r) / 2.0;
	        d_r = (tmp_l2->data[i][j].r + tmp_l2->data[i][j + 1].r) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        if(between(d_r, d_m2, d_m) == 0)    break;
	        I_2 = min3(fabs(d_l - d_m2), fabs(d_m - d_m2), fabs(d_r - d_m2));              /*特に制約なし*/

	    }else if(j == 0){

	        d_r = (tmp_l2->data[i][j].r + tmp_l2->data[i][j + 1].r) / 2.0;
	        if(between(d_m, d_m2, d_r) == 0)    break;
	        I_2 = min(fabs(d_m - d_m2), fabs(d_r - d_m2));                                /*Lの左端が出る*/

	    }else{

	        d_l = (tmp_l2->data[i][j].r + tmp_l2->data[i][j - 1].r) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        I_2 = min(fabs(d_l - d_m2), fabs(d_m - d_m2));                                /*Lの右端が出る*/
	    }
	    RED = 0;
	}
    IR = min(I_1, I_2);

    ///////////////////for G of RGB.
    I_1 = 0;
    I_2 = 0;
    while(GREEN){
	    d_m = tmp_l2->data[i][j].g;
	    d_m2 =  tmp_r2->data[i][j - x].g;

	    if((j - x - 1 >= 0) && (j - x + 1 <= tmp_r2->width - 1)){

	        d_l2 = (tmp_r2->data[i][j - x].g + tmp_r2->data[i][j - x - 1].g) / 2.0;
	        d_r2 = (tmp_r2->data[i][j - x].g + tmp_r2->data[i][j - x + 1].g) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        if(between(d_r2, d_m, d_m2) == 0)    break;
	        I_1 = min3(fabs(d_m - d_l2), fabs(d_m - d_m2), fabs(d_m - d_r2));              /*特に制約なし*/

	    }else if(j - x - 1 < 0){

	        d_r2 = (tmp_r2->data[i][j - x].g + tmp_r2->data[i][j - x + 1].g) / 2.0;
	        if(between(d_m2, d_m, d_r2) == 0)    break;
	        I_1 = min(fabs(d_m - d_m2), fabs(d_m - d_r2));                                /*Rの左端が出る*/

	    }else{

	        d_l2 = (tmp_r2->data[i][j - x].g + tmp_r2->data[i][j - x - 1].g) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        I_1 = min(fabs(d_m - d_l2), fabs(d_m - d_m2));                                /*Rの右端が出る*/
	    }

	    if((j != 0) && (j != tmp_l2->width - 1)){

	        d_l = (tmp_l2->data[i][j].g + tmp_l2->data[i][j - 1].g) / 2.0;
	        d_r = (tmp_l2->data[i][j].g + tmp_l2->data[i][j + 1].g) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        if(between(d_r, d_m2, d_m) == 0)    break;
	        I_2 = min3(fabs(d_l - d_m2), fabs(d_m - d_m2), fabs(d_r - d_m2));              /*特に制約なし*/

	    }else if(j == 0){

	        d_r = (tmp_l2->data[i][j].g + tmp_l2->data[i][j + 1].g) / 2.0;
	        if(between(d_m, d_m2, d_r) == 0)    break;
	        I_2 = min(fabs(d_m - d_m2), fabs(d_r - d_m2));                                /*Lの左端が出る*/

	    }else{

	        d_l = (tmp_l2->data[i][j].g + tmp_l2->data[i][j - 1].g) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        I_2 = min(fabs(d_l - d_m2), fabs(d_m - d_m2));                                /*Lの右端が出る*/
	    }
	    GREEN = 0;
	}

    IG = min(I_1, I_2);

       ///////////////////for B of RGB.
    I_1 = 0;
    I_2 = 0;
    while(BLUE){
	    d_m = tmp_l2->data[i][j].b;
	    d_m2 =  tmp_r2->data[i][j - x].b;

	    if((j - x - 1 >= 0) && (j - x + 1 <= tmp_r2->width - 1)){

	        d_l2 = (tmp_r2->data[i][j - x].b + tmp_r2->data[i][j - x - 1].b) / 2.0;
	        d_r2 = (tmp_r2->data[i][j - x].b + tmp_r2->data[i][j - x + 1].b) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        if(between(d_r2, d_m, d_m2) == 0)    break;
	        I_1 = min3(fabs(d_m - d_l2), fabs(d_m - d_m2), fabs(d_m - d_r2));              /*特に制約なし*/

	    }else if(j - x - 1 < 0){

	        d_r2 = (tmp_r2->data[i][j - x].b + tmp_r2->data[i][j - x + 1].b) / 2.0;
	        if(between(d_m2, d_m, d_r2) == 0)    break;
	        I_1 = min(fabs(d_m - d_m2), fabs(d_m - d_r2));                                /*Rの左端が出る*/

	    }else{

	        d_l2 = (tmp_r2->data[i][j - x].b + tmp_r2->data[i][j - x - 1].b) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        I_1 = min(fabs(d_m - d_l2), fabs(d_m - d_m2));                                /*Rの右端が出る*/
	    }

	    if((j != 0) && (j != tmp_l2->width - 1)){

	        d_l = (tmp_l2->data[i][j].b + tmp_l2->data[i][j - 1].b) / 2.0;
	        d_r = (tmp_l2->data[i][j].b + tmp_l2->data[i][j + 1].b) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        if(between(d_r, d_m2, d_m) == 0)    break;
	        I_2 = min3(fabs(d_l - d_m2), fabs(d_m - d_m2), fabs(d_r - d_m2));              /*特に制約なし*/

	    }else if(j == 0){

	        d_r = (tmp_l2->data[i][j].b + tmp_l2->data[i][j + 1].b) / 2.0;
	        if(between(d_m, d_m2, d_r) == 0)    break;
	        I_2 = min(fabs(d_m - d_m2), fabs(d_r - d_m2));                                /*Lの左端が出る*/

	    }else{

	        d_l = (tmp_l2->data[i][j].b + tmp_l2->data[i][j - 1].b) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        I_2 = min(fabs(d_l - d_m2), fabs(d_m - d_m2));                                /*Lの右端が出る*/
	    }
	    BLUE = 0;
	}
    IB = min(I_1, I_2);


    I = IR + IG + IB;

    if(image == 1) return min(IR * IR + IG * IG + IB * IB, 16 * 16);
    else if(image == 2) return min(IR * IR + IG * IG + IB * IB, 16 * 16);
    else if(image == 3) return min(I, 16 * 3);

    //I = min(I, 20);
    //I *= I;
    //return I * I;                                   /*Dを2乗する*/
    //return C_D * (I);

}



double energy(Graph *G, int *label, int *I, double T, int lambda) {
    int i;
    double energy = 0;
    // Dterm
    for (i = 1; i <= G->n - 2; i++) {
        energy += data(I, i, label[i]);
    }
    // Vterm
    for (i = 1; i <= G->m - 2 * (G->n - 2); i++) {
        energy += pairwise(label[G->tail[i]] , label[G->head[i]], T, lambda);
    }
    return energy;
}


int update_labels(Graph *G, int alpha, int *label, int *t, int *I, double T, int lambda) {
    int i, res;
    int *temp;
    double prev_energy, new_energy;
    if ((temp = (int *) malloc(sizeof(int) * (G->n - 1))) == NULL) {
        fprintf(stderr, "update_labels(): ERROR [temp = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // alpha-expansion
    for (i = 1; i <= G->n - 2; i++) {
        if (t[i] == 1) temp[i] = alpha;
        else temp[i] = label[i];
    }

    // binary
    // for (i = 1; i <= G->n - 2; i++) {
    //     temp[i] = t[i];
    // }

    prev_energy = energy(G, label, I, T, lambda);
    new_energy = energy(G, temp, I, T, lambda);

    if (new_energy < prev_energy) {
        for (i = 1; i <= G->n - 2; i++) label[i] = temp[i];
        printf("Energy %lf ---> %lf\n", prev_energy, new_energy);
        res = 1;
    } else if (prev_energy == new_energy) {
        res = 0;
    } else {
        fprintf(stderr, "【 WARNING 】 energy(G, newLabel, I) > energy(G, currentLabel, I)\n");
        printf("【 WARNING 】 Energy %lf ---> %lf\n", prev_energy, new_energy);
        res = -1;
    }
    free(temp);
    return res;
}
void set_capacity(Graph *G, int alpha, int *label, int *I, double T, int lambda) {
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
        //labelがalphaでない->capacity設定
        //labelがalpha->枝を削除
        if (label[G->head[i]] == alpha) {
            deleteAdjEdge(G, i);
            deleteAdjEdge(G, i + grids_node);
        }
    }
    // set Dterm
     for(i = 1; i < G->n - 1; i++){
        if(pairwise(I[i],label[i], T, lambda) > pairwise(I[i],alpha, T, lambda)){
            G->capa[s2i_begin -1 + i] +=  (data(I, i, label[i]) - data(I, i, alpha));
        }else{
            G->capa[i2t_begin - 1  + i] += (data(I, i,alpha) - data(I, i, label[i]));
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
