#include "bmp.h"
#include "a_estr.h"
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

double energy(Graph *G, int *label,  double T, int lambda, int width, int *I_left, int *I_right) {
    int i;
    double energy = 0;
    //* Dterm
    for (i = 1; i <= G->n - 2; i++) {
        energy += data(i, label[i], width, I_left, I_right);
        // energy += data(I_left[i], label[i]);
    }
    // */
    // Vterm
    for (i = 1; i <= G->m - 2 * (G->n - 2); i++) {
        energy += pairwise(label[G->tail[i]], label[G->head[i]], T, lambda);
    }
    return energy;
}

int update_labels(Graph *G, int *label, int alpha, int width, double T, int lambda, int *t, int *I_left, int *I_right) {
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

    prev_energy = energy(G, label, width, T, lambda, I_left, I_right);
    new_energy = energy(G, temp, width, T, lambda, I_left, I_right);

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
void set_capacity(Graph *G, int *label, int width, int alpha, double T, int lambda, int *I_left, int *I_right) {
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
        if(data(i, label[i], width, I_left, I_right) > data(i, alpha, width, I_left, I_right)){
            G->capa[s2i_begin - 1 + i] +=  (data(i, label[i], width, I_left, I_right) - data(i, alpha, width, I_left, I_right));
        }else{
            G->capa[i2t_begin - 1  + i] += (data(i, alpha, width, I_left, I_right) - data(i, label[i], width, I_left, I_right));
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