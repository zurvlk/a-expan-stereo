#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bmp.h"
#include "a_estr.h"
#include "graph.h"
#include "ford_fulkerson.h"

#define _CHK_ENERGY_ 0      // ラベル更新時にエネルギーの減少を確認す　0:確認しない 1:確認する
#define _OUTPUT_INFO_ 0     // デバッグ情報出力 0:出力しない 1:出力する
#define _OUTPUT_GRAPH_ 0    // グラフ情報出力 　0:出力しない 1:出力する
#define _OUTPUT_PROGRESS_ 0 // 処理過程ファイル出力 　0:出力しない 1:出力する

int main(int argc, char *argv[]) {
    long int i, j, node, edge, grids_node;
    int scale = 16, alpha, label_size, label_max, last;
    int *I_left, *I_right, *t, *label;
    int lambda = 1;
    // I->入力画像の輝度, t->2値変数, label->ラベル付け
    double flag, temp, error, T = 9;
    double *f;
    char output_file[100];
    clock_t start;
    img left, right, truth, output;
    Graph G;

    char imgleft[100];
    char imgright[100];
    char imgtruth[100];

#if _OUTPUT_INFO_
    double maxflow;
#endif

#if _OUTPUT_PROGRESS_
    int k = 0;
    char pf[100];
    system("rm output/*.bmp &> /dev/null");
#endif

    if (argc != 2 && argc != 3 && argc != 5 && argc != 6) {
        printf("Usage: %s <input_file> <output_file(option)> <scale(option)> <lambda (option)> <T (oprion)>\n", argv[0]);
        return 1;
    }

    if (argc == 2) strcpy(output_file, "/dev/null");
    else strcpy(output_file, argv[2]);
    if (argc == 5 || argc == 6) {
        scale = atoi(argv[3]);
        lambda = atoi(argv[4]);
        if(argc == 6) T = atof(argv[5]);
    }


    label_size = 256 / scale;
    label_max = label_size - 1;


    label_size = 255 / scale;

    if (argc == 2) strcpy(output_file, "/dev/null");
    else strcpy(output_file, argv[2]);

   printf("----------------------------------------------\n");
    printf("input_file: %s\n", argv[1]);
    printf("output_file: %s\n", output_file);
    printf("label_size: %d\n", label_size);
    printf("lambda: %d\n", lambda);
    printf("T: %.2f\n", T);

    strcpy(imgleft, argv[1]);
    strcpy(imgright, argv[1]);
    strcpy(imgtruth, argv[1]);

    strcat(imgleft, "left.bmp");
    strcat(imgright, "right.bmp");
    strcat(imgtruth, "truth.bmp");

    ReadBmp(imgleft, &left);
    ReadBmp(imgright, &right);
    ReadBmp(imgtruth, &truth);
    ReadBmp(imgtruth, &output);

    Gray(&left, &left);
    Gray(&right, &right);
    Gray(&truth, &truth);

    grids_node = left.height * left.width;

    if ((I_left = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->I_left]\n");
        exit(EXIT_FAILURE);
    }
    if ((I_right = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->I_right]\n");
        exit(EXIT_FAILURE);
    }

    printf("height %ld, width %ld\n", left.height, left.width);
    for (i = 0; i <  left.height; i++) {
        for (j = 0; j < left.width; j++) {
            I_left[i * left.width + j + 1] = left.data[i][j].r / scale;
        }
    }
    for (i = 0; i <  right.height; i++) {
        for (j = 0; j < right.width; j++) {
            I_right[i * right.width + j + 1] = right.data[i][j].r / scale;
        }
    }

    node = grids_node + 2;
    edge = (left.height - 1) * left.width + left.height * (left.width - 1) + 2 * grids_node;

    // グラフ初期設定
    newGraph(&G, node, edge);
    set_all_edge(&G, left.height, left.width);
    initAdjList(&G);

    if ((f = (double *) malloc(sizeof(double) * (G.m + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [f = malloc()]\n");
        return (EXIT_FAILURE);
    }
    if ((t = (int *) malloc(sizeof(int) * (G.n + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [t = malloc()]\n");
        return (EXIT_FAILURE);
    }
    if ((label = (int *) malloc(sizeof(int) * (G.n - 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [label = malloc()]\n");
        return (EXIT_FAILURE);
    }
    // 輝度から初期ラベル設定
    for (i = 1; i <= G.n - 2; i++) label[i] = 0;
    printf("Energy (before): %lf\n", energy(&G, label, T, lambda, left.width, I_left, I_right));

    start = clock();
    last = label_size;
    flag = 1;
    while (flag > 0) {
        temp = energy(&G, label, T, lambda, left.width, I_left, I_right);
        for (alpha = 0; alpha <= label_size; alpha++) {
            if(last == alpha) break;
            for (i = 0; i <= G.m + 1; i++) {
                f[i] = 0;
                G.capa[i] = 0;
            }

            // capacity設定
            set_capacity(&G, label, left.width, alpha, T, lambda, I_left, I_right);
            // capacity(&G, label, I, alpha);

#if _OUTPUT_INFO_
            printf("Maximum Flow: %5.2lf\n", boykov_kolmogorov(G, f, t));
            printf("after bk-max-f\n");
            showAdjList(&G);
#else
            boykov_kolmogorov(G, f, t);
#endif
            //削除された枝はvoykov_kolmogorov(G, f, t)で再度リストに追加されている

            // tを基にラベル更新
#if _CHK_ENERGY_
            if(update_labels(&G, alpha, label, t, I)) last = alpha;
#else
            for (i = 1; i <= G.n - 2; i++) {
                if(t[i] == 1) {
                    label[i] = alpha;
                    last = alpha;
                }
            }
#endif

#if _OUTPUT_PROGRESS_
            for (i = 0; i <  image.height; i++) {
                for (j = 0; j < image.width; j++) {
                    output.data[i][j].r = label[i * image.width + j + 1] * scale;
                    output.data[i][j].g = output.data[i][j].r;
                    output.data[i][j].b = output.data[i][j].r;
                }
            }
            sprintf(pf, "output/image_%04d.bmp", k);
            WriteBmp(pf, &output);
            k++;
#endif
        }
        flag = temp - energy(&G, label, T, lambda, left.width, I_left, I_right);
        printf("Energy %lf\n", energy(&G, label, T, lambda, left.width, I_left, I_right));
    }

    printf("Energy (after): %lf\n", energy(&G, label, T, lambda, left.width, I_left, I_right));
    printf("Run time[%.2lf]\n", (double) (clock() - start) / CLOCKS_PER_SEC);

#if _OUTPUT_GRAPH_
    printf("---Graph information---\n");
    showGraph(&G);
#endif

    // output to bitmap file
    for (i = 0; i <  left.height; i++) {
        for (j = 0; j < left.width; j++) {
            output.data[i][j].r = label[i * left.width + j + 1] * scale;
            output.data[i][j].g = output.data[i][j].r;
            output.data[i][j].b = output.data[i][j].r;
        }
    }
    WriteBmp(output_file, &output);

#if _OUTPUT_PROGRESS_
    if (system("cd output && yes |./mkmovie.sh") == -1) {
      printf("動画の出力に失敗しました\n");
    }
#endif

    if (strcmp(output_file, "/dev/null") != 0){
        ReadBmp(output_file, &output);
        // Gray(&truth, &truth);
        Gray(&output, &output);
        error = err_rate(output, truth, scale);
        printf("Error rate : %lf\n", error);
    }
    delGraph(&G);
    free(I_left);
    free(I_right);
    free(f);
    free(t);
    free(label);

    return 0;
}
