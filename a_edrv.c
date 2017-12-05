#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bmp.h"
#include "a-expan.h"
#include "graph.h"
#include "ford_fulkerson.h"

#define _CHK_ENERGY_ 0      // ラベル更新時にエネルギーの減少を確認す　0:確認しない 1:確認する
#define _OUTPUT_INFO_ 0     // デバッグ情報出力 0:出力しない 1:出力する
#define _OUTPUT_GRAPH_ 0    // グラフ情報出力 　0:出力しない 1:出力する
#define _OUTPUT_PROGRESS_ 0 // 処理過程ファイル出力 　0:出力しない 1:出力する

int main(int argc, char *argv[]) {
    long int i, j, node, edge;
    int scale = 16, alpha, label_size, last, iteration;
    int *I, *t, *label;
    int lambda = 1;
    // I->入力画像の輝度, t->2値変数, label->ラベル付け
    double flag, temp, T = 9;
    double *f;
    char output_file[100];
    clock_t start;
    img image, output;
    Graph G;

#if _OUTPUT_INFO_
    double maxflow;
#endif

#if _OUTPUT_PROGRESS_
    int k = 0;
    char pf[100];
    system("rm output/*.bmp &> /dev/null");
#endif


    label_size = 255 / scale;

    if (argc != 2 && argc != 3) {
        printf("Usage: %s <input_file> <output_file(option)>\n", argv[0]);
        return 1;
    }

    if (argc == 2) strcpy(output_file, "/dev/null");
    else strcpy(output_file, argv[2]);

    ReadBmp(argv[1], &image);
    ReadBmp(argv[1], &output);

    if ((I = (int *)malloc(sizeof(int) * (image.height * image.width + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->I]\n");
        exit(EXIT_FAILURE);
    }

    printf("height %ld, width %ld\n", image.height, image.width);

    for (i = 0; i <  image.height; i++) {
        for (j = 0; j < image.width; j++) {
            I[i * image.width + j + 1] = image.data[i][j].r / scale;
        }
    }

    node = image.height * image.width + 2;
    edge = (image.height - 1) * image.width + image.height * (image.width - 1) + 2 * image.height * image.width;

    // グラフ初期設定
    newGraph(&G, node, edge);
    set_all_edge(&G, image.height, image.width);
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
    for (i = 1; i <= G.n - 2; i++) label[i] = I[i];
    printf("Energy (before): %lf\n", energy(&G, label, I, T, lambda));

    start = clock();
    last = label_size;
    flag = 1;
    iteration = 0;
    while (flag > 0) {
        temp = energy(&G, label, I, T, lambda);
        for (alpha = 0; alpha <= label_size; alpha++) {
            if(last == alpha) break;
            for (i = 0; i <= G.m + 1; i++) {
                f[i] = 0;
                G.capa[i] = 0;
            }

            // capacity設定
            set_capacity(&G, alpha, label, I, T, lambda);
            // capacity(&G, label, I, alpha);
#if _OUTPUT_INFO_
            printf("Maximum Flow: %5.2lf\n", boykov_kolmogorov(G, f, t));
            printf("after bk-max-f\n");
            showAdjList(&G);
#else
            boykov_kolmogorov(G, f, t);
#endif
            //削除された枝はvoykov_kolmogorov(G, f, t)で再度リストに追加されている
            iteration++;

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
        flag = temp - energy(&G, label, I, T, lambda);
        // printf("Energy %lf\n", energy(&G, label, I, T, lambda));
    }

    printf("Energy (after): %lf\n", energy(&G, label, I, T, lambda));
    printf("Run time[%.2lf]\n", (double) (clock() - start) / CLOCKS_PER_SEC);
    printf("iteration: %d\n", iteration);

#if _OUTPUT_GRAPH_
    printf("---Graph information---\n");
    showGraph(&G);
#endif

    // output to bitmap file
    for (i = 0; i <  image.height; i++) {
        for (j = 0; j < image.width; j++) {
            output.data[i][j].r = label[i * image.width + j + 1] * scale;
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

    delGraph(&G);
    free(I);
    free(f);
    free(t);
    free(label);

    return 0;
}
