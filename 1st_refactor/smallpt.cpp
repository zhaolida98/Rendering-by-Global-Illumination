#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include "struct.h"
#include "function.h"
void sample(unsigned short x, int y, Vec cx, Vec cy, Vec r, Vec *c,int h,int w, unsigned short Xi[], Ray cam);
int main(int argc, char *argv[]) {
    int w = 1280, h = 960; // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r, *c = new Vec[w * h];
    #pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for (int y = 0; y < h; y++) {                       // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%",4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, (unsigned short)pow(y,3)}; x < w; x++)   // Loop cols
            sample(x,y, cx, cy, r, c, h, w, Xi, cam);
    }
    FILE *f = fopen("image_100.ppm", "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
