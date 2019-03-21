#include <stdio.h>
#include <math.h>

bool intersect(const Ray &r, double &t, int &id);
int toInt(double x);
double clamp(double x);
Vec radiance(const Ray &r, int depth, unsigned short *Xi);
void sample(unsigned short x, int y, Vec cx, Vec cy, Vec r, Vec *c,int h,int w, unsigned short Xi[], Ray cam);