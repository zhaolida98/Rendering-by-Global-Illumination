#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -origin smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2

struct Vec          // Usage: time ./smallpt 5000 && xv image.ppm
{
    double x, y, z;                  // position, also color (r,g,b)
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) {
        x = x_;
        y = y_;
        z = z_;
    }

    Vec operator+(const Vec &b) const {
        return Vec(x + b.x, y + b.y, z + b.z);
    }

    Vec operator-(const Vec &b) const {
        return Vec(x - b.x, y - b.y, z - b.z);
    }

    Vec operator*(double b) const {
        return Vec(x * b, y * b, z * b);
    }

    Vec mult(const Vec &b) const {
        return Vec(x * b.x, y * b.y, z * b.z);
    }

    Vec &norm() {
        return *this = *this * (1 / sqrt(x * x + y * y + z * z));
    }

    double dot(const Vec &b) const {
        return x * b.x + y * b.y + z * b.z;    // cross:
    }

    Vec operator%(Vec &b) {
        return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};

struct Ray {
    Vec origin, direction;

    Ray(Vec o_, Vec d_) : origin(o_), direction(d_) {}
};

enum Refl_t {
    DIFF, SPEC, REFR
};  // material types, used in radiance()
struct Sphere {
    double rad;       // radius
    Vec position, emission, color;      // position, emission, color
    Refl_t reflection_type;      // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) :
            rad(rad_), position(p_), emission(e_), color(c_), reflection_type(refl_) {}

    double intersect(const Ray &r) const   // returns distance, 0 if nohit
    {
        Vec op = position - r.origin; // op is vector from ray orignal point to sphere's center
        // Solve t^2*direction.direction + 2*t*(origin-position).direction + (origin-position).(origin-position)-R^2 = 0
        double t;
        double eps = 1e-4;// threshold for distance
        double b = op.dot(
                r.direction); // r.direction is norm, so op.dot(r.direction) is |op|*cos(theta), the length of op's projection in the ray direction
        double det = b * b - op.dot(op) + rad * rad; // Pythagorean theorem
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t
                                                              : 0); // two intersection, choose one according to the threshold
    }
};

Sphere spheres[] =  //Scene: radius, position, emission, color, material
        {
                Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),//Left
                Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF),//Rght
                Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),//Back
                Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),//Frnt
                Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),//Botm
                Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF),//Top
                Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),//Mirr
                Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR),//Glas
                Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF) //Lite
        };

inline double clamp(double x) {
    return x < 0 ? 0 : x > 1 ? 1 : x; // merge to 0~1
}

inline int toInt(double x) {
    return int(pow(clamp(x), 1 / 2.2) * 255 + .5); //transfer 0~1 to color int
}

inline bool intersect(const Ray &r, double &t, int &id) {
    //find the nearlest spheres intersecting with ray
    double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
    for (int i = int(n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t) {
            t = d;
            id = i;
        }
    return t < inf;
}

Vec radiance(const Ray &ray, int depth, unsigned short *Xi) {
    double d_rayo_insec;                               // distance to intersection
    int id = 0;                               // id of intersected object
    if (!intersect(ray, d_rayo_insec, id))
        return Vec(); // if miss, return black
    const Sphere &obj = spheres[id];        // the hit object
    Vec x = ray.origin + ray.direction * d_rayo_insec; //intersection point to zero point
    Vec n = (x - obj.position).norm();
    Vec nl = n.dot(ray.direction) < 0 ? n : n * -1, f = obj.color;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max reflection_type color
    if (++depth > 5)
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return obj.emission; //R.R.
    if (obj.reflection_type == DIFF)                   // Ideal DIFFUSE reflection
    {
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.emission + f.mult(radiance(Ray(x, d), depth, Xi));
    } else if (obj.reflection_type == SPEC)              // Ideal SPECULAR reflection
        return obj.emission + f.mult(radiance(Ray(x, ray.direction - n * 2 * n.dot(ray.direction)), depth, Xi));
    Ray reflRay(x, ray.direction - n * 2 * n.dot(ray.direction));     // Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;                // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = ray.direction.dot(nl), cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)    // Total internal reflection
        return obj.emission + f.mult(radiance(reflRay, depth, Xi));
    Vec tdir = (ray.direction * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return obj.emission + f.mult(depth > 2 ? (erand48(Xi) < P ?   // Russian roulette
                                              radiance(reflRay, depth, Xi) * RP : radiance(Ray(x, tdir), depth, Xi) *
                                                                                  TP) :
                                 radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
}

int main(int argc, char *argv[]) {
    int width = 1024, height = 768, samps_num = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
    Vec x_incr = Vec(width * .5135 / height);
    Vec y_incr = (x_incr % cam.direction).norm() * .5135;  // .5135 camera angle
    Vec r;
    Vec *c = new Vec[width * height];
    #pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for (int y = 0; y < height; y++)                        // Loop over image rows
    {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps_num * 4, 100. * y / (height - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < width; x++)  // Loop cols
            for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++)     // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())         // 2x2 subpixel cols
                {
                    for (int s = 0; s < samps_num; s++) {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = x_incr * (((sx + .5 + dx) / 2 + x) / width - .5) +
                                y_incr * (((sy + .5 + dy) / 2 + y) / height - .5) + cam.direction;
                        r = r + radiance(Ray(cam.origin + d * 140, d.norm()), 0, Xi) * (1. / samps_num);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
    }
    FILE *f = fopen("image.ppm", "width");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
    for (int i = 0; i < width * height; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
