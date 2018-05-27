#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

const double Pi = acos(-1.0);

struct vec{
    double x, y ,z;
    vec(double x_ = 0, double y_ = 0, double z_ = 0): x(x_), y(y_), z(z_) {}
    vec operator+(const vec &b) const { return vec(x + b.x, y + b.y, z + b.z);}
    vec operator-(const vec &b) const { return vec(x - b.x, y - b.y, z - b.z);}
    vec operator*(double b) const { return vec(x*b, y*b, z*b);}
    vec mult(const vec &b) const { return vec(x*b.x, y*b.y, z*b.z);}
    double dot(const vec &b) const { return x * b.x + y * b.y + z * b.z;}
    double dis(){ return sqrt(x*x + y*y + z*z);}
    vec norm(){ return (*this)*(1.0/this->dis());}
};

enum Refl_t {mirror, highlgt, normal};

inline vec cross(vec r, vec s){
    return vec(r.y*s.z - r.z*s.y, r.z*s.x - r.x*s.z, r.x*s.y - r.y*s.x);
}

inline double clamp(double x){
    return x < 0 ? 0 : x > 1 ? 1 : x;
}

inline int toInt(double x){
    return int(pow(clamp(x), 1/2.2) * 255 + 0.5);
}

struct ray{
    vec e, d;
    ray(vec e_, vec d_): e(e_), d(d_) {}
};

struct Sphere{
    double radius;
    vec p, c, e;//position, color, emission
    Refl_t refl;
    Sphere(double radius_, vec p_, vec c_, vec e_, Refl_t refl_): radius(radius_), p(p_), c(c_), e(e_), refl(refl_) {}
    double intersect(const ray &r){
        double eps = 1e-4;
        vec tmp = r.e-p;
        double delta = pow(r.d.dot(tmp), 2) - tmp.dot(tmp) + radius * radius;
        if(delta < eps) return 0;
        double lft = -r.d.dot(tmp);
        double rgt = sqrt(delta);
        double t;
        double result = (t=lft-rgt) > eps ? t : (t=(lft+rgt)) > eps ? t : 0;
        return result;
    }
};

struct Rect{
    vec a,b,c,n,color;// (a-b)垂直于(c-b)
    Rect(vec a_, vec b_, vec c_, vec color_): a(a_), b(b_), c(c_), color(color_) {
        vec r = c-b, s = a-b;
        n = cross(r,s).norm();
    }
    double intersect(const ray& r){
        double t = n.dot(a-r.e) / n.dot(r.d);
        vec p = r.e + r.d * t;
        double alpha = (a-b).norm().dot(p-b), beta = (c-b).norm().dot(p-b);
        if(t > 1e-4 && alpha > 0 && alpha < (a-b).dis() && beta > 0 && beta < (c-b).dis()){
            return t;
        }
        return 0;
    }
};

Rect rects[]={
    Rect(vec(0,0,100), vec(0,0,0), vec(0,100,0), vec(0.75,0.25,0.25)),//left
    Rect(vec(100,100,0), vec(100,0,0), vec(100,0,100), vec(0.25,0.25,0.75)),//right
    Rect(vec(0,100,0), vec(100,100,0), vec(100,100,100), vec(0,0,0)),//front
    Rect(vec(0,0,0), vec(0,0,100), vec(100,0,100), vec(0.75,0.75,0.75)),//back
    Rect(vec(100,100,100), vec(100,0,100), vec(0,0,100), vec(0.75,0.75,0.75)),//top
    Rect(vec(100,100,0), vec(0,100,0), vec(0,0,0),vec(0.75,0.75,0.75))//bottom
};

Sphere spheres[] = {
    Sphere(15, vec(35,30,15), vec(0.99,0.99,0.99), vec(0,0,0), mirror),
    //Sphere(15, vec(70,45,15), vec(0.4,0.4,0.7), vec(0,0,0), highlgt),
    Sphere(30, vec(50,35,128), vec(1,1,1), vec(50,50,50), normal)
};

inline bool intersect(const ray &r, int &rect_id, int &sphere_id, double &t){
    rect_id = sphere_id = -1;
    t = 1e20;
    double tmp;
    for(size_t i = 0; i < sizeof(rects) / sizeof(Rect); i++){
        tmp = rects[i].intersect(r);
        if(tmp > 1e-4 && tmp < t){
            t = tmp;
            rect_id = i;
        }
    }
    for(size_t i =0; i < sizeof(spheres) / sizeof(Sphere); i++){
        tmp = spheres[i].intersect(r);
        if(tmp > 1e-4 && tmp < t){
            t = tmp;
            sphere_id = i;
        }
    }
    if(rect_id != -1 || sphere_id != -1)
        return true;
    return false;
}

int m = 16, maxdepth = 10;
double cl = 0.8, eps = 1e-4;
vec ca;

int rect_id, sphere_id;

vec radiance(ray r, unsigned short int *xi, int depth){

    //printf("e3 %d\n", depth);
    //getchar();

    if(depth > maxdepth) return vec();

    double t;
    if(intersect(r, rect_id, sphere_id, t) == 0) return vec();

    //printf("e4 %d\n",rect_id);
    //printf("e8 %d\n",sphere_id);

    vec p = r.e + r.d * t;
    if(sphere_id != -1){
        Sphere& obj = spheres[sphere_id];

        vec n = (p - obj.p).norm();
        if(obj.refl == normal){
            return obj.e;
        }
        /*else if(obj.refl == highlgt){//high light
            vec h = (r.d + l).norm();
            return vec(1,1,1) * cl * pow(h.mul(n), m);
        }*/
        else {//mirror
            vec lgt = (r.d - n * (n.dot(r.d) * 2)).norm();
            return obj.c.mult(radiance(ray(p, lgt), xi, depth + 1));
        }
    }
    else {
        Rect obj = rects[rect_id];
        double phi = 2 * Pi * erand48(xi), dist_2 = erand48(xi), dist = sqrt(dist_2);

        //printf("e5 %lf %lf %lf %lf\n", Pi, phi, dist_2, dist);

        double ul = sin(phi) * dist, vl = cos(phi) * dist, wl = sqrt(1 - dist_2);
        vec w = obj.n, u = cross(w.x == 0 ? vec(1,0,0): vec(0,1,0), w), v = cross(u,w);

        /*printf("e7 %lf %lf %lf\n", w.x, w.y, w.z);
        printf("e7 %lf %lf %lf\n", u.x, u.y, u.z);
        printf("e7 %lf %lf %lf\n", v.x, v.y, v.z);*/

        ray R = ray(p, (w * wl + u * ul + v * vl).norm());

        //printf("e6 %lf %lf %lf\n%lf %lf %lf\n", R.e.x, R.e.y, R.e.z, R.d.x, R.d.y, R.d.z);

        return obj.color.mult(radiance(R, xi, depth + 1));
    }
}

int main()
{
    int wi = 1024, h = 768, samps = 500;
    vec *c = new vec[wi*h], e = vec(50,100-(1e-4),15), g = vec(0,-1,0), up = vec(0,0,1), s, p;
    vec w = (vec()-g).norm(), u = cross(up, w).norm(), v = cross(w, u);

    double su, sv;

    #pragma omp parallel for schedule(dynamic, 1) private(su, sv)

    for(int j = 0; j < h; j++){
        if(j % 100 == 0) printf("%lf%%\n", j * 1.0 / h);
        for(int i = 0; i < wi; i++){
            su = 1 - (i + 0.5) * 2.0 / wi;
            sv = 1.5 - (j + 0.5) * 2.0 / h;
            s = e - w + u * su + v * sv;
            ray r = ray(e, (s-e).norm());
            unsigned short int seed[3] = {0,0,(unsigned short)time(NULL)};
            vec rad = vec();
            for(int s = 0; s < samps; s++){
                rad = rad + radiance(r, seed, 0);
                //printf("e1 %lf %lf %lf\n", rad.x, rad.y, rad.z);
                //getchar();
            }
            rad = rad * (1.0 / samps);

            //printf("e2 %lf %lf %lf\n", rad.x, rad.y, rad.z);
            //getchar();

            c[j*wi+i] = c[j*wi+i] + vec(clamp(rad.x), clamp(rad.y), clamp(rad.z));
        }
    }

    FILE *fp = fopen("image1.ppm","w");
    fprintf(fp, "P3\n%d %d\n%d\n", wi, h, 255);
    for(int i=0; i<wi*h; i++){
        fprintf(fp, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }
    return 0;
}
