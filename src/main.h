#ifndef MAIN_H
#define MAIN_H

#include "math.h"

typedef struct Scene
{
    Camera camera;

    int32_t  material_count;
    Material *materials;

    int32_t sphere_count;
    Sphere  *spheres;

    int32_t plane_count;
    Plane   *planes;
} Scene;

typedef struct Job
{
    Scene   *scene;
    Image3f *image;
    Rect2i   region;
} Job;

#define ARRAY_SIZE(x) ((sizeof(x) / sizeof(x[0])))

#endif