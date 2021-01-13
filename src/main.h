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

typedef struct Ray
{
    Vec3f origin;
    Vec3f direction;

    Vec3f color;
    uint32_t hit_count;
} Ray;

typedef struct Ray_Result
{
    bool hit;

    // infos géométriques
    Vec3f intersection_point;
    Vec3f surface_normal;
    float t_min;

    // infos graphiques
    int32_t material_id;
} Ray_Result;

typedef struct Job
{
    Scene   *scene;
    Image3f *image;
    Rect2i   region;
} Job;

#define ARRAY_SIZE(x) ((sizeof(x) / sizeof(x[0])))

#endif