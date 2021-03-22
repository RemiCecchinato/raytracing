#ifndef MAIN_H
#define MAIN_H

#include "math.h"
#include "fast_obj.h"

typedef struct Scene
{
    Camera camera;

    int32_t  material_count;
    Material *materials;

    int32_t sphere_count;
    Sphere  *spheres;

    int32_t plane_count;
    Plane   *planes;

    int32_t light_count;
    Light  *lights;

    int32_t mesh_count;
    Scene_Mesh *meshes;
} Scene;

typedef struct Ray
{
    Vec3f origin;
    Vec3f direction;
} Ray;

typedef enum Object_Type
{
    SPHERE,
    PLANE,
    LIGHT,
    TRIANGLE,
} Object_Type;

typedef struct Ray_Result
{
    bool hit;
    Object_Type object_type;
    uint32_t object_id;
    uint32_t mesh_id;

    // infos géométriques
    Vec3f intersection_point;
    Vec3f surface_normal;
    float t_min;

    float triangle_alpha;
    float triangle_beta;
    float triangle_gamma;
    Vec2f triangle_uv;


    bool enter_shape;

    // infos graphiques
    int32_t material_id;
} Ray_Result;

typedef struct Job
{
    Scene   *scene;
    Image3f *image;
    Rect2i   region;
    uint32_t ray_per_pixel;
    Random_Serie serie;
} Job;

#define ARRAY_SIZE(x) ((sizeof(x) / sizeof(x[0])))

#define THREAD_COUNT 15

#define IMAGE_COUNT 1

#define IMAGE_SIZE 1024
#define RAY_PER_PIXEL 10000

#endif