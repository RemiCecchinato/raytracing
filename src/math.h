#ifndef MATH_H
#define MATH_H

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

typedef union Vec2i
{
    struct {
        uint32_t x, y;
    };

    struct {
        uint32_t width, height;
    };

    struct {
        uint32_t coord[2];
    };
} Vec2i;

typedef union Vec2f
{
    struct {
        float x, y;
    };

    struct {
        float coord[2];
    };
} Vec2f;

typedef struct Rect2i
{
    Vec2i point;
    Vec2i size;
} Rect2i;

typedef union Vec3f
{
    struct {
        float x, y, z;
    };

    struct {
        float r, g, b;
    };

    struct {
        float coord[3];
    };
} Vec3f;

typedef struct Camera
{
    Vec3f position;
    Vec3f lookAt;
    Vec3f up;
    float field_of_view; // in degrees
} Camera;

typedef struct Image3f
{
    Vec2i  size;
    float *r;
    float *g;
    float *b;
} Image3f;

typedef struct Material
{
    Vec3f color;
} Material;

typedef struct Sphere
{
    Vec3f center;
    float radius;
    int32_t material_id;
} Sphere;

typedef struct Plane
{
    Vec3f normal;
    float d;
    int32_t material_id;
} Plane;

#endif