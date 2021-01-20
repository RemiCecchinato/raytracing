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
    Vec3f *pixels;
} Image3f;

typedef struct Light
{
    Vec3f position;
    Vec3f color;
    float intensity;
} Light;

typedef struct Material
{
    // reflection information
    Vec3f diffuse_color;
    Vec3f mirror_color;

    // transparancy information
    bool transparent;
    float n;
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

const Vec3f BLACK = {0.0f, 0.0f, 0.0f};
const Vec3f WHITE = {1.0f, 1.0f, 1.0f};
const Vec3f RED   = {1.0f, 0.0f, 0.0f};
const Vec3f GREEN = {0.0f, 1.0f, 0.0f};
const Vec3f BLUE  = {0.0f, 0.0f, 1.0f};
const Vec3f PINK  = {1.0f, 0.0f, 1.0f};

#endif