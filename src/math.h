#ifndef MATH_H
#define MATH_H

#include <stdint.h>
#include <immintrin.h>
#include <wmmintrin.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

typedef struct Random_Serie
{
    union
    {
        uint8_t bytes[32];
        __m128i xmm[2];
    };
} Random_Serie;

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

typedef union Vec4f
{
    struct {
        float x, y, z, w;
    };

    struct {
        float coord[4];
    };
} Vec4f;

typedef struct Vertex
{
    uint32_t vertexId;
    uint32_t texcoordId;
    uint32_t normalId;
} Vertex;

typedef struct Triangle
{
    Vertex ids[3];
} Triangle;

typedef struct Bounding_Box
{
    Vec3f min;
    Vec3f max;
} Bounding_Box;

typedef struct Mesh_Node
{
    Bounding_Box bbox;

    struct Mesh_Node *childs;

    uint32_t first_triangle_id;
    uint32_t last_triangle_id;
} Mesh_Node;

typedef struct Mesh
{
    uint32_t point_count;
    Vec3f   *points;

    uint32_t texcoord_count;
    Vec2f   *texcoords;

    uint32_t normal_count;
    Vec3f   *normals;

    uint32_t  face_count;
    Triangle *faces;

    uint32_t tree_depth;
    Mesh_Node root;
} Mesh;

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

typedef struct Material
{
    // reflection information
    Vec3f diffuse_color;

    bool mirror;
    Vec3f mirror_color;

    // transparancy information
    bool transparent;
    float n;

    // texture information
    bool textured;
    Image3f *texture;
} Material;

typedef struct Sphere
{
    Vec3f center;
    float radius;
    int32_t material_id;
} Sphere;

typedef struct Light
{
    Vec3f center;
    float radius;
    Vec3f color;
    float albedo;
} Light;

typedef struct Plane
{
    Vec3f normal;
    float d;
    int32_t material_id;
} Plane;

const Vec3f BLACK  = {0.0f, 0.0f, 0.0f};
const Vec3f WHITE  = {1.0f, 1.0f, 1.0f};
const Vec3f RED    = {1.0f, 0.0f, 0.0f};
const Vec3f GREEN  = {0.0f, 1.0f, 0.0f};
const Vec3f BLUE   = {0.0f, 0.0f, 1.0f};
const Vec3f PINK   = {1.0f, 0.0f, 1.0f};
const Vec3f YELLOW = {1.0f, 1.0f, 0.0f};

#endif