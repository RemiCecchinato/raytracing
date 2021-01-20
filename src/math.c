#include <math.h>

float dot3f(Vec3f left, Vec3f right)
{
    float result = 0.0f;

    for (int i = 0; i < 3; i++) {
        result += left.coord[i] * right.coord[i];
    }

    return result;
}

float norm2(Vec3f vec)
{
    float result = dot3f(vec, vec);

    return result;
}

float norm(Vec3f vec)
{
    float vec_norm2 = norm2(vec);
    float result = sqrtf(vec_norm2);

    return result;
}

Vec3f normalize(Vec3f vec)
{
    float vec_norm = norm(vec);

    Vec3f result = {
        .x = vec.x / vec_norm,
        .y = vec.y / vec_norm,
        .z = vec.z / vec_norm,
    };

    return result;
}

Vec3f add3f(Vec3f left, Vec3f right)
{
    Vec3f result = {
        .x = left.x + right.x,
        .y = left.y + right.y,
        .z = left.z + right.z,
    };

    return result;
}

Vec3f sub3f(Vec3f left, Vec3f right)
{
    Vec3f result = {
        .x = left.x - right.x,
        .y = left.y - right.y,
        .z = left.z - right.z,
    };

    return result;
}

Vec3f mul3f(Vec3f left, Vec3f right)
{
    Vec3f result = {
        .x = left.x * right.x,
        .y = left.y * right.y,
        .z = left.z * right.z,
    };

    return result;
}

Vec3f min3f(Vec3f left, Vec3f right)
{
    Vec3f result = {
        .x = MIN(left.x, right.x),
        .y = MIN(left.y, right.y),
        .z = MIN(left.z, right.z),
    };

    return result;
}

Vec3f pow3f(Vec3f left, float power)
{
    Vec3f result = {
        .x = powf(left.x, power),
        .y = powf(left.y, power),
        .z = powf(left.z, power),
    };

    return result;
}

Vec3f cross3f(Vec3f left, Vec3f right)
{
    Vec3f result = {
        .x = left.y * right.z - left.z * right.y,
        .y = left.z * right.x - left.x * right.z,
        .z = left.x * right.y - left.y * right.x,
    };

    return result;
}

Vec3f scale3f(Vec3f vec, float scale)
{
    Vec3f result = {
        .x = vec.x * scale,
        .y = vec.y * scale,
        .z = vec.z * scale,
    };

    return result;
}

Vec3f project3f(Vec3f onto, Vec3f projectee)
{
    Vec3f result = scale3f(onto, dot3f(onto, projectee));

    return result;
}