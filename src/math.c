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

Vec3f div3f(Vec3f left, Vec3f right)
{
    Vec3f result = {
        .x = left.x / right.x,
        .y = left.y / right.y,
        .z = left.z / right.z,
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

Vec3f max3f(Vec3f left, Vec3f right)
{
    Vec3f result = {
        .x = MAX(left.x, right.x),
        .y = MAX(left.y, right.y),
        .z = MAX(left.z, right.z),
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

Random_Serie create_random_serie()
{
    Random_Serie serie;
    for (int i = 0; i < 32; i++) {
        serie.bytes[i] = rand() % 256;
    }

    return serie;
}

void random_advance(Random_Serie *serie)
{
    serie->xmm[0] = _mm_aesdec_si128(serie->xmm[0], serie->xmm[1]);
    serie->xmm[1] = _mm_aesdec_si128(serie->xmm[0], serie->xmm[1]);
}

Vec4f random4f_uniform(Random_Serie *serie, float min, float max)
{
    __m128i mantissa_mask  = _mm_set1_epi32((1 << 23) - 1);
    __m128i exponent_field = _mm_set1_epi32(127 << 23);

    float scale1 = max - min;
    float step1 = -max + 2 * min;

    __m128 scale = _mm_load1_ps(&scale1);
    __m128 step  = _mm_load1_ps(&step1);

    union {
        __m128i xmmi;
        __m128  xmmf;
    } xmm;

    random_advance(serie);

    xmm.xmmi = _mm_xor_si128(serie->xmm[0], serie->xmm[1]);
    xmm.xmmi = _mm_and_si128(xmm.xmmi, mantissa_mask);
    xmm.xmmi = _mm_or_si128(xmm.xmmi, exponent_field);

    xmm.xmmf = _mm_fmadd_ps(xmm.xmmf, scale, step);

    Vec4f result;
    _mm_storeu_ps(result.coord, xmm.xmmf);

    return result;
}

Vec3f random3f_uniform(Random_Serie *serie, float min, float max)
{
    Vec4f vec4f = random4f_uniform(serie, min, max);
    Vec3f result = {
        .x = vec4f.x,
        .y = vec4f.y,
        .z = vec4f.z,
    };

    return result;
}

Vec2f random2f_uniform(Random_Serie *serie, float min, float max)
{
    Vec4f vec4f = random4f_uniform(serie, min, max);
    Vec2f result = {
        .x = vec4f.x,
        .y = vec4f.y,
    };

    return result;
}

float random_uniform(Random_Serie *serie, float min, float max)
{
    Vec4f vec4f = random4f_uniform(serie, min, max);
    float result = vec4f.x;

    return result;
}

Vec2f random2f_gaussian(Random_Serie *serie, float mean, float sigma)
{
    Vec4f u = random4f_uniform(serie, 0, 1);

    Vec2f result = {
        .x = mean + sigma * cosf(2.0f *  M_PI * u.x) * sqrtf(-2 * logf(u.y)),
        .y = mean + sigma * sinf(2.0f *  M_PI * u.x) * sqrtf(-2 * logf(u.y)),
    };

    return result;
}

Vec4f random4f_gaussian(Random_Serie *serie, float mean, float sigma)
{
    Vec4f u = random4f_uniform(serie, 0, 1);

    float log1 = logf(u.y);
    float log2 = logf(u.w);

    Vec4f result = {
        .x = mean + sigma * cosf(2.0f *  M_PI * u.x) * sqrtf(-2 * log1),
        .y = mean + sigma * sinf(2.0f *  M_PI * u.x) * sqrtf(-2 * log1),
        .z = mean + sigma * cosf(2.0f *  M_PI * u.z) * sqrtf(-2 * log2),
        .w = mean + sigma * sinf(2.0f *  M_PI * u.z) * sqrtf(-2 * log2),
    };

    return result;
}

Vec3f random3f_sphere(Random_Serie *serie, float radius)
{
    Vec3f result;

    do
    {
        result = random3f_uniform(serie, -1, 1);
    } while (norm2(result) > 1);

    result = normalize(result);
    result = scale3f(result, radius);

    return result;
}

Vec3f sample_texture(Image3f *texture, Vec2f uv)
{
    int ui = uv.x * (texture->size.x - 0.5);
    int vi = uv.y * (texture->size.y - 0.5);

    Vec3f color = texture->pixels[vi * texture->size.x + ui];
    return color;
}

Vec3f interpolate_triangle3f(Vec3f a, Vec3f b, Vec3f c, Vec3f d)
{
    Vec3f result = {
        .x = a.x * d.x + b.x * d.y + c.x * d.z,
        .y = a.y * d.x + b.y * d.y + c.y * d.z,
        .z = a.z * d.x + b.z * d.y + c.z * d.z,
    };

    return result;
}

Vec2f interpolate_triangle2f(Vec2f a, Vec2f b, Vec2f c, Vec3f d)
{
    Vec2f result = {
        .x = a.x * d.x + b.x * d.y + c.x * d.z,
        .y = a.y * d.x + b.y * d.y + c.y * d.z,
    };

    return result;
}

Ray transform_ray_to_mesh_space(Ray original_ray, Transform transform)
{
    Vec3f invRotX = {transform.rotX.x, transform.rotY.x, transform.rotZ.x};
    Vec3f invRotY = {transform.rotX.y, transform.rotY.y, transform.rotZ.y};
    Vec3f invRotZ = {transform.rotX.z, transform.rotY.z, transform.rotZ.z};

    Ray result = {
        .origin = interpolate_triangle3f(invRotX, invRotY, invRotZ, sub3f(original_ray.origin, transform.position)),
        .direction = interpolate_triangle3f(invRotX, invRotY, invRotZ, original_ray.direction),
    };

    return result;
}