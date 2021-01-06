#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Mes headers
#include "main.h"

// Mes sources
#include "math.c"

Image3f create_blank_image(uint32_t width, uint32_t height)
{
    Image3f image = {};

    uint32_t pixel_count = width * height;
    uint32_t bpp = sizeof(*image.r) * 3;

    float *memory = calloc(pixel_count, bpp);
    if (!memory) return image;

    image.size.width  = width;
    image.size.height = height;
    image.r = memory;
    image.g = memory + pixel_count;
    image.b = memory + 2 * pixel_count;

    return image;
};

bool save_image(char *filename, Image3f *image)
{
    if (!filename) return false;
    if (!image) return false;

    uint8_t *pixels = malloc(image->size.width * image->size.height * 3);

    uint8_t *pixel = pixels;
    for (int j = 0; j < image->size.height; j++) {
        for (int i = 0; i < image->size.width; i++) {
            *pixel++ = (uint8_t)(image->r[j * image->size.width + i] * 255.0f);
            *pixel++ = (uint8_t)(image->g[j * image->size.width + i] * 255.0f);
            *pixel++ = (uint8_t)(image->b[j * image->size.width + i] * 255.0f);
        }
    }

    stbi_write_png(filename, image->size.width, image->size.height, 3, pixels, 0);

    free(pixels);

    return true;
}

void raytrace_ray(Job *job, Vec2i pixel_coord, Vec3f direction)
{
    Scene *scene = job->scene;
    Camera *camera = &scene->camera;

    float t_min = FLT_MAX;
    int32_t material_id = -1;
    Vec3f intersection_point = {};
    Vec3f surface_normal = {};

    // On teste toutes les sphères
    for (int i = 0; i < scene->sphere_count; i++) {
        Sphere *sphere = scene->spheres + i;

        Vec3f OC = sub3f(camera->position, sphere->center);

        // On calcule un determinant réduit qui évite les multiplications par 4 puis les divions par 4
        float b = dot3f(direction, OC);
        float c = norm2(OC) - sphere->radius * sphere->radius;

        float det = b * b - c;

        if (det < 0.0f) continue; // Pas d'intersection

        float sqrt_det = sqrtf(det);

        // Calcul du t
        float t = -b - sqrt_det;
        
        // si t < 0, alors on est à l'intérieur de la sphère (bizard ...)
        // donc on prends le second point, plus loin
        if (t < 0.0f) t = -b + sqrt_det;

        // On est pas la sphère la plus proche
        if (t > t_min) continue;

        t_min = t;
        material_id = sphere->material_id;

        intersection_point = add3f(camera->position, scale3f(direction, t));
        surface_normal = normalize(sub3f(intersection_point, sphere->center));
    }

    // On teste tous les plans
    for (int i = 0; i < scene->plane_count; i++) {
        Plane *plane = scene->planes + i;

        float t = (plane->d - dot3f(camera->position, plane->normal)) / dot3f(direction, plane->normal);

        if (t < 0) continue;

        if (t > t_min) continue;

        t_min = t;
        material_id = plane->material_id;

        intersection_point = add3f(camera->position, scale3f(direction, t));
        surface_normal = plane->normal;
    }

    if (material_id != -1) {
        Image3f *image = job->image;
        Material *material = scene->materials + material_id;

        image->r[image->size.width * pixel_coord.y + pixel_coord.x] = material->color.r;
        image->g[image->size.width * pixel_coord.y + pixel_coord.x] = material->color.g;
        image->b[image->size.width * pixel_coord.y + pixel_coord.x] = material->color.b;
    }
}

void raytrace_region(Job *job)
{
    Scene *scene = job->scene;
    Camera *camera = &scene->camera;
    Image3f *image = job->image;

    uint32_t x_start = job->region.point.x;
    uint32_t y_start = job->region.point.y;
    uint32_t x_end = job->region.point.x + job->region.size.x;
    uint32_t y_end = job->region.point.y + job->region.size.y;

    Vec3f camera_direction = normalize(sub3f(camera->lookAt, camera->position));
    Vec3f camera_up = normalize(camera->up);
    Vec3f camera_right = cross3f(camera_up, camera_direction);

    float fov_h = camera->field_of_view / (float)image->size.width;
    float fov_v = (float)image->size.height / (float)image->size.width * fov_h;

    Vec2f pixel_center = {
        .x = (float)image->size.width / 2,
        .y = (float)image->size.height / 2,
    };

    for (uint32_t y = y_start; y < y_end; y++) {
        for (uint32_t x = x_start; x < x_end; x++) {
            Vec2i pixel_coord = {.x = x, .y = y};

            Vec2f pixel_coord_f = {.x = (float)x, .y = (float)y};
            float angle_h = fov_h * (pixel_coord_f.x - pixel_center.x);
            float angle_v = fov_v * (pixel_coord_f.y - pixel_center.y);

            float scale_h = sinf(angle_h / 180.0f * M_PI);
            float scale_v = sinf(angle_v / 180.0f * M_PI);

            Vec3f direction = camera_direction;
            direction = add3f(direction, scale3f(camera_right, scale_h));
            direction = add3f(direction, scale3f(camera_up, scale_v));
            direction = normalize(direction);

            raytrace_ray(job, pixel_coord, direction);
        }
    }
}

int main()
{
    int width = 512, height = 512;

    Image3f dest_image = create_blank_image(width, height);

    Material materials[] = {
        {
            .color = {1.0f, 0.0f, 0.0f},
        },
        {
            .color = {1.0f, 1.0f, 1.0f},
        }
    };

    Sphere spheres[] = {
        {
            .center = {0.0f, 0.0f, 0.0f},
            .radius = 1.0f,
            .material_id = 0,
        },
    };

    Plane planes[] = {
        {
            .normal = {1.0f, 0.0f, 0.0f},
            .d = 1000.0f,
            .material_id = 1,
        }
    };

    Scene scene = {
        .camera = {
            .position = {-5.0f, 0.0f, 0.0f},
            .lookAt = {0.0f, 0.0f, 0.0f},
            .up = {0.0f, 0.0f, 1.0f},
            .field_of_view = 90,
        },

        .material_count = ARRAY_SIZE(materials),
        .materials = materials,

        .sphere_count = ARRAY_SIZE(spheres),
        .spheres = spheres,

        .plane_count = ARRAY_SIZE(planes),
        .planes = planes,
    };

    Job job = {
        .scene = &scene,
        .image = &dest_image,
        .region = {
            .point = {0.0f, 0.0f},
            .size = {width, height},
        },
    };

    raytrace_region(&job);

    save_image("image.png", &dest_image);

    return 0;
}