#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <pthread.h>

#if !defined(USE_THREADS)
    #define USE_THREADS 0
#endif

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#if 0
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#endif

// Mes headers
#include "main.h"

// Mes sources
#include "math.c"

Image3f create_blank_image(uint32_t width, uint32_t height)
{
    Image3f image = {};

    uint32_t pixel_count = width * height;
    uint32_t bpp = sizeof(*image.pixels);

    // @Leak
    Vec3f *memory = (Vec3f*)((uint64_t)(calloc(pixel_count * bpp + 64, 1) + 63) & ~64LL);

    if (!memory) return image;

    image.size.width  = width;
    image.size.height = height;
    image.pixels = memory;

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
            Vec3f pixel_color = min3f(image->pixels[j * image->size.width + i], (Vec3f)WHITE);
            pixel_color = pow3f(pixel_color, 1.0f / 2.2f);
            pixel_color = scale3f(pixel_color, 255.0f);

            *pixel++ = (uint8_t)pixel_color.r;
            *pixel++ = (uint8_t)pixel_color.g;
            *pixel++ = (uint8_t)pixel_color.b;
        }
    }

    stbi_write_png(filename, image->size.width, image->size.height, 3, pixels, 0);

    free(pixels);

    return true;
}

Ray_Result find_ray_hit(Scene *scene, Ray ray)
{
    Ray_Result result = {
        .hit = false,

        .intersection_point = {},
        .surface_normal = {},
        .t_min = FLT_MAX,
        
        .enter_shape = false,
        
        .material_id = -1,
    };

    // On teste toutes les sphères
    for (int i = 0; i < scene->sphere_count; i++) {
        Sphere *sphere = scene->spheres + i;

        Vec3f OC = sub3f(ray.origin, sphere->center);

        // On calcule un determinant réduit qui évite les multiplications par 4 puis les divions par 4
        float b = dot3f(ray.direction, OC);
        float c = norm2(OC) - sphere->radius * sphere->radius;

        float det = b * b - c;

        if (det < 0.0f) continue; // Pas d'intersection

        float sqrt_det = sqrtf(det);

        // Calcul du t
        float t = -b - sqrt_det;
        
        // si t < 0, alors on est à l'intérieur de la sphère (bizard ...)
        // donc on prends le second point, plus loin
        if (t < 1e-3) {
            t = -b + sqrt_det;
        } else {
            result.enter_shape = true;
        }

        // On est pas la sphère la plus proche
        if (t > result.t_min) continue;
        if (t < 1e-3) continue;

        result.hit = true;
        result.t_min = t;
        result.material_id = sphere->material_id;

        result.intersection_point = add3f(ray.origin, scale3f(ray.direction, t));
        result.surface_normal = normalize(sub3f(result.intersection_point, sphere->center));
    }

    // On teste tous les plans
    for (int i = 0; i < scene->plane_count; i++) {
        Plane *plane = scene->planes + i;

        float t = (plane->d - dot3f(ray.origin, plane->normal)) / dot3f(ray.direction, plane->normal);

        if (t < 1e-3) continue;

        if (t > result.t_min) continue;

        result.hit = true;
        result.t_min = t;
        result.material_id = plane->material_id;

        result.intersection_point = add3f(ray.origin, scale3f(ray.direction, t));
        result.surface_normal = plane->normal;
    }

    return result;
}

Vec3f raytrace_ray(Scene *scene, Ray ray)
{
    Camera *camera = &scene->camera;
    
    Vec3f color_scale = WHITE;
    Vec3f accumulated_color = BLACK;

    for (int ray_count = 0; ray_count < 64; ray_count++) {
        Ray_Result ray_result = find_ray_hit(scene, ray);

        if (!ray_result.hit) break;

        Material *material = scene->materials + ray_result.material_id;

        Light light = {
            .position = {-5000.0, -10000.0f, 1000.0f},
            .color = WHITE,
            .intensity = 6e8f,
        };

        Ray ray_to_light = {
            .origin = ray_result.intersection_point,
            .direction = normalize(sub3f(light.position, ray_result.intersection_point)),
        };

        Ray_Result light_ray_result = find_ray_hit(scene, ray_to_light);

        if (!light_ray_result.hit) {
            Vec3f light_path = sub3f(ray_result.intersection_point, light.position);
            Vec3f light_direction = normalize(light_path);

            float dist2 = norm2(light_path);
            float power = light.intensity / (4.0f * M_PI * dist2);

            // Vec3f max_light_direction = add3f(light_direction, scale3f(ray_result.surface_normal, - 2 * dot3f(ray_result.surface_normal, light_direction)));
            // float scale = - dot3f(max_light_direction, direction);
            float scale = - dot3f(light_direction, ray_result.surface_normal);
            if (scale < 0) scale = 0;

            Vec3f final_light_color = scale3f(light.color, power * scale);

            Vec3f color = mul3f(final_light_color, material->diffuse_color);
            
            accumulated_color = add3f(accumulated_color, mul3f(color_scale, color));
        }

        color_scale = mul3f(color_scale, material->mirror_color);

        // Si on ajoute plus assez de lumière après chaque rebond, on arrête
        if (norm2(color_scale) < 1e-3) break;

        ray.origin = ray_result.intersection_point;

        Vec3f dir_normal = scale3f(ray_result.surface_normal, dot3f(ray_result.surface_normal, ray.direction));
        ray.direction = add3f(ray.direction, scale3f(dir_normal, -2.0f));
    }

    return accumulated_color;
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
    Vec3f camera_right = cross3f(camera_direction, camera_up);

    float fov_h = camera->field_of_view;
    float fov_v = (float)image->size.height / (float)image->size.width * fov_h;

    Vec3f pixel_right = scale3f(camera_right, sinf(fov_h * 0.5f) / (float)image->size.width * 2);
    Vec3f pixel_up    = scale3f(camera_up, sinf(fov_v * 0.5f) / (float)image->size.height * 2);

    Vec2f pixel_center = {
        .x = (float)image->size.width / 2,
        .y = (float)image->size.height / 2,
    };

    for (uint32_t y = y_start; y < y_end; y++) {
        for (uint32_t x = x_start; x < x_end; x++) {
            Vec2i pixel_coord = {.x = x, .y = y};

            Vec2f pixel_coord_f = {.x = (float)x, .y = (float)y};

            Vec3f direction = camera_direction;
            direction = add3f(direction, scale3f(pixel_right, pixel_coord_f.x - pixel_center.x));
            direction = add3f(direction, scale3f(pixel_up,  -(pixel_coord_f.y - pixel_center.y)));
            direction = normalize(direction);

            Ray ray = {
                .origin = camera->position,
                .direction = direction,
            };

            Vec3f ray_color = raytrace_ray(scene, ray);

            image->pixels[image->size.width * pixel_coord.y + pixel_coord.x] = ray_color;
        }
    }
}

int main()
{
    int width = 512, height = 512;

    Image3f dest_image = create_blank_image(width, height);

    Material materials[] = {
        {
            .diffuse_color = RED,
            .mirror_color = scale3f(WHITE, 0.1),
        },
        {
            .diffuse_color = WHITE,
            .mirror_color = BLACK,
        },
        {
            .diffuse_color = GREEN,
            .mirror_color = scale3f(WHITE, 0.5),
        },
    };

    Sphere spheres[] = {
        {
            .center = {0.0f, 0.0f, 0.0f},
            .radius = 1.0f,
            .material_id = 0,
        },
#if 1
        {
            .center = {2.0f, 0.0f, 0.5f},
            // .center = {0.0f, -1.0f, 0.0f},
            .radius = 0.7f,
            .material_id = 2,
        }
#endif
    };

    Plane planes[] = {
        {
            .normal = {0.0f, -1.0f, 0.0f},
            .d = -10000.0f,
            .material_id = 1,
        }
    };

    Scene scene = {
        .camera = {
            .position = {0.0f, -3.0f, 0.0f},
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

    if(USE_THREADS) {
        pthread_t threads[16] = {};
        Job jobs[16] = {};

        uint32_t block_width = width / 4;
        uint32_t block_height = height / 4;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                jobs[4 * i + j].scene = &scene;
                jobs[4 * i + j].image = &dest_image;
                jobs[4 * i + j].region = (Rect2i){
                    .point = {i * block_width, j * block_height},
                    .size = {block_width, block_height}
                };

                int code = pthread_create(&threads[4 * i + j], NULL, (void*)raytrace_region, &jobs[4 * i + j]);
                if (0 != code) {
                    printf("Error creating thread.\n");
                    exit(-1);
                }
            }
        }

        for (int i = 0; i < 16; i++) {
            pthread_join(threads[i], NULL);
        }
    } else {
        Job job = {
            .scene = &scene,
            .image = &dest_image,
            .region = {
                .point = {0, 0},
                .size = {width, height},
            },
        };

        raytrace_region(&job);
    }

    save_image("image.png", &dest_image);

    return 0;
}