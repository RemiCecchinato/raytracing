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

Vec3f raytrace_ray(Job *job, Scene *scene, Ray ray)
{
    Camera *camera = &scene->camera;

    Vec3f color_scale = WHITE;
    Vec3f accumulated_color = BLACK;

    // Keeps tracks of whether we are inside a transparent object or not.
    // With this method it is impossible to have transparent objects intersecting.
    bool within_transparent_object = false;

    for (int ray_count = 0; ray_count < 64; ray_count++) {
        Ray_Result ray_result = find_ray_hit(scene, ray);

        if (!ray_result.hit) break;

        Material *material = scene->materials + ray_result.material_id;

        if (material->transparent) {
            Vec3f incident_normal = project3f(ray_result.surface_normal, ray.direction);
            Vec3f incident_tangent = sub3f(ray.direction, incident_normal);

            float index_ratio = within_transparent_object ? (material->n) : (1.0f / material->n);

            Vec3f transmited_tangent = scale3f(incident_tangent, index_ratio);

            if (norm2(transmited_tangent) > 1) {
                // Incomming angle is too big, simply reflect on next ray.

                ray.origin = ray_result.intersection_point;
                ray.direction = sub3f(incident_tangent, incident_normal);
            } else {
                // Transmit the ray.
                within_transparent_object = !within_transparent_object;

                Vec3f transmited_normal = scale3f(normalize(incident_normal), sqrtf(1 - norm2(transmited_tangent)));
                ray.origin = ray_result.intersection_point;
                ray.direction = add3f(transmited_normal, transmited_tangent);
            }
        } else { // Non transparent material.
            Light light = {
                .position = {-10, 20, 40},
                .color = WHITE,
                .intensity = 5e4f,
            };

            Ray ray_to_light = {
                .origin = ray_result.intersection_point,
                .direction = normalize(sub3f(light.position, ray_result.intersection_point)),
            };

            Ray_Result light_ray_result = find_ray_hit(scene, ray_to_light);

            Vec3f light_path = sub3f(ray_result.intersection_point, light.position);
            Vec3f light_direction = normalize(light_path);

            float dist2 = norm2(light_path);

            if (!light_ray_result.hit || light_ray_result.t_min * light_ray_result.t_min > dist2) {
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

            Vec2f r = random2f_uniform(&job->serie, 0, 1);
            Vec3f random_direction = {
                .x = cosf(2.0f * M_PI * r.x) * sqrtf(1 - r.y),
                .y = sinf(2.0f * M_PI * r.x) * sqrtf(1 - r.y),
                .z = sqrtf(r.y),
            };

            Vec3f dir_normal = ray_result.surface_normal; // scale3f(ray_result.surface_normal, dot3f(ray_result.surface_normal, ray.direction));
            Vec3f dir_tan0 = {};
            if (fabsf(dir_normal.x) >= fabsf(dir_normal.y) && fabsf(dir_normal.x) >= fabsf(dir_normal.z)) {
                dir_tan0.x =   dir_normal.y;
                dir_tan0.y = - dir_normal.x;
                dir_tan0.z = 0;
            } else if (fabsf(dir_normal.y) >= fabsf(dir_normal.x) && fabsf(dir_normal.y) >= fabsf(dir_normal.z)) {
                dir_tan0.x =   dir_normal.y;
                dir_tan0.y = - dir_normal.x;
                dir_tan0.z = 0;
            } else if (fabsf(dir_normal.z) >= fabsf(dir_normal.x) && fabsf(dir_normal.z) >= fabsf(dir_normal.y)) {
                dir_tan0.x = 0;
                dir_tan0.y = - dir_normal.z;
                dir_tan0.z = dir_normal.y;
            }
            Vec3f dir_tan1 = cross3f(dir_normal, dir_tan0);

            Vec3f new_dir = add3f(scale3f(dir_normal, random_direction.z),
                            add3f(scale3f(dir_tan0, - random_direction.y),
                                  scale3f(dir_tan1, - random_direction.x)));

            new_dir = normalize(new_dir);

            color_scale = scale3f(color_scale, dot3f(new_dir, dir_normal));

            ray.direction = new_dir; // add3f(ray.direction, scale3f(new_dir, -2.0f));
        }
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
    Vec3f camera_right = normalize(cross3f(camera_direction, camera->up));
    Vec3f camera_up = cross3f(camera_right, camera_direction);

    for (uint32_t y = y_start; y < y_end; y++) {
        for (uint32_t x = x_start; x < x_end; x++) {
            Vec3f total_color = {};
            float sigma = 0.5f;

            for (int i = 0; i < job->ray_per_pixel; i++) {
                // Vec2f random_direction = random2f_uniform(&job->serie, -0.5f, 0.5f);
                Vec2f random_direction = random2f_gaussian(&job->serie, 0, 0.25f);
                Vec3f direction_right = scale3f(camera_right, (float)x - (float)image->size.width / 2.0f + 0.5f + random_direction.x);
                Vec3f direction_up    = scale3f(camera_up, -((float)y - (float)image->size.height / 2.0f + 0.5f + random_direction.y));
                Vec3f direction_main  = scale3f(camera_direction, (float)image->size.height / (2.0f * tanf(camera->field_of_view / 2.0f)));
                Vec3f direction = normalize(add3f(add3f(direction_right, direction_up), direction_main));

                Ray ray = {
                    .origin = camera->position,
                    .direction = direction,
                };

                Vec3f ray_color = raytrace_ray(job, scene, ray);

                // float p_x = expf(- random_direction.x * random_direction.x / (2.0f * sigma * sigma));
                // float p_y = expf(- random_direction.y * random_direction.y / (2.0f * sigma * sigma));
                // ray_color = scale3f(ray_color, 1.0f / (p_x * p_y));

                total_color = add3f(total_color, ray_color);
            }

            float factor = 1.0f / (sigma * sigma * 2.0f * M_PI * job->ray_per_pixel);
            image->pixels[image->size.width * y + x] = scale3f(total_color, factor);
        }
    }
}

int main()
{
#if 0
    Random_Serie serie = create_random_serie();

    float sum = 0;
    int N = 1000000;
    for (int i = 0; i < N; i++) {
        #if 0
        float sigma = 0.5;
        float x = random2f_gaussian(&serie, 0, 0.5f).x;
        if ((x < -M_PI_2) || (x > M_PI_2)) continue;

        float p = 1.0f / (sigma * sqrtf(2.0f * M_PI)) * expf(- x * x / (2.0f * sigma * sigma));
        float f_x = powf(cosf(x), 10.0f);
        #else
        float sigma = 1;
        Vec4f r4 = random4f_gaussian(&serie, 0, sigma);
        float x = r4.x + r4.y + r4.z + r4.w;

        if ((r4.x < -M_PI_2) || (r4.x > M_PI_2)) continue;
        if ((r4.y < -M_PI_2) || (r4.y > M_PI_2)) continue;
        if ((r4.z < -M_PI_2) || (r4.z > M_PI_2)) continue;
        if ((r4.w < -M_PI_2) || (r4.w > M_PI_2)) continue;

        float p_x = 1.0f / (sigma * sqrtf(2.0f * M_PI)) * expf(- r4.x * r4.x / (2.0f * sigma * sigma));
        float p_y = 1.0f / (sigma * sqrtf(2.0f * M_PI)) * expf(- r4.y * r4.y / (2.0f * sigma * sigma));
        float p_z = 1.0f / (sigma * sqrtf(2.0f * M_PI)) * expf(- r4.z * r4.z / (2.0f * sigma * sigma));
        float p_w = 1.0f / (sigma * sqrtf(2.0f * M_PI)) * expf(- r4.w * r4.w / (2.0f * sigma * sigma));
        float p = p_x * p_y * p_z * p_w;

        float f_x = powf(cosf(x), 2.0f);
        #endif

        sum += f_x / p;
    }
    sum /= N;

    printf("%f\n", sum);

    return 0;
#endif

    int width = 512, height = 512;

    Image3f dest_image = create_blank_image(width, height);

    Material materials[] = {
        {   // The sphere in the center of the scene.
            .diffuse_color = scale3f(WHITE, 0.5f),
            .mirror_color = scale3f(WHITE, 0.1f),
        },
        {   // The background wall.
            .diffuse_color = GREEN,
            .mirror_color = BLACK,
        },
        {   // The ground.
            .diffuse_color = BLUE,
            .mirror_color = BLACK,
        },
        {   // Left wall.
            .diffuse_color = PINK,
            .mirror_color = BLACK,
        },
        {   // Right Wall.
            .diffuse_color = YELLOW,
            .mirror_color = BLACK,
        },
        {   // Perfect mirror sphere.
            .diffuse_color = BLACK,
            .mirror_color = WHITE,
        },
        {   // Ceilling.
            .diffuse_color = RED,
            .mirror_color = BLACK,
        },
        {   // Transparent sphere.
            .diffuse_color = BLACK,
            .mirror_color = BLACK,

            .transparent = true,
            .n = 1.5f,
        },
    };

    Sphere spheres[] = {
        {   // The sphere in the center of the scene.
            .center = {0.0f, 0.0f, 0.0f},
            .radius = 10.0f,
            .material_id = 0,
        },
        {   // Left-most sphere, verry mirrory.
            .center = {-8.0f, 10.0f, 15.0f},
            .radius = 5.0f,
            .material_id = 5,
        },
        {   // Right-most sphere, transparent.
            .center = {-4.0f, 5.0f, 25.0f},
            .radius = 4.0f,
            .material_id = 7,
        },
#if 0
#endif
    };

    Plane planes[] = {
        {
            .normal = {0.0f, 0.0f, 1.0f},
            .d = -60.0f,
            .material_id = 1,
        },
        {
            .normal = {0.0f, 1.0f, 0.0f},
            .d = -10.0f,
            .material_id = 2,
        },
        {
            .normal = {1.0f, 0.0f, 0.0f},
            .d = -60.0f,
            .material_id = 3,
        },
        {
            .normal = {-1.0f, 0.0f, 0.0f},
            .d = -60.0f,
            .material_id = 4,
        },
        {
            .normal = {0.0f, -1.0f, 0.0f},
            .d = -90.0f,
            .material_id = 6,
        },
    };

    Scene scene = {
        .camera = {
            .position = {0.0f, 10.0f, 55.0f},
            .lookAt = {0.0f, 0.0f, 0.0f},
            .up = {0.0f, 1.0f, 0.0f},
            .field_of_view = 60.0f / 180.0f * M_PI,
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
                jobs[4 * i +j].ray_per_pixel = 256;
                jobs[4 * i +j].serie = create_random_serie();

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
            .ray_per_pixel = 1,
            .serie = create_random_serie(),
        };

        raytrace_region(&job);
    }

    save_image("image.png", &dest_image);

    return 0;
}