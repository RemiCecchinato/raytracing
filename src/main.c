#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <stdbool.h>
#include <pthread.h>

#if !defined(USE_THREADS)
    #define USE_THREADS 0
#endif

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
    uint32_t bpp = sizeof(*image.pixels);

    // @Leak
    Vec3f *memory = (Vec3f*)((uint64_t)(calloc(pixel_count * bpp + 64, 1) + 63) & ~64LL);

    if (!memory) return image;

    image.size.width  = width;
    image.size.height = height;
    image.pixels = memory;

    return image;
};

Image3f load_texture(char *filename)
{
    int channels = 0, width = 0, height = 0;
    uint8_t *pixels = stbi_load(filename, &width, &height, &channels, 0);

    Image3f texture = create_blank_image(width, height);

    uint32_t pixel_count = texture.size.x * texture.size.y;
    for (uint32_t j = 0; j < texture.size.y; j++) {
        for (uint32_t i = 0; i < texture.size.x; i++) {
            texture.pixels[texture.size.x * (texture.size.y - 1 - j) + i].r = powf(pixels[3 * (texture.size.x * j + i) + 0] / 255.0f, 2.2f);
            texture.pixels[texture.size.x * (texture.size.y - 1 - j) + i].g = powf(pixels[3 * (texture.size.x * j + i) + 1] / 255.0f, 2.2f);
            texture.pixels[texture.size.x * (texture.size.y - 1 - j) + i].b = powf(pixels[3 * (texture.size.x * j + i) + 2] / 255.0f, 2.2f);
        }
    }

    free(pixels);
    return texture;
}

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

void split_mesh_node(Mesh_Node *node, Mesh *mesh, int depth)
{
    mesh->tree_depth = MAX(mesh->tree_depth, depth);

    node->bbox.min.x =  FLT_MAX;
    node->bbox.min.y =  FLT_MAX;
    node->bbox.min.z =  FLT_MAX;
    node->bbox.max.x = -FLT_MAX;
    node->bbox.max.y = -FLT_MAX;
    node->bbox.max.z = -FLT_MAX;

    for (uint32_t i = node->first_triangle_id; i <= node->last_triangle_id; i++) {
        Triangle *triangle = mesh->faces + i;

        for (int j = 0; j < 3; j++) {
            uint32_t vertex_id = triangle->ids[j].vertexId;
            
            node->bbox.min.x = MIN(node->bbox.min.x, mesh->points[vertex_id].x);
            node->bbox.min.y = MIN(node->bbox.min.y, mesh->points[vertex_id].y);
            node->bbox.min.z = MIN(node->bbox.min.z, mesh->points[vertex_id].z);

            node->bbox.max.x = MAX(node->bbox.max.x, mesh->points[vertex_id].x);
            node->bbox.max.y = MAX(node->bbox.max.y, mesh->points[vertex_id].y);
            node->bbox.max.z = MAX(node->bbox.max.z, mesh->points[vertex_id].z);
        }
    }

    // Si il y a moins que 5 triangles, on ne divise pas le noeud.
    if (node->last_triangle_id - node->first_triangle_id + 1 <= 16) {
        return;
    }

    int maxdim = 0;
    if (node->bbox.max.coord[1] - node->bbox.min.coord[1] > node->bbox.max.coord[maxdim] - node->bbox.min.coord[maxdim]) {
        maxdim = 1;
    }
    if (node->bbox.max.coord[2] - node->bbox.min.coord[2] > node->bbox.max.coord[maxdim] - node->bbox.min.coord[maxdim]) {
        maxdim = 2;
    }

    float middle = (node->bbox.max.coord[maxdim] + node->bbox.min.coord[maxdim]) / 2;

    int pivot = node->first_triangle_id;
    for (;;) {
        for (uint32_t i = node->first_triangle_id; i <= node->last_triangle_id; i++) {
            Triangle *triangle = mesh->faces + i;

            float center = 0;
            for (int j = 0; j < 3; j++) {
                center += mesh->points[triangle->ids[j].vertexId].coord[maxdim] / 3;
            }

            if (center < middle) {
                Triangle tmp = mesh->faces[pivot];
                mesh->faces[pivot] = mesh->faces[i];
                mesh->faces[i] = tmp;

                pivot++;
            }
        }

        if (pivot == node->first_triangle_id) {
            middle = (middle + node->bbox.max.coord[maxdim]) / 2;
        } else if (pivot == node->last_triangle_id + 1) {
            middle = (node->bbox.min.coord[maxdim] + middle) / 2;
            pivot = node->first_triangle_id;
        } else {
            break;
        }
    }

    assert(pivot != node->first_triangle_id);
    assert(pivot != node->last_triangle_id+1);

    node->childs = calloc(1, 2 * sizeof(*node->childs));

    node->childs[0].first_triangle_id = node->first_triangle_id;
    node->childs[0].last_triangle_id  = pivot - 1;
    node->childs[1].first_triangle_id = pivot;
    node->childs[1].last_triangle_id  = node->last_triangle_id;

    split_mesh_node(node->childs + 0, mesh, depth + 1);
    split_mesh_node(node->childs + 1, mesh, depth + 1);
}

Mesh load_mesh(char *filename)
{
    printf("Chargement de %s...\n", filename);
    fastObjMesh* fastMesh = fast_obj_read(filename);
    assert(fastMesh);

    Mesh mesh = {};
    mesh.point_count    = fastMesh->position_count;
    mesh.points         = (Vec3f*)fastMesh->positions;
    mesh.texcoord_count = fastMesh->texcoord_count;
    mesh.texcoords      = (Vec2f*)fastMesh->texcoords;
    mesh.normal_count   = fastMesh->normal_count;
    mesh.normals        = (Vec3f*)fastMesh->normals;

    for (uint32_t i = 0; i < mesh.point_count; i++) {
        float t = mesh.points[i].y;
        mesh.points[i].y = mesh.points[i].z;
        mesh.points[i].z = t;
    }

    for (uint32_t i = 0; i < mesh.point_count; i++) {
        float t = mesh.points[i].x;
        mesh.points[i].x = mesh.points[i].z;
        mesh.points[i].z = t;
    }

    for (uint32_t i = 0; i < mesh.normal_count; i++) {
        float t = mesh.normals[i].y;
        mesh.normals[i].y = mesh.normals[i].z;
        mesh.normals[i].z = t;
    }

    for (uint32_t i = 0; i < mesh.normal_count; i++) {
        float t = mesh.normals[i].x;
        mesh.normals[i].x = mesh.normals[i].z;
        mesh.normals[i].z = t;
    }

    mesh.face_count = fastMesh->face_count * 2;
    mesh.faces      = malloc(sizeof(*mesh.faces) * mesh.face_count);

    assert(4 * fastMesh->face_count == array_size(fastMesh->indices));

    for (uint32_t i = 0; i < fastMesh->face_count; i++) {
        assert(fastMesh->face_vertices[i] == 4);

        mesh.faces[2 * i + 0].ids[0] = (Vertex){
            .vertexId   = fastMesh->indices[4 * i + 0].p,
            .texcoordId = fastMesh->indices[4 * i + 0].t,
            .normalId   = fastMesh->indices[4 * i + 0].n,
        };
        mesh.faces[2 * i + 0].ids[1] = (Vertex){
            .vertexId   = fastMesh->indices[4 * i + 1].p,
            .texcoordId = fastMesh->indices[4 * i + 1].t,
            .normalId   = fastMesh->indices[4 * i + 1].n,
        };
        mesh.faces[2 * i + 0].ids[2] = (Vertex){
            .vertexId   = fastMesh->indices[4 * i + 2].p,
            .texcoordId = fastMesh->indices[4 * i + 2].t,
            .normalId   = fastMesh->indices[4 * i + 2].n,
        };

        mesh.faces[2 * i + 1].ids[0] = (Vertex){
            .vertexId   = fastMesh->indices[4 * i + 0].p,
            .texcoordId = fastMesh->indices[4 * i + 0].t,
            .normalId   = fastMesh->indices[4 * i + 0].n,
        };
        mesh.faces[2 * i + 1].ids[1] = (Vertex){
            .vertexId   = fastMesh->indices[4 * i + 2].p,
            .texcoordId = fastMesh->indices[4 * i + 2].t,
            .normalId   = fastMesh->indices[4 * i + 2].n,
        };
        mesh.faces[2 * i + 1].ids[2] = (Vertex){
            .vertexId   = fastMesh->indices[4 * i + 3].p,
            .texcoordId = fastMesh->indices[4 * i + 3].t,
            .normalId   = fastMesh->indices[4 * i + 3].n,
        };
    }

    mesh.root.first_triangle_id = 0;
    mesh.root.last_triangle_id = mesh.face_count - 1;
    split_mesh_node(&mesh.root, &mesh, 1);

    // fast_obj_destroy(fastMesh);
    return mesh;
}

Ray_Result find_ray_hit(Scene *scene, Ray ray, float max_distance)
{
    Ray_Result result = {
        .hit = false,

        .intersection_point = {},
        .surface_normal = {},
        .t_min = max_distance,
        
        .enter_shape = false,
        
        .material_id = -1,
    };

    // On teste toutes les sources de lumière.
    // @CutNPaste des spheres
    for (int i = 0; i < scene->light_count; i++) {
        Light *light = scene->lights + i;

        Vec3f OC = sub3f(ray.origin, light->center);

        // On calcule un determinant réduit qui évite les multiplications par 4 puis les divions par 4
        float b = dot3f(ray.direction, OC);
        float c = norm2(OC) - light->radius * light->radius;

        float det = b * b - c;

        if (det < 0.0f) continue; // Pas d'intersection.

        float sqrt_det = sqrtf(det);

        // Calcul du t
        float t = -b - sqrt_det;
        
        // si t < 0, alors on est à l'intérieur de la lumière (bizard ...)
        // donc on prends le second point, plus loin.
        if (t < 1e-3) {
            t = -b + sqrt_det;
        } else {
            result.enter_shape = true;
        }

        // On est pas la lumière la plus proche.
        if (t > result.t_min) continue;
        if (t < 1e-3) continue;

        result.hit = true;
        result.object_type = LIGHT;
        result.object_id = i;
        result.t_min = t;
    }

    // On teste toutes les sphères.
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
        
        // si t < 0, alors on est à l'intérieur de la sphère
        // donc on prends le second point, plus loin.
        if (t < 1e-3) {
            t = -b + sqrt_det;
        } else {
            result.enter_shape = true;
        }

        // On est pas la sphère la plus proche
        if (t > result.t_min) continue;
        if (t < 1e-3) continue;

        result.hit = true;
        result.object_type = SPHERE;
        result.object_id = i;
        result.t_min = t;
        result.material_id = sphere->material_id;
    }

    // On teste tous les plans.
    for (int i = 0; i < scene->plane_count; i++) {
        Plane *plane = scene->planes + i;

        float t = (plane->d - dot3f(ray.origin, plane->normal)) / dot3f(ray.direction, plane->normal);

        if (t < 1e-3) continue;

        if (t > result.t_min) continue;

        result.hit = true;
        result.object_type = PLANE;
        result.object_id = i;
        result.t_min = t;
        result.material_id = plane->material_id;
    }

    // On teste tous les modèles.
    for (int i = 0; i < scene->mesh_count; i++) {
        Scene_Mesh *scene_mesh = scene->meshes + i;
        Mesh *mesh = scene_mesh->mesh;

        Ray mesh_ray = transform_ray_to_mesh_space(ray, scene_mesh->transform);

        {
            Vec3f t0 = div3f(sub3f(mesh->root.bbox.min, mesh_ray.origin), mesh_ray.direction);
            Vec3f t1 = div3f(sub3f(mesh->root.bbox.max, mesh_ray.origin), mesh_ray.direction);

            Vec3f tmin = min3f(t0, t1);
            Vec3f tmax = max3f(t0, t1);

            float t_min = MAX(MAX(tmin.x, tmin.y), tmin.z);
            float t_max = MIN(MIN(tmax.x, tmax.y), tmax.z);

            if (t_max < t_min) continue;

            if (t_min > result.t_min) continue;
        }

        int stack_head = 0;
        Mesh_Node **stack = alloca(sizeof(*stack) * (mesh->tree_depth + 1));
        stack[stack_head] = &mesh->root;

        while (stack_head >= 0) {
            Mesh_Node *node = stack[stack_head--];

            float t_min0, t_min1;
            bool hit0 = false, hit1 = false;
            if (node->childs) {
                {
                    Vec3f t0 = div3f(sub3f(node->childs[0].bbox.min, mesh_ray.origin), mesh_ray.direction);
                    Vec3f t1 = div3f(sub3f(node->childs[0].bbox.max, mesh_ray.origin), mesh_ray.direction);

                    Vec3f tmin = min3f(t0, t1);
                    Vec3f tmax = max3f(t0, t1);

                    float t_min = MAX(MAX(tmin.x, tmin.y), tmin.z);
                    float t_max = MIN(MIN(tmax.x, tmax.y), tmax.z);

                    if ((t_min < t_max) && (t_min < result.t_min)) {
                        t_min0 = t_min;
                        hit0 = true;
                    }
                }

                {
                    Vec3f t0 = div3f(sub3f(node->childs[1].bbox.min, mesh_ray.origin), mesh_ray.direction);
                    Vec3f t1 = div3f(sub3f(node->childs[1].bbox.max, mesh_ray.origin), mesh_ray.direction);

                    Vec3f tmin = min3f(t0, t1);
                    Vec3f tmax = max3f(t0, t1);

                    float t_min = MAX(MAX(tmin.x, tmin.y), tmin.z);
                    float t_max = MIN(MIN(tmax.x, tmax.y), tmax.z);

                    if ((t_min < t_max) && (t_min < result.t_min)) {
                        t_min1 = t_min;
                        hit1 = true;
                    }
                }

                if (hit0 && hit1) {
                    if (t_min0 > t_min1) {
                        stack[++stack_head] = node->childs + 0;
                        stack[++stack_head] = node->childs + 1; // LIFO
                    } else {
                        stack[++stack_head] = node->childs + 1;
                        stack[++stack_head] = node->childs + 0; // LIFO
                    }
                } else if (hit0) {
                    stack[++stack_head] = node->childs + 0;
                } else if (hit1) {
                    stack[++stack_head] = node->childs + 1;
                }
            } else {
                for (int j = node->first_triangle_id; j <= node->last_triangle_id; j++) {
                    Vec3f a = mesh->points[mesh->faces[j].ids[0].vertexId];
                    Vec3f b = mesh->points[mesh->faces[j].ids[1].vertexId];
                    Vec3f c = mesh->points[mesh->faces[j].ids[2].vertexId];

                    Vec3f e1 = sub3f(b, a);
                    Vec3f e2 = sub3f(c, a);

                    Vec3f o = mesh_ray.origin;
                    Vec3f u = mesh_ray.direction;

                    Vec3f oau = cross3f(sub3f(o, a), u);
                    Vec3f n = cross3f(e1, e2);

                    float t     = - dot3f(sub3f(o, a), n) / dot3f(u, n);
                    if (t < 1e-3) continue;
                    if (t > result.t_min) continue;

                    float beta  = - dot3f(e2, oau) / dot3f(u, n);
                    float gamma =   dot3f(e1, oau) / dot3f(u, n);

                    if ((beta < 0)  || (beta > 1))  continue;
                    if ((gamma < 0) || (gamma > 1)) continue;
                    if (beta + gamma > 1) continue;

                    result.hit = true;
                    result.t_min = t;
                    result.object_type = TRIANGLE;
                    result.object_id = j;
                    result.mesh_id = i;
                    result.triangle_beta  = beta;
                    result.triangle_gamma = gamma;
                    result.material_id = scene_mesh->material_id;
                }
            }
        }
    }

    if (result.hit) {
        switch (result.object_type) {
            case SPHERE: {
                Sphere *sphere = scene->spheres + result.object_id;

                result.intersection_point = add3f(ray.origin, scale3f(ray.direction, result.t_min));
                result.surface_normal = normalize(sub3f(result.intersection_point, sphere->center));
            } break;
            case PLANE: {
                Plane *plane = scene->planes + result.object_id;

                result.intersection_point = add3f(ray.origin, scale3f(ray.direction, result.t_min));
                result.surface_normal = plane->normal;
            } break;
            case TRIANGLE: {
                Scene_Mesh *scene_mesh = scene->meshes + result.mesh_id;
                Mesh *mesh = scene_mesh->mesh;
                Transform *transform = &scene_mesh->transform;
                Triangle *triangle = mesh->faces + result.object_id;

                result.triangle_alpha = 1 - result.triangle_beta - result.triangle_gamma;

                Vec3f alpha_beta_gamma = {
                    result.triangle_alpha,
                    result.triangle_beta,
                    result.triangle_gamma,
                };

                Vec3f normal0 = mesh->normals[triangle->ids[0].normalId];
                Vec3f normal1 = mesh->normals[triangle->ids[1].normalId];
                Vec3f normal2 = mesh->normals[triangle->ids[2].normalId];
                Vec3f normal = interpolate_triangle3f(normal0, normal1, normal2, alpha_beta_gamma);

                Vec2f uv0 = mesh->texcoords[triangle->ids[0].texcoordId];
                Vec2f uv1 = mesh->texcoords[triangle->ids[1].texcoordId];
                Vec2f uv2 = mesh->texcoords[triangle->ids[2].texcoordId];
                result.triangle_uv = interpolate_triangle2f(uv0, uv1, uv2, alpha_beta_gamma);

                result.intersection_point = add3f(ray.origin, scale3f(ray.direction, result.t_min));
                result.surface_normal = interpolate_triangle3f(transform->rotX, transform->rotY, transform->rotZ, normal);
            } break;
            case LIGHT: {
                // Rien à calculer.
            } break;
            default: {
                // On a rien à faire là.
            }
        }
    }

    return result;
}

typedef struct Light_Contribution
{
    Vec3f color;
    float light_hit;
} Light_Contribution;

Light_Contribution compute_light_contribution(Job *job, Scene *scene, Ray previous_ray, Ray_Result previous_ray_result)
{
    Light_Contribution contribution = {};

    for (int32_t i = 0; i < scene->light_count; i++) {
        Light *light = scene->lights + i;

        Vec2f r = random2f_uniform(&job->serie, 0, 1);
        Vec3f point_on_hemi_sphere = {
            .x = cosf(2 * M_PI * r.x) * sqrtf(2 * r.y - r.y * r.y),
            .y = sinf(2 * M_PI * r.x) * sqrtf(2 * r.y - r.y * r.y),
            .z = 1 - r.y,
        };

        Vec3f dir_normal = normalize(sub3f(light->center, previous_ray_result.intersection_point)); // scale3f(ray_result.surface_normal, dot3f(ray_result.surface_normal, ray.direction));
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
        dir_tan0 = normalize(dir_tan0);
        Vec3f dir_tan1 = cross3f(dir_normal, dir_tan0);

        Vec3f point_on_sphere = add3f(scale3f(dir_normal, - point_on_hemi_sphere.z),
                                add3f(scale3f(dir_tan0, - point_on_hemi_sphere.y),
                                scale3f(dir_tan1, - point_on_hemi_sphere.x)));

        point_on_sphere = scale3f(normalize(point_on_sphere), light->radius);

        Vec3f target_on_light = add3f(point_on_sphere, light->center);

        Vec3f light_path = sub3f(target_on_light, previous_ray_result.intersection_point);
        Vec3f light_direction = normalize(light_path);

        float scale = dot3f(light_direction, previous_ray_result.surface_normal);
        if (scale < 0) continue;

        Ray ray = {
            .origin = previous_ray_result.intersection_point,
            .direction = light_direction,
        };

        Ray_Result ray_result = find_ray_hit(scene, ray, norm(light_path));

        if (ray_result.hit) continue;

        float power = light->albedo / (4.0f * M_PI * norm2(light_path));
        Vec3f final_light_color = scale3f(light->color, power * scale);

        contribution.color = add3f(contribution.color, final_light_color);
        contribution.light_hit += 1;
    }

    return contribution;
}

Vec3f raytrace_ray(Job *job, Scene *scene, Ray ray)
{
    Camera *camera = &scene->camera;

    float light_hit_count = 0;
    Vec3f color_scale = WHITE;
    Vec3f accumulated_color = BLACK;

    // Keeps tracks of whether we are inside a transparent object or not.
    // With this method it is impossible to have transparent objects intersecting.
    bool within_transparent_object = false;

    for (int ray_count = 0; ray_count < 64; ray_count++) {
        Ray_Result ray_result = find_ray_hit(scene, ray, FLT_MAX);

        if (!ray_result.hit) break;

        Material *material = scene->materials + ray_result.material_id;

        if (ray_result.object_type == LIGHT) {
            light_hit_count += 1;

            Light *light = scene->lights + ray_result.object_id;
            float dist = ray_result.t_min;
            float scale = light->albedo / (4.0f * M_PI * ray_result.t_min * ray_result.t_min);

            Vec3f light_color = scale3f(light->color, scale);

            accumulated_color = add3f(accumulated_color, mul3f(color_scale, light_color));

            // Lorsque qu'on atteint une lumière de force, il n'y a plus de réflexions.
            break;
        } else if (material->mirror) {
            color_scale = mul3f(color_scale, material->mirror_color);

            Vec3f incident_normal = project3f(ray_result.surface_normal, ray.direction);
            Vec3f incident_tangent = sub3f(ray.direction, incident_normal);

            Vec3f new_dir = sub3f(incident_tangent, scale3f(incident_normal, 2.0f));

            ray.origin = ray_result.intersection_point;
            ray.direction = new_dir;
        } else if (material->transparent) {
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
            Light_Contribution light_contribution = compute_light_contribution(job, scene, ray, ray_result);
            Vec3f light_color = light_contribution.color;
            light_hit_count += light_contribution.light_hit;

            Vec3f diffuse_color = material->diffuse_color;
            if (material->textured) {
                diffuse_color = sample_texture(material->texture, ray_result.triangle_uv);
            }

            accumulated_color = add3f(accumulated_color, mul3f(color_scale, mul3f(light_color, diffuse_color)));

            color_scale = mul3f(color_scale, diffuse_color);

            // Si on ajoute plus assez de lumière après chaque rebond, on arrête
            if (norm2(color_scale) < 1e-6) break;

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
            dir_tan0 = normalize(dir_tan0);
            Vec3f dir_tan1 = cross3f(dir_normal, dir_tan0);

            Vec3f new_dir = add3f(scale3f(dir_normal, random_direction.z),
                            add3f(scale3f(dir_tan0, - random_direction.y),
                                  scale3f(dir_tan1, - random_direction.x)));

            new_dir = normalize(new_dir);

            color_scale = scale3f(color_scale, dot3f(new_dir, dir_normal));

            // Si on ajoute plus assez de lumière après chaque rebond, on arrête
            if (norm2(color_scale) < 1e-6) break;

            ray.direction = new_dir;
            // ray.direction = add3f(ray.direction, scale3f(new_dir, -2.0f));
        }
    }

    if (light_hit_count == 0) {
        return BLACK;
    }

    return scale3f(accumulated_color, 1 / light_hit_count);
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

const    uint32_t REGION_SPLIT_SIZE = 16;
volatile uint32_t region_count = 0;
volatile uint32_t finished = 0;

void thread_entry_point(Job *job)
{
    for (;;) {
        uint32_t region = __sync_fetch_and_add(&region_count, 1);

        if (region >= REGION_SPLIT_SIZE * REGION_SPLIT_SIZE) return;

        uint32_t x = region % REGION_SPLIT_SIZE;
        uint32_t y = region / REGION_SPLIT_SIZE;

        job->region.point.x = x * job->region.size.x;
        job->region.point.y = y * job->region.size.y;

        raytrace_region(job);

        uint32_t finished_regions = __sync_fetch_and_add(&finished, 1);
        printf("Finished %d regions out of %d.\n", finished_regions+1, REGION_SPLIT_SIZE * REGION_SPLIT_SIZE);
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

    Image3f texture = load_texture("models/Australian_Cattle_Dog_v1_L3.123c9c6a5764-399b-4e86-9897-6bcb08b5e8ed/Australian_Cattle_Dog_dif.jpg");

    int width = 512, height = 512;

    Image3f dest_image = create_blank_image(width, height);

    Material materials[] = {
        {   // The sphere in the center of the scene.
            .diffuse_color = scale3f(WHITE, 0.5f),
            .mirror = false,
        },
        {   // The background wall.
            .diffuse_color = GREEN,
            .mirror = false,
        },
        {   // The ground.
            .diffuse_color = BLUE,
            .mirror = false,
        },
        {   // Left wall.
            .diffuse_color = PINK,
            .mirror = false,
        },
        {   // Right Wall.
            .diffuse_color = YELLOW,
            .mirror = false,
        },
        {   // Perfect mirror sphere.
            .diffuse_color = BLACK,
            .mirror = true,
            .mirror_color = WHITE,
        },
        {   // Ceilling.
            .diffuse_color = RED,
            .mirror = false,
        },
        {   // Transparent sphere.
            .diffuse_color = BLACK,
            .mirror = false,

            .transparent = true,
            .n = 1.5f,
        },
        {
            .diffuse_color = WHITE,
            .mirror = false,

            .textured = true,
            .texture  = &texture,
        },
    };

    Sphere spheres[] = {
#if 0
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

    Light lights[] = {
        {
            .center = {-10.0f, 20.0f, 75.0f},
            .radius = 3.0f,
            .color = WHITE,
            .albedo = 1e5f,
        },
        {
            .center = {10, 17, 75},
            .radius = 2.0f,
            .color = {1, 0.1, 0.1},
            .albedo = 1e4f,
        },
    };

    Mesh dog_mesh = load_mesh("models/dog.obj");

    Scene_Mesh meshes[] = {
        {
            .mesh = &dog_mesh,
            .material_id = 5, // Chien mirror, parce que pourquoi pas ! (Cela crée plein de caustiques cependant !)
            .transform = {
                .position = {-25, -10, 0},
            
                .rotX = {-1, 0, 0},
                .rotY = {0, 1, 0},
                .rotZ = {0, 0, 1},
            },
        },
        {
            .mesh = &dog_mesh,
            .material_id = 8,
            .transform = {
                .position = {25, -10, 0},
            
                .rotX = {1, 0, 0},
                .rotY = {0, 1, 0},
                .rotZ = {0, 0, 1},
            },
        },
    };
    
    

    Scene scene = {
        .camera = {
            .position = {0.0f, 10.0f, 85.0f},
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

        .light_count = ARRAY_SIZE(lights),
        .lights = lights,

        .mesh_count = ARRAY_SIZE(meshes),
        .meshes = meshes,
    };

    if (USE_THREADS) {
        pthread_t threads[16] = {};
        Job jobs[16] = {};

        uint32_t block_width = width / 4;
        uint32_t block_height = height / 4;

#if 0
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                jobs[4 * i + j].scene = &scene;
                jobs[4 * i + j].image = &dest_image;
                jobs[4 * i + j].region = (Rect2i){
                    .point = {i * block_width, j * block_height},
                    .size = {block_width, block_height}
                };
                jobs[4 * i + j].ray_per_pixel = 64;
                jobs[4 * i + j].serie = create_random_serie();

                int code = pthread_create(&threads[4 * i + j], NULL, (void*)raytrace_region, &jobs[4 * i + j]);
                if (0 != code) {
                    printf("Error creating thread.\n");
                    exit(-1);
                }
            }
        }
#else
        for (int i = 0; i < 16; i++) {
            jobs[i].scene = &scene;
            jobs[i].image = &dest_image;
            jobs[i].region.size.x = width / REGION_SPLIT_SIZE;
            jobs[i].region.size.y = height / REGION_SPLIT_SIZE;
            jobs[i].ray_per_pixel = 2048;
            jobs[i].serie = create_random_serie();

            int code = pthread_create(&threads[i], NULL, (void*)thread_entry_point, &jobs[i]);
            if (0 != code) {
                printf("Error creating thread.\n");
                exit(-1);
            }
        }
#endif

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

#define FAST_OBJ_IMPLEMENTATION
#include "fast_obj.h"
