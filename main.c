#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FILE_NAME_LENGTH 256    // Max filename length
#define MAX_DEPTH 10                // Maximum recursion depth
#define EPSILON 0.01                // Intersections off-set
int debug = 0;                      // Flag for console debugging

typedef struct {
    double x, y, z;
} Vector;

typedef struct {
    Vector origin;
    Vector direction;
} Ray;

typedef struct {
    double r, g, b;
} Colour;

typedef struct {
    Vector center;
    double radius;
    Colour diffuse;
    Colour specular;
    Colour ambient;
    double reflection;
    double shininess;
    double refraction;
    double indexOfRefraction;
} Sphere;

typedef struct {
    Vector normal;
    double distance;
    Colour diffuse;
    Colour specular;
    Colour ambient;
    double reflection;
    double shininess;
    double refraction;
    double indexOfRefraction;
} Plane;

typedef struct {
    Vector T1, T2, T3;
    Vector normal;
    double distance;
    Colour diffuse;
    Colour specular;
    Colour ambient;
    double reflection;
    double shininess;
    double refraction;
    double indexOfRefraction;
} Triangle;

typedef struct {
    Vector direction;
    Colour diffuse;
    Colour specular;
} Light;

typedef enum {
    SPHERE,
    PLANE,
    TRIANGLE
} ObjectType;

typedef struct {
    ObjectType type;
    union {
        Sphere sphere;
        Plane plane;
        Triangle triangle;
    };
} Object;

Vector cameraPos;
int width, height;
Light lightSource;
double ambientIntensity;
int useSky;
Colour bgColour;
int numObjects;
Object* objects;

// Function for sky gradient background
Colour skyBackground(Ray ray) {
    double t = 0.5 * (ray.direction.y + 1.0);
    Colour result = (Colour) { 
        (1.0 - t) * 1.0 + t * 0.5, 
        (1.0 - t) * 1.0 + t * 0.7,
        (1.0 - t) * 1.0 + t * 1.0
    };
    //if (debug) printf("Sky gradient: (%f, %f, %f)\n", result.r, result.g, result.b);
    return result;
}

// Function to return appropriate background colour
Colour backgroundColour(Ray ray) {
    if (useSky == 1) {
        return skyBackground(ray);
    }
    else {
        return bgColour;
    }
}

// Function to add a value to a vector
Vector add2(double value, Vector v) {
    return (Vector) { v.x + value, v.y + value, v.z + value };
}

// Function to add two vectors
Vector add(Vector a, Vector b) {
    return (Vector) { a.x + b.x, a.y + b.y, a.z + b.z };
}

// Function to compute the cross product of two vectors
Vector cross(Vector a, Vector b) {
    return (Vector) { a.y* b.z - a.z * b.y, a.z* b.x - a.x * b.z, a.x* b.y - a.y * b.x };
}

// Function to compute the dot product of two vectors
double dot(Vector a, Vector b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Function to subtract two vectors
Vector subtract(Vector a, Vector b) {
    return (Vector) { a.x - b.x, a.y - b.y, a.z - b.z };
}

// Function to compute the length of a vector
double length(Vector v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Function to normalise a vector
Vector normalize(Vector v) {
    double len = length(v);
    if (len != 0.0) {
        return (Vector) { v.x / len, v.y / len, v.z / len };
    }
    return v;
}

// Function for checking for object reflection
int hasRefl(int hitobj) {
    switch (objects[hitobj].type) {
    case SPHERE:
        return objects[hitobj].sphere.reflection > 0.0;
    case PLANE:
        return objects[hitobj].plane.reflection > 0.0;
    case TRIANGLE:
        return objects[hitobj].triangle.reflection > 0.0;
    }
}

// Function to get the reflection coefficient for an object
double getReflectionCoefficient(int hitobj) {
    switch (objects[hitobj].type) {
    case SPHERE:
        return objects[hitobj].sphere.reflection;
    case PLANE:
        return objects[hitobj].plane.reflection;
    case TRIANGLE:
        return objects[hitobj].triangle.reflection;
    }
}

// Function for checking for object transmission
int hasTrans(int hitobj) {
    switch (objects[hitobj].type) {
    case SPHERE:
        return objects[hitobj].sphere.refraction > 0.0;
    case PLANE:
        return objects[hitobj].plane.refraction > 0.0;
    case TRIANGLE:
        return objects[hitobj].triangle.refraction > 0.0;
    }
}

// Function to get the refraction coefficient for an object
double getRefractionCoefficient(int hitobj) {
    switch (objects[hitobj].type) {
    case SPHERE:
        return objects[hitobj].sphere.refraction;
    case PLANE:
        return objects[hitobj].plane.refraction;
    case TRIANGLE:
        return objects[hitobj].triangle.refraction;
    }
}

// Function to get the index of refraction for an object
double getIndexOfRefraction(int hitobj) {
    switch (objects[hitobj].type) {
    case SPHERE:
        return objects[hitobj].sphere.indexOfRefraction;
    case PLANE:
        return objects[hitobj].plane.indexOfRefraction;
    case TRIANGLE:
        return objects[hitobj].triangle.indexOfRefraction;
    }
}

// Function to invert a vector
Vector invert(Vector v) {
    return (Vector) { -v.x, -v.y, -v.z };
}

// Function to multiply a vector by a scalar
Vector scalar(double scalar, Vector v) {
    return (Vector) { scalar * v.x, scalar * v.y, scalar * v.z };
}

// Function to multiply a vector by a scalar
Colour scaleColour(double scalar, Colour c) {
    return (Colour) { scalar* c.r, scalar* c.g, scalar* c.b };
}

// Function to get user selected input file
char* getInputFileFromUser() {
    char defaultFileName[] = "scenedata.txt";
    // Allocate memory for filename
    char* fileName = (char*)malloc(MAX_FILE_NAME_LENGTH * sizeof(char));

    if (fileName == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    int choice;

    printf("RAY TRACER\n");
    printf("Enter 1 for default scene input file or 0 to specify filename: ");
    if (scanf("%d", &choice) != 1) {
        // Handle input error
        printf("Error reading input.\n");
        exit(1);
    }

    // Clear the input buffer to handle newline characters
    int c;
    while ((c = getchar()) != '\n' && c != EOF);

    if (choice) {
        // Use the default filename
        strcpy(fileName, defaultFileName);
    }
    else {
        // Get custom filename
        printf("Enter the name of the input text file (ex.'scenedata.txt'): ");

        if (fgets(fileName, MAX_FILE_NAME_LENGTH, stdin) == NULL) {
            // Handle input error
            printf("Error reading input.\n");
            free(fileName);
            exit(1);
        }

        // Remove possible newline character
        size_t len = strlen(fileName);
        if (len > 0 && fileName[len - 1] == '\n') {
            fileName[len - 1] = '\0';
        }
    }
    printf("\n");
    return fileName;
}

// Function to print object data to console
void printObject(const Object* obj) {
    switch (obj->type) {
    case SPHERE:
        printf("SPHERE: centre(%.2lf, %.2lf, %.2lf), radius=%.2lf, diffuse(%.2lf, %.2lf, %.2lf), specular(%.2lf, %.2lf, %.2lf), ambient(%.2lf, %.2lf, %.2lf), reflection=%.2lf, shininess=%.2lf, refraction=%.2lf, indexOfRefraction=%.2lf\n",
            obj->sphere.center.x, obj->sphere.center.y, obj->sphere.center.z,
            obj->sphere.radius,
            obj->sphere.diffuse.r, obj->sphere.diffuse.g, obj->sphere.diffuse.b,
            obj->sphere.specular.r, obj->sphere.specular.g, obj->sphere.specular.b,
            obj->sphere.ambient.r, obj->sphere.ambient.g, obj->sphere.ambient.b,
            obj->sphere.reflection, obj->sphere.shininess, obj->sphere.refraction,
            obj->sphere.indexOfRefraction);
        break;

    case PLANE:
        printf("PLANE: normal(%.2lf, %.2lf, %.2lf), distance=%.2lf, diffuse(%.2lf, %.2lf, %.2lf), specular(%.2lf, %.2lf, %.2lf), ambient(%.2lf, %.2lf, %.2lf), reflection=%.2lf, shininess=%.2lf, refraction=%.2lf, indexOfRefraction=%.2lf\n",
            obj->plane.normal.x, obj->plane.normal.y, obj->plane.normal.z,
            obj->plane.distance,
            obj->plane.diffuse.r, obj->plane.diffuse.g, obj->plane.diffuse.b,
            obj->plane.specular.r, obj->plane.specular.g, obj->plane.specular.b,
            obj->plane.ambient.r, obj->plane.ambient.g, obj->plane.ambient.b,
            obj->plane.reflection, obj->plane.shininess, obj->plane.refraction,
            obj->plane.indexOfRefraction);
        break;

    case TRIANGLE:
        printf("TRIANGLE: T1(%.2lf, %.2lf, %.2lf), T2(%.2lf, %.2lf, %.2lf), T3(%.2lf, %.2lf, %.2lf), normal(%.2lf, %.2lf, %.2lf), distance=%.2lf, diffuse(%.2lf, %.2lf, %.2lf), specular(%.2lf, %.2lf, %.2lf), ambient(%.2lf, %.2lf, %.2lf), reflection=%.2lf, shininess=%.2lf, refraction=%.2lf, indexOfRefraction=%.2lf\n",
            obj->triangle.T1.x, obj->triangle.T1.y, obj->triangle.T1.z,
            obj->triangle.T2.x, obj->triangle.T2.y, obj->triangle.T2.z,
            obj->triangle.T3.x, obj->triangle.T3.y, obj->triangle.T3.z,
            obj->triangle.normal.x, obj->triangle.normal.y, obj->triangle.normal.z,
            obj->triangle.distance,
            obj->triangle.diffuse.r, obj->triangle.diffuse.g, obj->triangle.diffuse.b,
            obj->triangle.specular.r, obj->triangle.specular.g, obj->triangle.specular.b,
            obj->triangle.ambient.r, obj->triangle.ambient.g, obj->triangle.ambient.b,
            obj->triangle.reflection, obj->triangle.shininess, obj->triangle.refraction,
            obj->triangle.indexOfRefraction);
        break;

    }
}

// Function to print scene data to console for user
void printSceneData() {
    printf("SCENE DATA\n");
    printf("Number of Objects: %d\n", numObjects);

    printf("Camera Position: (%lf, %lf, %lf)\n", cameraPos.x, cameraPos.y, cameraPos.z);
    printf("Image properties: %dx%d\n", width, height);

    printf("Light Source Direction: (%lf, %lf, %lf)\n", lightSource.direction.x, lightSource.direction.y, lightSource.direction.z);
    printf("Light Source Diffuse Colour: (%lf, %lf, %lf)\n", lightSource.diffuse.r, lightSource.diffuse.g, lightSource.diffuse.b);
    printf("Light Source Specular Colour: (%lf, %lf, %lf)\n", lightSource.specular.r, lightSource.specular.g, lightSource.specular.b);

    printf("Ambient Intensity: %lf\n", ambientIntensity);
    printf("Use Sky Gradient Background (1-Yes; 0-No): %d\n", useSky);
    printf("Background Colour: (%lf, %lf, %lf)\n", bgColour.r, bgColour.g, bgColour.b);

    // Print each object
    for (int i = 0; i < numObjects; i++) {
        printObject(&objects[i]);
    }
    printf("\n");
}

// Function to read scene data from file
void readSceneData(const char* filename, int* numObjects, Vector* cameraPos, Light* lightSource, double* ambientIntensity, Colour* bgColor) {
    FILE* file;
    if ((file = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return;
    }

    // Read the number of objects
    fscanf(file, "%d", numObjects);

    // Read camera position
    fscanf(file, "%lf %lf %lf", &cameraPos->x, &cameraPos->y, &cameraPos->z);

    // Read image properties
    fscanf(file, "%d %d", &width, &height);

    // Read light source
    fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &lightSource->direction.x, &lightSource->direction.y, &lightSource->direction.z,
        &lightSource->diffuse.r, &lightSource->diffuse.g, &lightSource->diffuse.b,
        &lightSource->specular.r, &lightSource->specular.g, &lightSource->specular.b);

    // Read ambient intensity
    fscanf(file, "%lf", ambientIntensity);

    // Read if sky background
    fscanf(file, "%d", &useSky);

    // Read background color
    fscanf(file, "%lf %lf %lf", &bgColor->r, &bgColor->g, &bgColor->b);

    // Allocate memory for objects
    objects = (Object*)malloc(sizeof(Object) * (*numObjects));
    if (objects == NULL) {
        fprintf(stderr, "Memory allocation error for objects.\n");
        fclose(file);
        return;
    }

    // Read object data
    for (int i = 0; i < (*numObjects); i++) {
        int type;
        fscanf(file, "%d", &type);

        objects[i].type = (ObjectType)type;

        switch (objects[i].type) {
        case SPHERE:
            fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &objects[i].sphere.center.x, &objects[i].sphere.center.y, &objects[i].sphere.center.z,
                &objects[i].sphere.radius,
                &objects[i].sphere.diffuse.r, &objects[i].sphere.diffuse.g, &objects[i].sphere.diffuse.b,
                &objects[i].sphere.specular.r, &objects[i].sphere.specular.g, &objects[i].sphere.specular.b,
                &objects[i].sphere.ambient.r, &objects[i].sphere.ambient.g, &objects[i].sphere.ambient.b,
                &objects[i].sphere.reflection, &objects[i].sphere.shininess,
                &objects[i].sphere.refraction, &objects[i].sphere.indexOfRefraction);
            break;

        case PLANE:
            fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &objects[i].plane.normal.x, &objects[i].plane.normal.y, &objects[i].plane.normal.z,
                &objects[i].plane.distance,
                &objects[i].plane.diffuse.r, &objects[i].plane.diffuse.g, &objects[i].plane.diffuse.b,
                &objects[i].plane.specular.r, &objects[i].plane.specular.g, &objects[i].plane.specular.b,
                &objects[i].plane.ambient.r, &objects[i].plane.ambient.g, &objects[i].plane.ambient.b,
                &objects[i].plane.reflection, &objects[i].plane.shininess,
                &objects[i].plane.refraction, &objects[i].plane.indexOfRefraction);
            break;

        case TRIANGLE:
            fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                &objects[i].triangle.T1.x, &objects[i].triangle.T1.y, &objects[i].triangle.T1.z,
                &objects[i].triangle.T2.x, &objects[i].triangle.T2.y, &objects[i].triangle.T2.z,
                &objects[i].triangle.T3.x, &objects[i].triangle.T3.y, &objects[i].triangle.T3.z,
                &objects[i].triangle.normal.x, &objects[i].triangle.normal.y, &objects[i].triangle.normal.z,
                &objects[i].triangle.distance,
                &objects[i].triangle.diffuse.r, &objects[i].triangle.diffuse.g, &objects[i].triangle.diffuse.b,
                &objects[i].triangle.specular.r, &objects[i].triangle.specular.g, &objects[i].triangle.specular.b,
                &objects[i].triangle.ambient.r, &objects[i].triangle.ambient.g, &objects[i].triangle.ambient.b,
                &objects[i].triangle.reflection, &objects[i].triangle.shininess,
                &objects[i].triangle.refraction, &objects[i].triangle.indexOfRefraction);
            break;
        }
    }
    fclose(file);
}

// Function for checking ray-sphere intersection
int RaySphereIntersection(Ray ray, Sphere sphere, double* t, Vector* hitpt, Vector* hitnorm) {
    Vector oc = subtract(ray.origin, sphere.center);
    double a = dot(ray.direction, ray.direction);
    double b = 2.0 * dot(oc, ray.direction);
    double c = dot(oc, oc) - sphere.radius * sphere.radius;

    double discriminant = b * b - 4 * a * c;

    if (discriminant > 0) {
        double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
        double t2 = (-b + sqrt(discriminant)) / (2.0 * a);

        if (t1 > 0) {
            *t = t1;
        }
        else if (t2 > 0) {
            *t = t2;
        }
        else {
            return 0;   // No intersection
        }

        // Calculate the intersection point
        *hitpt = (Vector){
            ray.origin.x + *t * ray.direction.x,
            ray.origin.y + *t * ray.direction.y,
            ray.origin.z + *t * ray.direction.z
        };

        *hitnorm = normalize(subtract(*hitpt, sphere.center));

        return 1;   // Intersection found
    }
    else {
        return 0;   // No intersection
    }
}

// Function for checking ray-plane intersection
int RayPlaneIntersection(Ray ray, Plane plane, double* t, Vector* hitpt, Vector* hitnorm) {
    double denominator = dot(plane.normal, ray.direction);

    // Ensure that the ray and plane are not parallel
    if (fabs(denominator) > 1e-6) {
        *t = -(dot(ray.origin, plane.normal) + plane.distance) / denominator;

        if (*t >= 0) {
            // Calculate the intersection point
            *hitpt = (Vector){
                ray.origin.x + *t * ray.direction.x,
                ray.origin.y + *t * ray.direction.y,
                ray.origin.z + *t * ray.direction.z
            };

            *hitnorm = normalize(plane.normal);

            return 1;  // Intersection found
        }
    }

    return 0;  // No intersection
}

// Function to test if point inside triangle
// code from https://blackpawn.com/texts/pointinpoly/
int isPointInsideTriangle(Vector P, Vector T1, Vector T2, Vector T3) {
    // Calculate vectors
    Vector edge1 = subtract(T2, T1);
    Vector edge2 = subtract(T3, T1);
    Vector P_minus_T1 = subtract(P, T1);

    // Calculate dot products
    double dot11 = dot(edge1, edge1);
    double dot12 = dot(edge1, edge2);
    double dot22 = dot(edge2, edge2);
    double dotP1 = dot(P_minus_T1, edge1);
    double dotP2 = dot(P_minus_T1, edge2);

    // Calculate barycentric coordinates
    double invDenom = 1.0 / (dot11 * dot22 - dot12 * dot12);
    double alpha = (dot22 * dotP1 - dot12 * dotP2) * invDenom;
    double beta = (dot11 * dotP2 - dot12 * dotP1) * invDenom;

    // Check if P is inside the triangle
    return (alpha >= 0.0 && beta >= 0.0 && alpha + beta <= 1.0);
}

// Function to calculate normal and distance for triangle
void calculateTriangleNormal(Triangle* triangle) {
    Vector edge1 = subtract(triangle->T1, triangle->T2);
    Vector edge2 = subtract(triangle->T3, triangle->T1);
    triangle->normal = normalize(cross(edge1, edge2));
    triangle->distance = -dot(triangle->normal, triangle->T1);
}

// Function for checking ray-triangle intersection
int RayTriangleIntersection(Ray ray, Triangle triangle, double* t, Vector* hitpt, Vector* hitnorm) {
    double denominator = dot(triangle.normal, ray.direction);

    // Ensure that the ray and plane are not parallel
    if (fabs(denominator) > 1e-6) {
        *t = -(dot(ray.origin, triangle.normal) + triangle.distance) / denominator;

        if (*t > 0) {
            // Calculate the intersection point
            *hitpt = (Vector){
                ray.origin.x + *t * ray.direction.x,
                ray.origin.y + *t * ray.direction.y,
                ray.origin.z + *t * ray.direction.z
            };

            *hitnorm = normalize(triangle.normal);

            // Check if the intersection point is inside the triangle
            if (isPointInsideTriangle(*hitpt, triangle.T1, triangle.T2, triangle.T3)) {
                return 1;  // Intersection found
            }
        }
    }

    return 0;  // No intersection
}

// Function for checking ray-object intersection
int RayIntersection(Ray ray, int* hitobj, Vector* hitpt, Vector* hitnorm) {
    double closestT = INFINITY;

    for (int i = 0; i < numObjects; i++) {
        double t;
        Vector localHitpt, localHitnorm;

        switch (objects[i].type) {
        case SPHERE:
            if (RaySphereIntersection(ray, objects[i].sphere, &t, &localHitpt, &localHitnorm)) {
                if (t > EPSILON && t < closestT ) {
                    closestT = t;
                    *hitobj = i;
                    *hitpt = localHitpt;
                    *hitnorm = localHitnorm;
                }
            }
            break;

        case PLANE:
            if (RayPlaneIntersection(ray, objects[i].plane, &t, &localHitpt, &localHitnorm)) {
                if (t > EPSILON && t < closestT) {
                    closestT = t;
                    *hitobj = i;
                    *hitpt = localHitpt;
                    *hitnorm = localHitnorm;
                }
            }
            break;

        case TRIANGLE:
            if (RayTriangleIntersection(ray, objects[i].triangle, &t, &localHitpt, &localHitnorm)) {
                if (t > EPSILON && t < closestT) {
                    closestT = t;
                    *hitobj = i;
                    *hitpt = localHitpt;
                    *hitnorm = localHitnorm;
                }
            }
            break;
        }
    }

    return (closestT != INFINITY);
}

// Function to get the specular highlight coefficient for an object
double getSpecularExponent(int hitobj) {
    switch (objects[hitobj].type) {
    case SPHERE:
        return objects[hitobj].sphere.shininess;
    case PLANE:
        return objects[hitobj].plane.shininess;
    case TRIANGLE:
        return objects[hitobj].triangle.shininess;
    }
}

// Function to get the ambient surface coefficient for an object
Colour getAmbientSurfaceCoefficient(int hitobj) {
    switch (objects[hitobj].type) {
    case SPHERE:
        return objects[hitobj].sphere.ambient;
    case PLANE:
        return objects[hitobj].plane.ambient;
    case TRIANGLE:
        return objects[hitobj].triangle.ambient;
    }
}

// Function to get the diffuse coefficient for an object
Colour getDiffuseCoefficient(int hitobj) {
    switch (objects[hitobj].type) {
    case SPHERE:
        return objects[hitobj].sphere.diffuse;
    case PLANE:
        return objects[hitobj].plane.diffuse;
    case TRIANGLE:
        return objects[hitobj].triangle.diffuse;
    }
}

// Function to get the specular reflection coefficient for an object
Colour getSpecularCoefficient(int hitobj) {
    switch (objects[hitobj].type) {
    case SPHERE:
        return objects[hitobj].sphere.specular;
    case PLANE:
        return objects[hitobj].plane.specular;
    case TRIANGLE:
        return objects[hitobj].triangle.specular;
    }
}

// Function for checking shadow ray-object intersections
int RayShadowIntersection(Ray ray, int excludeObject) {
    double maxDist = length(subtract(lightSource.direction, ray.origin));

    for (int i = 0; i < numObjects; i++) {
        // Skip the excluded object
        if (i == excludeObject) {
            continue;
        }

        double t;
        Vector localHitpt, localHitnorm;

        switch (objects[i].type) {
        case SPHERE:
            if (RaySphereIntersection(ray, objects[i].sphere, &t, &localHitpt, &localHitnorm) && t > 0.0 && t < maxDist && objects[i].sphere.indexOfRefraction != 1) {
                return 1;  // In shadow
            }
            break;

        case PLANE:
            if (RayPlaneIntersection(ray, objects[i].plane, &t, &localHitpt, &localHitnorm) && t > 0.0 && t < maxDist && objects[i].plane.indexOfRefraction != 1) {
                return 1;  // In shadow
            }
            break;

        case TRIANGLE:
            if (RayTriangleIntersection(ray, objects[i].triangle, &t, &localHitpt, &localHitnorm) && t > 0.0 && t < maxDist && objects[i].triangle.indexOfRefraction != 1) {
                return 1;  // In shadow
            }
            break;

        }
    }

    return 0;  // Not in shadow
}

// Function to compute local lighting at hitpoint
Colour shade(int hitobj, Vector hitpt, Vector hitnorm) {
    // Ambient component
    Colour ambient = {
        getAmbientSurfaceCoefficient(hitobj).r * ambientIntensity,
        getAmbientSurfaceCoefficient(hitobj).g * ambientIntensity,
        getAmbientSurfaceCoefficient(hitobj).b * ambientIntensity
    };

    Vector lightDir = normalize(subtract(lightSource.direction,hitpt));
    double cosTheta = dot(normalize(hitnorm), lightDir);

    // Check for shadows
    Ray shadowRay = { add(hitpt,scalar(EPSILON, hitnorm)), lightDir};

    if (RayShadowIntersection(shadowRay, -1)) {
        return (Colour) { ambient.r, ambient.g, ambient.b };
    }

    if (cosTheta > 0) {
        // Diffuse reflection
        Colour diffuse = {
            getDiffuseCoefficient(hitobj).r * cosTheta,
            getDiffuseCoefficient(hitobj).g * cosTheta,
            getDiffuseCoefficient(hitobj).b * cosTheta
        };

        // Specular reflection
        Vector reflectDir = normalize(subtract(scalar(2.0 * cosTheta, hitnorm), lightDir));
        double cosB = dot(reflectDir, normalize(subtract(cameraPos, hitpt)));
        double specularTerm = pow(cosB, getSpecularExponent(hitobj));

        Colour specular = {
            lightSource.specular.r * getSpecularCoefficient(hitobj).r * specularTerm,
            lightSource.specular.g * getSpecularCoefficient(hitobj).g * specularTerm,
            lightSource.specular.b * getSpecularCoefficient(hitobj).b * specularTerm
        };

        return (Colour) {
            diffuse.r + ambient.r + specular.r,
                diffuse.g + ambient.g + specular.g,
                diffuse.b + ambient.b + specular.b
        };
    }
    else {
        return (Colour) { ambient.r, ambient.g, ambient.b }; // Only return ambient
    }
}

// Function for calculating reflection vector
Ray calcReflection(Ray ray, int hitobj, Vector hitpt, Vector hitnorm) {
    double cosI = dot(ray.direction, hitnorm);
    Vector reflVec = normalize(subtract(ray.direction, scalar(2.0 * cosI, hitnorm)));
    Vector reflHitpoint = add(hitpt, scalar(EPSILON, reflVec));
    return (Ray) { reflHitpoint, reflVec };
}

// Function for calculating transmission vector
Ray* calcTransmission(Ray ray, int hitobj, Vector hitpt, Vector hitnorm) {
    double eta_i, eta_t;  // Refractive indices of the two media

    // Determine the direction of the incident ray
    double cosI = dot(hitnorm, ray.direction);

    // Determine the refractive indices
    if (cosI > 0) {  
        // The ray is inside the object
        eta_i = 1.0;
        eta_t = getIndexOfRefraction(hitobj);
        hitnorm = invert(hitnorm);
    }
    else { 
        eta_i = getIndexOfRefraction(hitobj);
        eta_t = 1.0;
        cosI = -cosI;
    }

    double eta = eta_i / eta_t;
    double rootTerm = 1.0 + eta * eta * (cosI * cosI - 1.0);

    if (rootTerm < 0.0) {
        // Total internal reflection, no refraction
        if (debug) printf("Total internal reflection\n");
        return NULL;
    }
    else {
        // Calculate the transmission vector
        if (debug) printf("eta: %f/%f\n", eta_i, eta_t);
        Vector transVec = normalize(add(scalar(eta, ray.direction), scalar(eta * cosI - sqrt(rootTerm), hitnorm)));
        Vector transHitpt = add(hitpt, scalar(EPSILON, transVec));
        static Ray result;
        result = (Ray) { transHitpt, transVec };
        return &result;
    }
}

// Function for ray tracing
Colour RayTrace(Ray ray, int recDepth) {
    if (recDepth > MAX_DEPTH) {
        if (debug) printf("Max depth of %d reached, returning background colour.\n", recDepth - 1);
        return backgroundColour(ray);
    }
    else {
        int hitobj;
        Vector hitpt, hitnorm;

        if (RayIntersection(ray, &hitobj, &hitpt, &hitnorm)) {
            if (debug) printf("Intersection found with object %d\n", hitobj);
            if (debug) printf("Hitpoint: (%f, %f, %f)\n", hitpt.x, hitpt.y, hitpt.z);
            if (debug) printf("Normal_%d: (%f, %f, %f)\n", recDepth, hitnorm.x, hitnorm.y, hitnorm.z);
            Colour localColour = shade(hitobj, hitpt, hitnorm);

            Colour reflectColour = { 0.0, 0.0, 0.0 };
            if (hasRefl(hitobj)) {
                Ray reflRay = calcReflection(ray, hitobj, hitpt, hitnorm);                
                if (debug) printf("Reflection ray_%d origin: (%f, %f, %f)\n", recDepth, reflRay.origin.x, reflRay.origin.y, reflRay.origin.z);
                if (debug) printf("Reflection ray_%d direction : (% f, % f, % f)\n", recDepth, reflRay.direction.x, reflRay.direction.y, reflRay.direction.z);
                reflectColour = RayTrace(reflRay, recDepth + 1);
                double reflectionCoefficient = getReflectionCoefficient(hitobj);
                // Attenuate the reflection colour by the reflection coefficient
                reflectColour = scaleColour(reflectionCoefficient, reflectColour);
            }

            Colour transColour = { 0.0, 0.0, 0.0 };
            if (hasTrans(hitobj)) {
                Ray *transRay = calcTransmission(ray, hitobj, hitpt, hitnorm);
                if (transRay != NULL) {
                    if (debug) printf("Transmission ray_%d origin: (%f, %f, %f)\n", recDepth, transRay->origin.x, transRay->origin.y, transRay->origin.z);
                    if (debug) printf("Transmission ray_%d direction: (%f, %f, %f)\n", recDepth, transRay->direction.x, transRay->direction.y, transRay->direction.z);
                    transColour = RayTrace(*transRay, recDepth + 1);
                    // Attenuate the refraction colour by the refraction coefficient                
                    double refractionCoefficient = getRefractionCoefficient(hitobj);
                    transColour = scaleColour(refractionCoefficient, transColour);
                }
            }

            // Combine local, reflective, and refractive colours, clamped between 0 and 1
            Colour finalColour = {
                fmin(1.0, localColour.r + reflectColour.r + transColour.r),
                fmin(1.0, localColour.g + reflectColour.g + transColour.g),
                fmin(1.0, localColour.b + reflectColour.b + transColour.b)
            };

            return finalColour;
        }
        else {
            if (debug) printf("No intersection found, returning background colour.\n");
            return backgroundColour(ray);
        }
    }
}

// Function for pixel-wise debugging
void shootRayAtPixel(int pixelX, int pixelY) {
    double u = (double)pixelX / (width - 1);
    double v = (double)pixelY / (height - 1);

    Vector pixel = { 2 * u - 1, 1 - 2 * v, 1 };
    pixel.x *= width / (double)height;

    Ray ray = { cameraPos, pixel };
    ray.direction = normalize(ray.direction);

    printf("Ray origin: (%f, %f, %f)\n", ray.origin.x, ray.origin.y, ray.origin.z);
    printf("Ray direction: (%f, %f, %f)\n", ray.direction.x, ray.direction.y, ray.direction.z);

    Colour colour = RayTrace(ray, 1);
}

// Function to set up view and generate image file
void generateImage() {
    FILE* fp;

    // Open file
    if ((fp = fopen("output.ppm", "w")) == NULL) {
        fprintf(stderr, "Error opening the file.\n");
        return;
    }

    // Write PPM header
    fprintf(fp, "P3\n%d %d\n255\n", width, height);

    // Generate the image
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            // Normalize pixel coordinates to the range [0, 1]
            double u = (double)i / (width);
            double v = (double)j / (height);

            // Map normalized coordinates to the screen space
            Vector pixel = { 2 * u - 1, 1 - 2 * v, 1 };
            pixel.x *= width / (double)height;  // scale x based on aspect ratio

            // Generate a ray for each pixel
            Ray ray = { cameraPos, pixel };
            ray.direction = normalize(ray.direction);

            // Trace the ray and get the colour
            Colour colour = RayTrace(ray, 1);

            // Scale the colour values to [0, 255] and write to file
            fprintf(fp, "%d %d %d ", (int)(colour.r * 255), (int)(colour.g * 255), (int)(colour.b * 255));
        }
        fprintf(fp, "\n");
    }
    // Close the file
    fclose(fp);
    printf("Image file output.ppm saved in root application directory\n");
}

int main() {
    // Specify scene data input file
    char* fileName = getInputFileFromUser();

    // Read scene data
    readSceneData(fileName, &numObjects, &cameraPos, &lightSource, &ambientIntensity, &bgColour);

    // Calculate normal and distance for each triangle
    for (int i = 0; i < numObjects; i++) {
        if (objects[i].type == TRIANGLE) {
            calculateTriangleNormal(&objects[i].triangle);
        }
    }
    // Print scene data to console
    printSceneData(); 
   
    //debug = 1;
     
    if (debug) {
        shootRayAtPixel(400, 183);
    }
    else {
        generateImage();
    }

    free(objects); 
    free(fileName);

    return 0;
}
