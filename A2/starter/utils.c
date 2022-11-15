/*
   utils.c - F.J. Estrada, Dec. 9, 2010

   Utilities for the ray tracer. You will need to complete
   some of the functions in this file. Look for the sections
   marked "TO DO". Be sure to read the rest of the file and
   understand how the entire code works.

   HOWEVER: Note that there are a lot of incomplete functions
            that will only be used for the advanced ray tracer!
      So, read the handout carefully and implement only
      the code you need for the corresponding assignment.

   Last updated: Aug. 2017  -  F.J.E.
*/

/*****************************************************************************
 * COMPLETE THIS TEXT BOX:
 *
 * 1) Student Name: Alexia Tolentino
 * 2) Student Name: Reem Al Halabi
 *
 * 1) Student number: 1006293967
 * 2) Student number: 1006321927
 *
 * 1) UtorID tolent55
 * 2) UtorID alhala16
 *
 * We hereby certify that the work contained here is our own
 *
 * Alexia Tolentino                  Reem Al Halabi
 * (sign with your name)            (sign with your name)
 *
 ********************************************************************************/

#include "utils.h"

// A useful 4x4 identity matrix which can be used at any point to
// initialize or reset object transformations
double eye4x4[4][4] = {{1.0, 0.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0, 0.0},
                       {0.0, 0.0, 1.0, 0.0},
                       {0.0, 0.0, 0.0, 1.0}};

/////////////////////////////////////////////
// Primitive data structure section
/////////////////////////////////////////////
struct point3D *newPoint(double px, double py, double pz)
{
  // Allocate a new point structure, initialize it to
  // the specified coordinates, and return a pointer
  // to it.

  struct point3D *pt = (struct point3D *)calloc(1, sizeof(struct point3D));
  if (!pt)
    fprintf(stderr, "Out of memory allocating point structure!\n");
  else
  {
    pt->px = px;
    pt->py = py;
    pt->pz = pz;
    pt->pw = 1.0;
  }
  return (pt);
}

struct pointLS *newPLS(struct point3D *p0, double r, double g, double b)
{
  // Allocate a new point light sourse structure. Initialize the light
  // source to the specified RGB colour
  // Note that this is a point light source in that it is a single point
  // in space, if you also want a uniform direction for light over the
  // scene (a so-called directional light) you need to place the
  // light source really far away.

  struct pointLS *ls = (struct pointLS *)calloc(1, sizeof(struct pointLS));
  if (!ls)
    fprintf(stderr, "Out of memory allocating light source!\n");
  else
  {
    memcpy(&ls->p0, p0, sizeof(struct point3D)); // Copy light source location

    ls->col.R = r; // Store light source colour and
    ls->col.G = g; // intensity
    ls->col.B = b;
  }
  return (ls);
}

/////////////////////////////////////////////
// Ray and normal transforms
/////////////////////////////////////////////
inline void rayTransform(struct ray3D *ray_orig, struct ray3D *ray_transformed, struct object3D *obj)
{
  // Transforms a ray using the inverse transform for the specified object. This is so that we can
  // use the intersection test for the canonical object. Note that this has to be done carefully!

  ///////////////////////////////////////////
  // TO DO: Complete this function
  ///////////////////////////////////////////

  // Setting transformed ray to be original ray
  ray_transformed->p0 = ray_orig->p0;
  ray_transformed->d = ray_orig->d;

  // Transforming the ray's point by objects inverse transform
  matVecMult(obj->Tinv, &ray_transformed->p0);

  // Constructing the inverse transpose matrix
  double Tinv_l[4][4];
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      Tinv_l[i][j] = obj->Tinv[i][j];
      Tinv_l[i][3] = 0;
      Tinv_l[3][i] = 0;
    }
  }
  Tinv_l[3][3] = 1;

  // Transforming by linear part of matrix
  matVecMult(Tinv_l, &ray_transformed->d);
}
inline void normalTransform(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj)
{
  // Computes the normal at an affinely transformed point given the original normal and the
  // object's inverse transformation. From the notes:
  // n_transformed=A^-T*n normalized.

  ///////////////////////////////////////////
  // TO DO: Complete this function
  ///////////////////////////////////////////

  // Setting normal transformed to be the original normal
  *n_transformed = *n_orig;

  // Getting transpose of matrix
  double obj_T[4][4];

  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      obj_T[i][j] = obj->Tinv[j][i];
    }
  }

  // Applying transformation
  matVecMult(obj_T, n_transformed);

  // Normalize normal!
  normalize(n_transformed);
  n_transformed->pw = 1;
}

/////////////////////////////////////////////
// Object management section
/////////////////////////////////////////////
void insertObject(struct object3D *o, struct object3D **list)
{
  if (o == NULL)
    return;
  // Inserts an object into the object list.
  if (*(list) == NULL)
  {
    *(list) = o;
    (*(list))->next = NULL;
  }
  else
  {
    o->next = (*(list))->next;
    (*(list))->next = o;
  }
}

struct object3D *newPlane(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny)
{
  // Intialize a new plane with the specified parameters:
  // ra, rd, rs, rg - Albedos for the components of the Phong model
  // r, g, b, - Colour for this plane
  // alpha - Transparency, must be set to 1 unless you are doing refraction
  // r_index - Refraction index if you are doing refraction.
  // shiny - Exponent for the specular component of the Phong model
  //
  // The plane is defined by the following vertices (CCW)
  // (1,1,0), (-1,1,0), (-1,-1,0), (1,-1,0)
  // With normal vector (0,0,1) (i.e. parallel to the XY plane)

  struct object3D *plane = (struct object3D *)calloc(1, sizeof(struct object3D));

  if (!plane)
    fprintf(stderr, "Unable to allocate new plane, out of memory!\n");
  else
  {
    plane->alb.ra = ra;
    plane->alb.rd = rd;
    plane->alb.rs = rs;
    plane->alb.rg = rg;
    plane->col.R = r;
    plane->col.G = g;
    plane->col.B = b;
    plane->alpha = alpha;
    plane->r_index = r_index;
    plane->shinyness = shiny;
    plane->intersect = &planeIntersect;
    plane->surfaceCoords = &planeCoordinates;
    plane->randomPoint = &planeSample;
    plane->texImg = NULL;
    plane->photonMap = NULL;
    plane->normalMap = NULL;
    memcpy(&plane->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
    memcpy(&plane->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
    plane->textureMap = &texMap;
    plane->frontAndBack = 1;
    plane->photonMapped = 0;
    plane->normalMapped = 0;
    plane->isCSG = 0;
    plane->isLightSource = 0;
    plane->CSGnext = NULL;
    plane->next = NULL;
  }
  return (plane);
}

struct object3D *newSphere(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double r_index, double shiny)
{
  // Intialize a new sphere with the specified parameters:
  // ra, rd, rs, rg - Albedos for the components of the Phong model
  // r, g, b, - Colour for this plane
  // alpha - Transparency, must be set to 1 unless you are doing refraction
  // r_index - Refraction index if you are doing refraction.
  // shiny -Exponent for the specular component of the Phong model
  //
  // This is assumed to represent a unit sphere centered at the origin.
  //

  struct object3D *sphere = (struct object3D *)calloc(1, sizeof(struct object3D));

  if (!sphere)
    fprintf(stderr, "Unable to allocate new sphere, out of memory!\n");
  else
  {
    sphere->alb.ra = ra;
    sphere->alb.rd = rd;
    sphere->alb.rs = rs;
    sphere->alb.rg = rg;
    sphere->col.R = r;
    sphere->col.G = g;
    sphere->col.B = b;
    sphere->alpha = alpha;
    sphere->r_index = r_index;
    sphere->shinyness = shiny;
    sphere->intersect = &sphereIntersect;
    sphere->surfaceCoords = &sphereCoordinates;
    sphere->randomPoint = &sphereSample;
    sphere->texImg = NULL;
    sphere->photonMap = NULL;
    sphere->normalMap = NULL;
    memcpy(&sphere->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
    memcpy(&sphere->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
    sphere->textureMap = &texMap;
    sphere->frontAndBack = 0;
    sphere->photonMapped = 0;
    sphere->normalMapped = 0;
    sphere->isCSG = 0;
    sphere->isLightSource = 0;
    sphere->CSGnext = NULL;
    sphere->next = NULL;
  }
  return (sphere);
}

struct object3D *newCyl(double ra, double rd, double rs, double rg, double r, double g, double b, double alpha, double R_index, double shiny)
{
  ///////////////////////////////////////////////////////////////////////////////////////
  // TO DO:
  //	Complete the code to create and initialize a new cylinder object.
  ///////////////////////////////////////////////////////////////////////////////////////
  struct object3D *cyl = (struct object3D *)calloc(1, sizeof(struct object3D));
  if (!cyl)
    fprintf(stderr, "Unable to allocate new cylinder, out of memory!\n");
  else
  {
    cyl->alb.ra = ra;
    cyl->alb.rd = rd;
    cyl->alb.rs = rs;
    cyl->alb.rg = rg;
    cyl->col.R = r;
    cyl->col.G = g;
    cyl->col.B = b;
    cyl->alpha = alpha;
    cyl->r_index = R_index;
    cyl->shinyness = shiny;
    cyl->intersect = &cylIntersect;
    cyl->surfaceCoords = &cylCoordinates;
    cyl->randomPoint = &cylSample;
    cyl->texImg = NULL;
    cyl->photonMap = NULL;
    cyl->normalMap = NULL;
    memcpy(&cyl->T[0][0], &eye4x4[0][0], 16 * sizeof(double));
    memcpy(&cyl->Tinv[0][0], &eye4x4[0][0], 16 * sizeof(double));
    cyl->textureMap = &texMap;
    cyl->frontAndBack = 0;
    cyl->photonMapped = 0;
    cyl->normalMapped = 0;
    cyl->isCSG = 0;
    cyl->isLightSource = 0;
    cyl->CSGnext = NULL;
    cyl->next = NULL;
  }
  return (cyl);
}

///////////////////////////////////////////////////////////////////////////////////////
// TO DO:
//	Complete the functions that compute intersections for the canonical plane
//      and canonical sphere with a given ray. This is the most fundamental component
//      of the raytracer.
///////////////////////////////////////////////////////////////////////////////////////
void planeIntersect(struct object3D *plane, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
  // Computes and returns the value of 'lambda' at the intersection
  // between the specified ray and the specified canonical plane.

  /////////////////////////////////
  // TO DO: Complete this function.
  /////////////////////////////////
  // Clone of the ray to avoid modifying original
  struct ray3D rayclone;
  struct point3D origin;
  struct point3D og_n;
  rayclone.p0 = ray->p0;
  rayclone.d = ray->d;

  // Origin
  origin.px = 0;
  origin.py = 0;
  origin.pz = 0;
  origin.pw = 1;

  // Normal Vector
  og_n.px = 0;
  og_n.py = 0;
  og_n.pz = 1;
  og_n.pw = 1;

  // Intersection point
  struct point3D intersection;

  // Transform the ray wrt plane's transformations
  rayTransform(ray, &rayclone, plane);

  // p-a stored in origin
  subVectors(&rayclone.p0, &origin);
  *lambda = dot(&origin, &og_n) / dot(&rayclone.d, &og_n);

  // Update the ray's intersection
  rayPosition(&rayclone, *lambda, &intersection);

  // Check if modified intersection in plane
  if (abs(intersection.px) >= 1 || abs(intersection.py) >= 1)
  {
    *lambda = -1;
    return;
  }

  // IF INTERSECTION IN VIEW PLANE?
  // Plane is in view plane
  if (*lambda > 0)
  {
    // Set normal wrt plane
    normalTransform(&og_n, n, plane);
    normalize(n);

    // Set a and b wrt plane
    *a = (intersection.px + 1) / 2;
    *b = (intersection.py + 1) / 2;

    // Update ray position wrt found lambda
    rayPosition(ray, *lambda, p);
  }
  // Plane is not in the viewplane
  else
  {
    *lambda = -1;
    *a = -1;
    *b = -1;
  }
}

void sphereIntersect(struct object3D *sphere, struct ray3D *ray, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
  // Computes and returns the value of 'lambda' at the intersection
  // between the specified ray and the specified canonical sphere.

  /////////////////////////////////
  // TO DO: Complete this function.
  /////////////////////////////////

  // Initialize lambda
  *lambda = -1;

  // Set normal
  struct point3D og_n;
  og_n.px = 0;
  og_n.py = 0;
  og_n.pz = 0;
  og_n.pw = 1;

  // Set origin
  struct point3D origin;
  origin.px = 0;
  origin.py = 0;
  origin.pz = 0;
  origin.pw = 1;

  // Clone of the ray to avoid modifying original
  struct ray3D rayclone;
  rayclone.p0 = ray->p0;
  rayclone.d = ray->d;

  // Transform ray wrt sphere's transformation
  rayTransform(ray, &rayclone, sphere);

  // Set variables for quadratic formula
  double A = dot(&rayclone.d, &rayclone.d);
  double B = dot(&rayclone.p0, &rayclone.d);
  double C = dot(&rayclone.p0, &rayclone.p0) - 1;
  double disc = pow(B, 2) - (A * C);

  // No solutions
  if (disc < 0)
  {
    *lambda = -1;
    return;
  }
  // One solution
  else if (disc == 0)
  {
    *lambda = ((-1) * B + sqrt(disc)) / A;
  }
  // Two solutions
  else if (disc > 0)
  {
    // Setting lambda1 and lambda2
    double lambda1 = ((-1) * B + sqrt(disc)) / A;
    double lambda2 = ((-1) * B - sqrt(disc)) / A;

    // Checking is lambda is in viewplane
    // Assigning smallest lambda
    if (lambda1 < 0 && lambda2 < 0)
    {
      *lambda = -1;
    }
    else if (lambda1 > 0 && lambda2 < 0)
    {
      *lambda = lambda1;
    }
    else if (lambda1 < 0 && lambda2 > 0)
    {
      *lambda = lambda2;
    }
    else
    {
      *lambda = ((lambda1 < lambda2) ? (lambda1) : (lambda2));
    }
  }

  // If assigned lambda is valid
  if (*lambda > 0)
  {
    // Update ray position wrt new lambda
    rayPosition(ray, *lambda, p);
    rayPosition(&rayclone, *lambda, &og_n);

    // Normalizing and setting normal
    normalTransform(&og_n, n, sphere);
    normalize(n);

    // Setting a and b
    *a = acos(p->pz) / (2 * PI);
    *b = (atan(p->py / p->px) + PI / 2) / PI;
  }
  else
  {
    *lambda = -1;
    *a = -1;
    *b = -1;
  }
}

void cylIntersect(struct object3D *cylinder, struct ray3D *r, double *lambda, struct point3D *p, struct point3D *n, double *a, double *b)
{
  // Computes and returns the value of 'lambda' at the intersection
  // between the specified ray and the specified canonical cylinder.

  /////////////////////////////////
  // TO DO: Complete this function.
  /////////////////////////////////

  // Setting lambda = -1 to indicate no changes made
  *lambda = -1;

  // Declaring necessary objects
  struct ray3D rayclone;
  struct point3D circ_norm;
  rayclone.p0 = r->p0;
  rayclone.d = r->d;

  // Normal Vector
  circ_norm.px = 0;
  circ_norm.py = 0;
  circ_norm.pz = 0;
  circ_norm.pw = 1;

  // Transform the ray wrt the cylinder
  rayTransform(r, &rayclone, cylinder);

  // FINDING INTERSECTION W QUADRATIC WALL (unit circles st z is in the range)
  // Using quadratic formula
  double A = pow(rayclone.d.px, 2) + pow(rayclone.d.py, 2);
  double B = rayclone.p0.px * rayclone.d.px + rayclone.p0.py * rayclone.d.py;
  double C = pow(rayclone.p0.px, 2) + pow(rayclone.p0.py, 2) - 1;
  double disc = pow(B, 2) - (A * C);

  struct point3D intersected;
  double z;
  double opt_lambda;
  double lambda1;
  double lambda2;
  // Checking number of solutions using discriminant
  // TWO SOLUTIONS
  if (disc > 0)
  {
    lambda1 = ((-1) * B + sqrt(disc)) / A;
    lambda2 = ((-1) * B - sqrt(disc)) / A;

    // Checking is lambda is in viewplane
    // Assigning smallest lambda
    if (lambda1 < 0 && lambda2 < 0)
    {
      *lambda = -1;
    }
    else if (lambda1 > 0 && lambda2 < 0)
    {
      // Check if z for lambda 1 is within the height of unit cyl
      z = rayclone.p0.pz + lambda1 * rayclone.d.pz;
      if (abs(z) <= 1 && lambda1 > 0)
      { // if height is valid and we're in viewplane
        // See if this is the smallest lambda
        // Or if lambda hasn't been assigned
        if ((lambda1 < *lambda) || (*lambda == -1))
        {
          *lambda = lambda1;
        }
      }
    }
    else if (lambda1 < 0 && lambda2 > 0)
    {
      // Check if z for lambda 2 is within the height of unit cyl
      z = rayclone.p0.pz + lambda2 * rayclone.d.pz;
      if (abs(z) <= 1 && lambda2 > 0)
      { // if height is valid and we're in viewplane
        // See if this is the smallest lambda in comparison to 1
        // Or if lambda hasn't been assigned
        if ((lambda2 < *lambda) || (*lambda == -1))
        {
          *lambda = lambda2;
        }
      }
    }
    else
    {
      opt_lambda = ((lambda1 < lambda2) ? (lambda1) : (lambda2));
      // Check if z for lambda 2 is within the height of unit cyl
      z = rayclone.p0.pz + opt_lambda * rayclone.d.pz;
      if (abs(z) <= 1 && opt_lambda > 0)
      { // if height is valid and we're in viewplane
        // See if this is the smallest lambda in comparison to 1
        // Or if lambda hasn't been assigned
        if ((*lambda == -1))
        {
          *lambda = opt_lambda;
        }
      }
    }
  }
  // ONE SOLUTION
  else if (disc == 0)
  {
    opt_lambda = (-1 * B + sqrt(disc)) / (A);

    // Check if z for opt_lambda is within the height of unit cyl
    z = rayclone.p0.pz + opt_lambda * rayclone.d.pz;
    if (abs(z) <= 1 && opt_lambda > 0)
    { // if height is valid and we're in viewplane
      // See if this is the smallest lambda
      // Or if lambda hasn't been assigned
      if ((opt_lambda < *lambda) || (*lambda == -1))
      {
        *lambda = opt_lambda;
      }
    }
  } // Otherwsie NO SOLUTIONS so lambda stays -1
  else
  {
    *lambda = -1;
  }

  // If assigned lambda is valid
  if (*lambda > 0)
  {
    // Update ray position wrt new lambda
    rayPosition(r, *lambda, p);
    rayPosition(&rayclone, *lambda, &circ_norm);

    // Normalizing and setting normal
    normalTransform(&circ_norm, n, cylinder);
    normalize(n);

    // Setting a and b
    *a = acos(p->pz) / (2 * PI);
    *b = (atan(p->py / p->px) + PI / 2) / PI;
    return;
  }
  else
  {
    *lambda = -1;
    *a = -1;
    *b = -1;
  }

  // FIND INTERSECTION BETWEEN BASES OF CYL (plane st z = 0 or 1)
  // Declaring necessary objects
  struct point3D origin, og_n, intersection;

  // Origin
  origin.px = 0;
  origin.py = 0;
  origin.pz = 0;
  origin.pw = 1;

  // CHECK INTERSECTION FOR PLANE 1
  lambda1 = (1 - rayclone.p0.pz) / rayclone.d.pz;
  // Update the ray's position according to lambda of first plane
  rayPosition(&rayclone, lambda1, &intersection);

  // Check if modified pos in circular plane: if x^2 + y^2 <= 1
  if (pow(intersection.px, 2) + pow(intersection.py, 2) <= 1)
  {
    // See if this is the smallest lambda
    // Or if lambda hasn't been assigned
    if ((lambda1 < *lambda) || (*lambda == -1))
    {
      *lambda = lambda1;
    }
  }

  // CHECK INTERSECTION FOR PLANE 2
  lambda2 = (-1 - rayclone.p0.pz) / rayclone.d.pz;
  // Update the ray's position according to lambda of first plane
  rayPosition(&rayclone, lambda2, &intersection);

  // Check if modified pos in circular plane: if x^2 + y^2 <= 1
  if (pow(intersection.px, 2) + pow(intersection.py, 2) <= 1)
  {
    // See if this is the smallest lambda
    // Or if lambda hasn't been assigned
    if ((lambda2 < *lambda) || (*lambda == -1))
    {
      *lambda = lambda2;
    }
  }

  // Normal Vector Update before transform
  og_n.px = 0;
  og_n.py = 0;
  og_n.pz = rayclone.p0.pz + ((*lambda) * (rayclone.d.pz));
  og_n.pw = 1;
  normalize(&og_n);

  // CHECKING IF FINAL LAMBDA IN VIEW PLANE
  if (*lambda > 0)
  {
    // Plane is in the viewplane
    p->px = intersection.px;
    p->py = intersection.py;
    p->pz = intersection.pz;
    p->pw = intersection.pw;

    // Update position due to smallest lambda
    rayPosition(&rayclone, *lambda, &intersection);

    // Update normal
    normalTransform(&og_n, n, cylinder);

    // Setting a and b
    *a = (atan(intersection.py / intersection.px) + PI / 2) / PI;
    *b = intersection.pz;
  }
  else
  {
    // plane is not in the viewplane
    *lambda = -1;
    *a = -1;
    *b = -1;
  }
}

/////////////////////////////////////////////////////////////////
// Surface coordinates & random sampling on object surfaces
/////////////////////////////////////////////////////////////////
void planeCoordinates(struct object3D *plane, double a, double b, double *x, double *y, double *z)
{
  // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b in [0,1].
  // 'a' controls displacement from the left side of the plane, 'b' controls displacement from the
  // bottom of the plane.
  struct point3D original; // Original starts in the original left boottm corner
  original.px = -1 + 2 * a;
  original.py = -1 + 2 * b;
  original.pz = 0;
  original.pw = 1;

  // plane's coordinates (after transformed)
  matVecMult(plane->T, &original);

  // Setting x,y,z
  *x = original.px;
  *y = original.py;
  *z = original.pz;
}

void sphereCoordinates(struct object3D *sphere, double a, double b, double *x, double *y, double *z)
{
  // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
  // 'a' in [0, 2*PI] corresponds to the spherical coordinate theta
  // 'b' in [-PI/2, PI/2] corresponds to the spherical coordinate phi

  /////////////////////////////////
  // TO DO: Complete this function.
  /////////////////////////////////

  // from center we move to the transformed position
  struct point3D center;
  center.px = sin(b) * cos(a);
  center.py = sin(b) * sin(a);
  center.pz = cos(b);
  // sphere's coordinates (after transformed)
  matVecMult(sphere->T, &center);

  // starting from center of the circle (assuming radius is 1) we adjust based on the spherical angles
  *x = center.px;
  *y = center.py;
  *z = center.pz;
}

void cylCoordinates(struct object3D *cyl, double a, double b, double *x, double *y, double *z)
{
  // Return in (x,y,z) the coordinates of a point on the plane given by the 2 parameters a,b.
  // 'a' in [0, 2*PI] corresponds to angle theta around the cylinder
  // 'b' in [0, 1] corresponds to height from the bottom

  // Getting starting point
  struct point3D *original;
  original->px = cos(a);
  original->py = sin(a);
  original->pz = b;

  // cylinder's coordinates (after transformed)
  matVecMult(cyl->T, original);

  // Setting x,y,z
  *x = original->px;
  *y = original->py;
  *z = original->pz;
}

void planeSample(struct object3D *plane, double *x, double *y, double *z)
{
  // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the plane
  // Sapling should be uniform, meaning there should be an equal change of gedtting
  // any spot on the plane

  // Get random value for a and b
  double a = drand48();
  double b = drand48();

  // send random a and b to the planeCoordinates function
  planeCoordinates(plane, a, b, x, y, z);
}

void sphereSample(struct object3D *sphere, double *x, double *y, double *z)
{
  // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the sphere
  // Sampling should be uniform - note that this is tricky for a sphere, do some
  // research and document in your report what method is used to do this, along
  // with a reference to your source.

  // Generate random angles
  double a = drand48() * 2 * PI;            // get random angle from 0 to 2pi
  double b = (drand48() - 0.5) * PI;        // get random angle from -pi to pi
  sphereCoordinates(sphere, a, b, x, y, z); // get adjusted according to a and b
}

void cylSample(struct object3D *cyl, double *x, double *y, double *z)
{
  // Returns the 3D coordinates (x,y,z) of a randomly sampled point on the cylinder
  // Sampling should be uniform over the cylinder.

  // Get random value for a and b
  double a = drand48() * 2 * PI;
  double b = drand48();

  // send random a and b to the cylCoordinates function
  cylCoordinates(cyl, a, b, x, y, z);
}

/////////////////////////////////
// Texture mapping functions
/////////////////////////////////
void loadTexture(struct object3D *o, const char *filename, int type, struct textureNode **t_list)
{
  // Load a texture or normal map image from file and assign it to the
  // specified object.
  // type:   1  ->  Texture map  (RGB, .ppm)
  //         2  ->  Normal map   (RGB, .ppm)
  //         3  ->  Alpha map    (grayscale, .pgm)
  // Stores loaded images in a linked list to avoid replication
  struct image *im;
  struct textureNode *p;

  if (o != NULL)
  {
    // Check current linked list
    p = *(t_list);
    while (p != NULL)
    {
      if (strcmp(&p->name[0], filename) == 0)
      {
        // Found image already on the list
        if (type == 1)
          o->texImg = p->im;
        else if (type == 2)
          o->normalMap = p->im;
        else
          o->alphaMap = p->im;
        return;
      }
      p = p->next;
    }

    // Load this texture image
    if (type == 1 || type == 2)
      im = readPPMimage(filename);
    else if (type == 3)
      im = readPGMimage(filename);

    // Insert it into the texture list
    if (im != NULL)
    {
      p = (struct textureNode *)calloc(1, sizeof(struct textureNode));
      strcpy(&p->name[0], filename);
      p->type = type;
      p->im = im;
      p->next = NULL;
      // Insert into linked list
      if ((*(t_list)) == NULL)
        *(t_list) = p;
      else
      {
        p->next = (*(t_list))->next;
        (*(t_list))->next = p;
      }
      // Assign to object
      if (type == 1)
        o->texImg = im;
      else if (type == 2)
        o->normalMap = im;
      else
        o->alphaMap = im;
    }

  } // end if (o != NULL)
}

void texMap(struct image *img, double a, double b, double *R, double *G, double *B)
{
  /*
   Function to determine the colour of a textured object at
   the normalized texture coordinates (a,b).

   a and b are texture coordinates in [0 1].
   img is a pointer to the image structure holding the texture for
    a given object.

   The colour is returned in R, G, B. Uses bi-linear interpolation
   to determine texture colour.
  */

  //////////////////////////////////////////////////
  // TO DO (Assignment 4 only):
  //
  //  Complete this function to return the colour
  // of the texture image at the specified texture
  // coordinates. Your code should use bi-linear
  // interpolation to obtain the texture colour.
  //////////////////////////////////////////////////

  // Makes sure the coordinates are within 
 if (a >= 0 && a <= 1 && b >= 0 && b <= 1) {
  
  // Calulate the width and height scaling factors
  int width_factor = img->sx - 1;
  int height_factor = img->sy - 1;

  double u = a;
  double v = b;


  // Put texture colour on the object
  struct colourRGB p1, p2, p3, p4;
  double *rgbTex;
  rgbTex = (double *) img->rgbdata;

  double x = a * width_factor;
  double y = b * height_factor;

  int x_tl = floor(x);
  int y_tl = floor(y);

  int x_tr = ceil(x);
  int y_tr = floor(y);

  int x_bl = floor(x);
  int y_bl = ceil(y);

  int x_br = ceil(x);
  int y_br = ceil(y);

  p1.R = rgbTex[3 * (y_tl * height_factor + (width_factor - x_tl))];
  p1.G = rgbTex[3 * (y_tl * height_factor + (width_factor - x_tl)) + 1];
  p1.B = rgbTex[3 * (y_tl * height_factor + (width_factor - x_tl)) + 2];
  p2.R = rgbTex[3 * (y_tr * height_factor + (width_factor - x_tr))];
  p2.G = rgbTex[3 * (y_tr * height_factor + (width_factor - x_tr)) + 1];
  p2.B = rgbTex[3 * (y_tr * height_factor + (width_factor - x_tr)) + 2];
  p3.R = rgbTex[3 * (y_bl * height_factor + (width_factor - x_bl))];
  p3.G = rgbTex[3 * (y_bl * height_factor + (width_factor - x_bl)) + 1];
  p3.B = rgbTex[3 * (y_bl * height_factor + (width_factor - x_bl)) + 2];
  p4.R = rgbTex[3 * (y_br * height_factor + (width_factor - x_br))];
  p4.G = rgbTex[3 * (y_br * height_factor + (width_factor - x_br)) + 1];
  p4.B = rgbTex[3 * (y_br * height_factor + (width_factor - x_br)) + 2];

  // Bilinear interpolation
  double a_top, a_bot;
  // R
  a_top = u * p2.R + (1 - u) * p1.R;
  a_bot = u * p4.R + (1 - u) * p3.R;
  *(R) = v * a_bot + (1 - v) * a_top;
  // G
  a_top = u * p2.G + (1 - u) * p1.G;
  a_bot = u * p4.G + (1 - u) * p3.G;
  *(G) = v * a_bot + (1 - v) * a_top;
  // B
  a_top = u * p2.B + (1 - u) * p1.B;
  a_bot = u * p4.B + (1 - u) * p3.B;
  *(B) = v * a_bot + (1 - v) * a_top;
 }
  return;
}

void alphaMap(struct image *img, double a, double b, double *alpha)
{
  // Just like texture map but returns the alpha value at a,b,
  // notice that alpha maps are single layer grayscale images, hence
  // the separate function.

  //////////////////////////////////////////////////
  // TO DO (Assignment 4 only):
  //
  //  Complete this function to return the alpha
  // value from the image at the specified texture
  // coordinates. Your code should use bi-linear
  // interpolation to obtain the texture colour.
  //////////////////////////////////////////////////

  int a_im, b_im;
  double a_frac, b_frac;
  double *ptr_im=(double *)img->rgbdata; 

  if (a<0) a=0; 
  if (a>1) a=1;
  if (b<0) b=0;
  if (b>1) b=1;  

  a_frac=(a*img->sx);
  b_frac=(b*img->sy);
  
  a_im=(int)a_frac;
  b_im=(int)b_frac;

 *(alpha)=*(ptr_im+a_im+b_im);	// Returns 1 which means fully opaque. Replace
 return;	// with your code if implementing alpha maps.
}

/////////////////////////////
// Light sources
/////////////////////////////
void insertPLS(struct pointLS *l, struct pointLS **list)
{
  if (l == NULL)
    return;
  // Inserts a light source into the list of light sources
  if (*(list) == NULL)
  {
    *(list) = l;
    (*(list))->next = NULL;
  }
  else
  {
    l->next = (*(list))->next;
    (*(list))->next = l;
  }
}

void addAreaLight(double sx, double sy, double nx, double ny, double nz,
                  double tx, double ty, double tz, int N,
                  double r, double g, double b, struct object3D **o_list, struct pointLS **l_list, int obj_type)
{
  /*
    This function sets up and inserts a rectangular area light source
    with size (sx, sy)
    orientation given by the normal vector (nx, ny, nz)
    centered at (tx, ty, tz)
    consisting of (N) point light sources (uniformly sampled)
    and with colour/intensity (r,g,b)

    Note that the light source must be visible as a uniformly colored rectangle which
    casts no shadows. If you require a lightsource to shade another, you must
    make it into a proper solid box with a back and sides of non-light-emitting
    material
  */

  /////////////////////////////////////////////////////
  // TO DO: (Assignment 4!)
  // Implement this function to enable area light sources
  /////////////////////////////////////////////////////

  // NOTE: The best way to implement area light sources is to random sample from the
  //       light source's object surface within rtShade(). This is a bit more tricky
  //       but reduces artifacts significantly. If you do that, then there is no need
  //       to insert a series of point lightsources in this function.
  struct object3D *areaLS;

  // Initiate all 
  if (obj_type == 1){
    areaLS=newPlane(1,1,1,1,r,g,b,1,1,0);
  }
  else if (obj_type == 2){
    areaLS = newSphere(1,1,1,1,r,g,b,1,1,0);
  }
  else{
    areaLS = newCyl(1,1,1,1,r,g,b,1,1,0);
  }


  areaLS->isLightSource = 1;
  Scale(areaLS,sx,sy,1);
  Translate(areaLS,tx,ty,tz);
  RotateX(areaLS,nx);
  RotateY(areaLS,ny);
  RotateZ(areaLS,nz);
  invert(*areaLS->T,*areaLS->Tinv);	
  insertObject(areaLS,o_list);

  //Origin of the lightsource
  for (int i = 0; i < N; i++) {
    struct point3D origin;
    origin.pw = 1;   

    if (obj_type == 1) {
      planeSample(areaLS, &origin.px, &origin.py, &origin.pz);
    } else if (obj_type == 2) {
      sphereSample(areaLS, &origin.px, &origin.py, &origin.pz);
    } else{
      cylSample(areaLS, &origin.px, &origin.py, &origin.pz);
    }
      
    struct pointLS *l;
    l=newPLS(&origin,r,g,b);
    insertPLS(l,l_list);
  }
}

///////////////////////////////////
// Geometric transformation section
///////////////////////////////////

void invert(double *T, double *Tinv)
{
  // Computes the inverse of transformation matrix T.
  // the result is returned in Tinv.

  double *U, *s, *V, *rv1;
  int singFlag, i;

  // Invert the affine transform
  U = NULL;
  s = NULL;
  V = NULL;
  rv1 = NULL;
  singFlag = 0;

  SVD(T, 4, 4, &U, &s, &V, &rv1);
  if (U == NULL || s == NULL || V == NULL)
  {
    fprintf(stderr, "Error: Matrix not invertible for this object, returning identity\n");
    memcpy(Tinv, eye4x4, 16 * sizeof(double));
    return;
  }

  // Check for singular matrices...
  for (i = 0; i < 4; i++)
    if (*(s + i) < 1e-9)
      singFlag = 1;
  if (singFlag)
  {
    fprintf(stderr, "Error: Transformation matrix is singular, returning identity\n");
    memcpy(Tinv, eye4x4, 16 * sizeof(double));
    return;
  }

  // Compute and store inverse matrix
  InvertMatrix(U, s, V, 4, Tinv);

  free(U);
  free(s);
  free(V);
}

void RotateXMat(double T[4][4], double theta)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that rotates the object theta *RADIANS* around the
  // X axis.

  double R[4][4];
  memset(&R[0][0], 0, 16 * sizeof(double));

  R[0][0] = 1.0;
  R[1][1] = cos(theta);
  R[1][2] = -sin(theta);
  R[2][1] = sin(theta);
  R[2][2] = cos(theta);
  R[3][3] = 1.0;

  matMult(R, T);
}

void RotateX(struct object3D *o, double theta)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that rotates the object theta *RADIANS* around the
  // X axis.

  double R[4][4];
  memset(&R[0][0], 0, 16 * sizeof(double));

  R[0][0] = 1.0;
  R[1][1] = cos(theta);
  R[1][2] = -sin(theta);
  R[2][1] = sin(theta);
  R[2][2] = cos(theta);
  R[3][3] = 1.0;

  matMult(R, o->T);
}

void RotateYMat(double T[4][4], double theta)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that rotates the object theta *RADIANS* around the
  // Y axis.

  double R[4][4];
  memset(&R[0][0], 0, 16 * sizeof(double));

  R[0][0] = cos(theta);
  R[0][2] = sin(theta);
  R[1][1] = 1.0;
  R[2][0] = -sin(theta);
  R[2][2] = cos(theta);
  R[3][3] = 1.0;

  matMult(R, T);
}

void RotateY(struct object3D *o, double theta)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that rotates the object theta *RADIANS* around the
  // Y axis.

  double R[4][4];
  memset(&R[0][0], 0, 16 * sizeof(double));

  R[0][0] = cos(theta);
  R[0][2] = sin(theta);
  R[1][1] = 1.0;
  R[2][0] = -sin(theta);
  R[2][2] = cos(theta);
  R[3][3] = 1.0;

  matMult(R, o->T);
}

void RotateZMat(double T[4][4], double theta)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that rotates the object theta *RADIANS* around the
  // Z axis.

  double R[4][4];
  memset(&R[0][0], 0, 16 * sizeof(double));

  R[0][0] = cos(theta);
  R[0][1] = -sin(theta);
  R[1][0] = sin(theta);
  R[1][1] = cos(theta);
  R[2][2] = 1.0;
  R[3][3] = 1.0;

  matMult(R, T);
}

void RotateZ(struct object3D *o, double theta)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that rotates the object theta *RADIANS* around the
  // Z axis.

  double R[4][4];
  memset(&R[0][0], 0, 16 * sizeof(double));

  R[0][0] = cos(theta);
  R[0][1] = -sin(theta);
  R[1][0] = sin(theta);
  R[1][1] = cos(theta);
  R[2][2] = 1.0;
  R[3][3] = 1.0;

  matMult(R, o->T);
}

void TranslateMat(double T[4][4], double tx, double ty, double tz)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that translates the object by the specified amounts.

  double tr[4][4];
  memset(&tr[0][0], 0, 16 * sizeof(double));

  tr[0][0] = 1.0;
  tr[1][1] = 1.0;
  tr[2][2] = 1.0;
  tr[0][3] = tx;
  tr[1][3] = ty;
  tr[2][3] = tz;
  tr[3][3] = 1.0;

  matMult(tr, T);
}

void Translate(struct object3D *o, double tx, double ty, double tz)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that translates the object by the specified amounts.

  double tr[4][4];
  memset(&tr[0][0], 0, 16 * sizeof(double));

  tr[0][0] = 1.0;
  tr[1][1] = 1.0;
  tr[2][2] = 1.0;
  tr[0][3] = tx;
  tr[1][3] = ty;
  tr[2][3] = tz;
  tr[3][3] = 1.0;

  matMult(tr, o->T);
}

void ScaleMat(double T[4][4], double sx, double sy, double sz)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that scales the object as indicated.

  double S[4][4];
  memset(&S[0][0], 0, 16 * sizeof(double));

  S[0][0] = sx;
  S[1][1] = sy;
  S[2][2] = sz;
  S[3][3] = 1.0;

  matMult(S, T);
}

void Scale(struct object3D *o, double sx, double sy, double sz)
{
  // Multiply the current object transformation matrix T in object o
  // by a matrix that scales the object as indicated.

  double S[4][4];
  memset(&S[0][0], 0, 16 * sizeof(double));

  S[0][0] = sx;
  S[1][1] = sy;
  S[2][2] = sz;
  S[3][3] = 1.0;

  matMult(S, o->T);
}

void printmatrix(double mat[4][4])
{
  fprintf(stderr, "Matrix contains:\n");
  fprintf(stderr, "%f %f %f %f\n", mat[0][0], mat[0][1], mat[0][2], mat[0][3]);
  fprintf(stderr, "%f %f %f %f\n", mat[1][0], mat[1][1], mat[1][2], mat[1][3]);
  fprintf(stderr, "%f %f %f %f\n", mat[2][0], mat[2][1], mat[2][2], mat[2][3]);
  fprintf(stderr, "%f %f %f %f\n", mat[3][0], mat[3][1], mat[3][2], mat[3][3]);
}

/////////////////////////////////////////
// Camera and view setup
/////////////////////////////////////////
struct view *setupView(struct point3D *e, struct point3D *g, struct point3D *up, double f, double wl, double wt, double wsize)
{
  /*
    This function sets up the camera axes and viewing direction as discussed in the
    lecture notes.
    e - Camera center
    g - Gaze direction
    up - Up vector
    fov - Fild of view in degrees
    f - focal length
  */
  struct view *c;
  struct point3D *u, *v;

  u = v = NULL;

  // Allocate space for the camera structure
  c = (struct view *)calloc(1, sizeof(struct view));
  if (c == NULL)
  {
    fprintf(stderr, "Out of memory setting up camera model!\n");
    return (NULL);
  }

  // Set up camera center and axes
  c->e.px = e->px; // Copy camera center location, note we must make sure
  c->e.py = e->py; // the camera center provided to this function has pw=1
  c->e.pz = e->pz;
  c->e.pw = 1;

  // Set up w vector (camera's Z axis). w=-g/||g||
  c->w.px = -g->px;
  c->w.py = -g->py;
  c->w.pz = -g->pz;
  c->w.pw = 1;
  normalize(&c->w);

  // Set up the horizontal direction, which must be perpenticular to w and up
  u = cross(&c->w, up);
  normalize(u);
  c->u.px = u->px;
  c->u.py = u->py;
  c->u.pz = u->pz;
  c->u.pw = 1;

  // Set up the remaining direction, v=(u x w)  - Mind the signs
  v = cross(&c->u, &c->w);
  normalize(v);
  c->v.px = v->px;
  c->v.py = v->py;
  c->v.pz = v->pz;
  c->v.pw = 1;

  // Copy focal length and window size parameters
  c->f = f;
  c->wl = wl;
  c->wt = wt;
  c->wsize = wsize;

  // Set up coordinate conversion matrices
  // Camera2World matrix (M_cw in the notes)
  // Mind the indexing convention [row][col]
  c->C2W[0][0] = c->u.px;
  c->C2W[1][0] = c->u.py;
  c->C2W[2][0] = c->u.pz;
  c->C2W[3][0] = 0;

  c->C2W[0][1] = c->v.px;
  c->C2W[1][1] = c->v.py;
  c->C2W[2][1] = c->v.pz;
  c->C2W[3][1] = 0;

  c->C2W[0][2] = c->w.px;
  c->C2W[1][2] = c->w.py;
  c->C2W[2][2] = c->w.pz;
  c->C2W[3][2] = 0;

  c->C2W[0][3] = c->e.px;
  c->C2W[1][3] = c->e.py;
  c->C2W[2][3] = c->e.pz;
  c->C2W[3][3] = 1;

  // World2Camera matrix (M_wc in the notes)
  // Mind the indexing convention [row][col]
  c->W2C[0][0] = c->u.px;
  c->W2C[1][0] = c->v.px;
  c->W2C[2][0] = c->w.px;
  c->W2C[3][0] = 0;

  c->W2C[0][1] = c->u.py;
  c->W2C[1][1] = c->v.py;
  c->W2C[2][1] = c->w.py;
  c->W2C[3][1] = 0;

  c->W2C[0][2] = c->u.pz;
  c->W2C[1][2] = c->v.pz;
  c->W2C[2][2] = c->w.pz;
  c->W2C[3][2] = 0;

  c->W2C[0][3] = -dot(&c->u, &c->e);
  c->W2C[1][3] = -dot(&c->v, &c->e);
  c->W2C[2][3] = -dot(&c->w, &c->e);
  c->W2C[3][3] = 1;

  free(u);
  free(v);
  return (c);
}

// helper function to set up view if cam is null
void nullSetupView(struct view *c, struct point3D *e, struct point3D *g, struct point3D *up, struct point3D u, struct point3D v, struct point3D w)
{
  /*
    This function sets up the camera axes and viewing direction as discussed in the
    lecture notes.
    e - Camera center
    g - Gaze direction
    up - Up vector
    fov - Fild of view in degrees
    f - focal length
  */

  // Setting up camera to world matrix for the camera
  c->C2W[0][0] = u.px;
  c->C2W[1][0] = u.py;
  c->C2W[2][0] = u.pz;
  c->C2W[3][0] = 0;
  c->C2W[0][1] = v.px;
  c->C2W[1][1] = v.py;
  c->C2W[2][1] = v.pz;
  c->C2W[3][1] = 0;
  c->C2W[0][2] = w.px;
  c->C2W[1][2] = w.py;
  c->C2W[2][2] = w.pz;
  c->C2W[3][2] = 0;
  c->C2W[0][3] = e->px;
  c->C2W[1][3] = e->py;
  c->C2W[2][3] = e->pz;
  c->C2W[3][3] = 1;

  // Setting up world to camera matrix by transposing about with e
  double mwc[4][4];
  c->W2C[0][0] = u.px;
  c->W2C[1][0] = v.px;
  c->W2C[2][0] = w.px;
  c->W2C[3][0] = 0;
  c->W2C[0][1] = u.py;
  c->W2C[1][1] = v.py;
  c->W2C[2][1] = w.py;
  c->W2C[3][1] = 0;
  c->W2C[0][2] = u.pz;
  c->W2C[1][2] = v.pz;
  c->W2C[2][2] = w.pz;
  c->W2C[3][2] = 0;
  c->W2C[0][3] = -e->px;
  c->W2C[1][3] = -e->py;
  c->W2C[2][3] = -e->pz;
  c->W2C[3][3] = 1;
}



/////////////////////////////////////////
// Image I/O section
/////////////////////////////////////////
struct image *readPPMimage(const char *filename)
{
  // Reads an image from a .ppm file. A .ppm file is a very simple image representation
  // format with a text header followed by the binary RGB data at 24bits per pixel.
  // The header has the following form:
  //
  // P6
  // # One or more comment lines preceded by '#'
  // 340 200
  // 255
  //
  // The first line 'P6' is the .ppm format identifier, this is followed by one or more
  // lines with comments, typically used to inidicate which program generated the
  // .ppm file.
  // After the comments, a line with two integer values specifies the image resolution
  // as number of pixels in x and number of pixels in y.
  // The final line of the header stores the maximum value for pixels in the image,
  // usually 255.
  // After this last header line, binary data stores the RGB values for each pixel
  // in row-major order. Each pixel requires 3 bytes ordered R, G, and B.
  //
  // NOTE: Windows file handling is rather crotchetty. You may have to change the
  //       way this file is accessed if the images are being corrupted on read
  //       on Windows.
  //
  // readPPMdata converts the image colour information to floating point. This is so that
  // the texture mapping function doesn't have to do the conversion every time
  // it is asked to return the colour at a specific location.
  //

  FILE *f;
  struct image *im;
  char line[1024];
  int sizx, sizy;
  int i;
  unsigned char *tmp;
  double *fRGB;
  int tmpi;
  char *tmpc;

  im = (struct image *)calloc(1, sizeof(struct image));
  if (im != NULL)
  {
    im->rgbdata = NULL;
    f = fopen(filename, "rb+");
    if (f == NULL)
    {
      fprintf(stderr, "Unable to open file %s for reading, please check name and path\n", filename);
      free(im);
      return (NULL);
    }
    tmpc = fgets(&line[0], 1000, f);
    if (strcmp(&line[0], "P6\n") != 0)
    {
      fprintf(stderr, "Wrong file format, not a .ppm file or header end-of-line characters missing\n");
      free(im);
      fclose(f);
      return (NULL);
    }
    fprintf(stderr, "%s\n", line);
    // Skip over comments
    tmpc = fgets(&line[0], 511, f);
    while (line[0] == '#')
    {
      fprintf(stderr, "%s", line);
      tmpc = fgets(&line[0], 511, f);
    }
    sscanf(&line[0], "%d %d\n", &sizx, &sizy); // Read file size
    fprintf(stderr, "nx=%d, ny=%d\n\n", sizx, sizy);
    im->sx = sizx;
    im->sy = sizy;

    tmpc = fgets(&line[0], 9, f); // Read the remaining header line
    fprintf(stderr, "%s\n", line);
    tmp = (unsigned char *)calloc(sizx * sizy * 3, sizeof(unsigned char));
    fRGB = (double *)calloc(sizx * sizy * 3, sizeof(double));
    if (tmp == NULL || fRGB == NULL)
    {
      fprintf(stderr, "Out of memory allocating space for image\n");
      free(im);
      fclose(f);
      return (NULL);
    }

    tmpi = fread(tmp, sizx * sizy * 3 * sizeof(unsigned char), 1, f);
    fclose(f);

    // Conversion to floating point
    for (i = 0; i < sizx * sizy * 3; i++)
      *(fRGB + i) = ((double)*(tmp + i)) / 255.0;
    free(tmp);
    im->rgbdata = (void *)fRGB;

    return (im);
  }

  fprintf(stderr, "Unable to allocate memory for image structure\n");
  return (NULL);
}

struct image *readPGMimage(const char *filename)
{
  // Just like readPPMimage() except it is used to load grayscale alpha maps. In
  // alpha maps, a value of 255 corresponds to alpha=1 (fully opaque) and 0
  // correspondst to alpha=0 (fully transparent).
  // A .pgm header of the following form is expected:
  //
  // P5
  // # One or more comment lines preceded by '#'
  // 340 200
  // 255
  //
  // readPGMdata converts the image grayscale data to double floating point in [0,1].

  FILE *f;
  struct image *im;
  char line[1024];
  int sizx, sizy;
  int i;
  unsigned char *tmp;
  double *fRGB;
  int tmpi;
  char *tmpc;

  im = (struct image *)calloc(1, sizeof(struct image));
  if (im != NULL)
  {
    im->rgbdata = NULL;
    f = fopen(filename, "rb+");
    if (f == NULL)
    {
      fprintf(stderr, "Unable to open file %s for reading, please check name and path\n", filename);
      free(im);
      return (NULL);
    }
    tmpc = fgets(&line[0], 1000, f);
    if (strcmp(&line[0], "P5\n") != 0)
    {
      fprintf(stderr, "Wrong file format, not a .pgm file or header end-of-line characters missing\n");
      free(im);
      fclose(f);
      return (NULL);
    }
    // Skip over comments
    tmpc = fgets(&line[0], 511, f);
    while (line[0] == '#')
      tmpc = fgets(&line[0], 511, f);
    sscanf(&line[0], "%d %d\n", &sizx, &sizy); // Read file size
    im->sx = sizx;
    im->sy = sizy;

    tmpc = fgets(&line[0], 9, f); // Read the remaining header line
    tmp = (unsigned char *)calloc(sizx * sizy, sizeof(unsigned char));
    fRGB = (double *)calloc(sizx * sizy, sizeof(double));
    if (tmp == NULL || fRGB == NULL)
    {
      fprintf(stderr, "Out of memory allocating space for image\n");
      free(im);
      fclose(f);
      return (NULL);
    }

    tmpi = fread(tmp, sizx * sizy * sizeof(unsigned char), 1, f);
    fclose(f);

    // Conversion to double floating point
    for (i = 0; i < sizx * sizy; i++)
      *(fRGB + i) = ((double)*(tmp + i)) / 255.0;
    free(tmp);
    im->rgbdata = (void *)fRGB;

    return (im);
  }

  fprintf(stderr, "Unable to allocate memory for image structure\n");
  return (NULL);
}

struct image *newImage(int size_x, int size_y)
{
  // Allocates and returns a new image with all zeros. Assumes 24 bit per pixel,
  // unsigned char array.
  struct image *im;

  im = (struct image *)calloc(1, sizeof(struct image));
  if (im != NULL)
  {
    im->rgbdata = NULL;
    im->sx = size_x;
    im->sy = size_y;
    im->rgbdata = (void *)calloc(size_x * size_y * 3, sizeof(unsigned char));
    if (im->rgbdata != NULL)
      return (im);
  }
  fprintf(stderr, "Unable to allocate memory for new image\n");
  return (NULL);
}

void imageOutput(struct image *im, const char *filename)
{
  // Writes out a .ppm file from the image data contained in 'im'.
  // Note that Windows typically doesn't know how to open .ppm
  // images. Use Gimp or any other seious image processing
  // software to display .ppm images.
  // Also, note that because of Windows file format management,
  // you may have to modify this file to get image output on
  // Windows machines to work properly.
  //
  // Assumes a 24 bit per pixel image stored as unsigned chars
  //

  FILE *f;

  if (im != NULL)
    if (im->rgbdata != NULL)
    {
      f = fopen(filename, "wb+");
      if (f == NULL)
      {
        fprintf(stderr, "Unable to open file %s for output! No image written\n", filename);
        return;
      }
      fprintf(f, "P6\n");
      fprintf(f, "# Output from RayTracer.c\n");
      fprintf(f, "%d %d\n", im->sx, im->sy);
      fprintf(f, "255\n");
      fwrite((unsigned char *)im->rgbdata, im->sx * im->sy * 3 * sizeof(unsigned char), 1, f);
      fclose(f);
      return;
    }
  fprintf(stderr, "imageOutput(): Specified image is empty. Nothing output\n");
}

void deleteImage(struct image *im)
{
  // De-allocates memory reserved for the image stored in 'im'
  if (im != NULL)
  {
    if (im->rgbdata != NULL)
      free(im->rgbdata);
    free(im);
  }
}

void cleanup(struct object3D *o_list, struct pointLS *l_list, struct textureNode *t_list)
{
  // De-allocates memory reserved for the object list and the point light source
  // list. Note that *YOU* must de-allocate any memory reserved for images
  // rendered by the raytracer.
  struct object3D *p, *q;
  struct pointLS *r, *s;
  struct textureNode *t, *u;

  p = o_list; // De-allocate all memory from objects in the list
  while (p != NULL)
  {
    q = p->next;
    if (p->photonMap != NULL) // If object is photon mapped, free photon map memory
    {
      if (p->photonMap->rgbdata != NULL)
        free(p->photonMap->rgbdata);
      free(p->photonMap);
    }
    free(p);
    p = q;
  }

  r = l_list; // Delete light source list
  while (r != NULL)
  {
    s = r->next;
    free(r);
    r = s;
  }

  t = t_list; // Delete texture Images
  while (t != NULL)
  {
    u = t->next;
    if (t->im->rgbdata != NULL)
      free(t->im->rgbdata);
    free(t->im);
    free(t);
    t = u;
  }
}

double distance(struct point3D p1, struct point3D p2)
{
  return sqrt(pow((p1.px - p2.px), 2) + pow((p1.py - p2.py), 2) + pow((p1.pz - p2.pz), 2));
}
