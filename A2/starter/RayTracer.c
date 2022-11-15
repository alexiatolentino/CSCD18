/*
  CSC D18 - RayTracer code.

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO" - remember to check what
  functionality is actually needed for the corresponding
  assignment!

  Last updated: Aug. 2017   - F.J.E.
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
 ********************************************************************************/

#include "utils.h" // <-- This includes RayTracer.h

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct textureNode *texture_list;
int MAX_DEPTH;

void buildScene(void)
{
#include "buildscene.c" // <-- Import the scene definition!
}

/*Helper function for phone model
  E = raIa + rdId(max(0, ndots)) + rsIs(max(0,cdotm)) + rgIspec
  This function sets...
  the ambient component: raIa
  the diffuse component: rdId(max(0, ndots))
  and
  the specular component: rsIs(max(0,cdotm))

  which is used by the regular rtShade
*/
void phong(struct object3D *obj, struct point3D *p, struct point3D *n, 
          struct ray3D *ray_ls, struct ray3D *ray, double *amb, 
          struct colourRGB *diff, struct colourRGB *spec, double k)
{
  // Setting Ambient component
  *amb = obj->alb.ra;

  normalize(&ray_ls->d);

  // Set m  = 2(dot(n, &ray_ls->d)) * (n) - ray_ls->d
  struct point3D m_dir;
  m_dir.px = 2 * (dot(n, &ray_ls->d)) * n->px - ray_ls->d.px;
  m_dir.py = 2 * (dot(n, &ray_ls->d)) * n->py - ray_ls->d.py;
  m_dir.pz = 2 * (dot(n, &ray_ls->d)) * n->pz - ray_ls->d.pz;
  m_dir.pw = 1;
  normalize(&m_dir);

  // Set c = -d
  struct point3D c_dir;
  c_dir = ray->d;
  normalize(&c_dir);
  c_dir.px = -1 * c_dir.px;
  c_dir.py = -1 * c_dir.py;
  c_dir.pz = -1 * c_dir.pz;
  c_dir.pw = 1;
  normalize(&c_dir);

  if (obj->frontAndBack == 0)
  {
    // Setting Diffuse component
    diff->R = obj->alb.rd * light_list->col.R * max(0, dot(n, &ray_ls->d)) * k;
    diff->G = obj->alb.rd * light_list->col.G * max(0, dot(n, &ray_ls->d)) * k;
    diff->B = obj->alb.rd * light_list->col.B * max(0, dot(n, &ray_ls->d)) * k;

    // Setting Specular component
    spec->R = obj->alb.rs * light_list->col.R * pow(max(0, dot(&c_dir, &m_dir)), obj->shinyness) * k;
    spec->G = obj->alb.rs * light_list->col.G * pow(max(0, dot(&c_dir, &m_dir)), obj->shinyness) * k;
    spec->B = obj->alb.rs * light_list->col.B * pow(max(0, dot(&c_dir, &m_dir)), obj->shinyness) * k;
  }
  else
  {
    // Setting diffuse component
    diff->R = obj->alb.rd * light_list->col.R * max(0, abs(dot(n, &ray_ls->d))) * k;
    diff->G = obj->alb.rd * light_list->col.G * max(0, abs(dot(n, &ray_ls->d))) * k;
    diff->B = obj->alb.rd * light_list->col.B * max(0, abs(dot(n, &ray_ls->d))) * k;

    // Setting Specular component
    spec->R = obj->alb.rs * light_list->col.R * pow(max(0, abs(dot(&c_dir, &m_dir))), obj->shinyness) * k;
    spec->G = obj->alb.rs * light_list->col.G * pow(max(0, abs(dot(&c_dir, &m_dir))), obj->shinyness) * k;
    spec->B = obj->alb.rs * light_list->col.B * pow(max(0, abs(dot(&c_dir, &m_dir))), obj->shinyness) * k;
  }
}


void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
  // This function implements the shading model as described in lecture. It takes
  // - A pointer to the first object intersected by the ray (to get the colour properties)
  // - The coordinates of the intersection point (in world coordinates)
  // - The normal at the point
  // - The ray (needed to determine the reflection direction to use for the global component, as well as for
  //   the Phong specular component)
  // - The current racursion depth
  // - The (a,b) texture coordinates (meaningless unless texture is enabled)
  //
  // Returns:
  // - The colour for this ray (using the col pointer)
  //

  struct colourRGB tmp_col; // Accumulator for colour components
  double R, G, B;           // Colour for the object in R G and B

  // This will hold the colour as we process all the components of
  // the Phong illumination model
  tmp_col.R = 0;
  tmp_col.G = 0;
  tmp_col.B = 0;

  if (obj->texImg == NULL) // Not textured, use object colour
  {
    R = obj->col.R;
    G = obj->col.G;
    B = obj->col.B;
  }
  else
  {
    // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
    // for the object. Note that we will use textures also for Photon Mapping.
    obj->textureMap(obj->texImg, a, b, &R, &G, &B);
  }

  //////////////////////////////////////////////////////////////
  // TO DO: Implement this function. Refer to the notes for
  // details about the shading model.
  //////////////////////////////////////////////////////////////

  // Setting Light Source
  struct pointLS *lightsource = light_list;

  // Be sure to update 'col' with the final colour computed here!
  // Initialize ray from intersection to light source
  struct ray3D ray_ls;
  

  // Test for intersection btwn ray and other objects in the scene
  double lambda = -1;
  struct object3D *hitobj;
  struct point3D p_int, n_int;
  double a_int, b_int;

  // Setting local and global components
  struct colourRGB local, global;
  double ambient, distance;
  struct colourRGB diffuse, specular;

  while(lightsource != NULL){
    struct point3D ls_dir = lightsource->p0;
    subVectors(p, &ls_dir);

    distance = length(&ls_dir);
    normalize(&ls_dir);

    initRay(&ray_ls, p, &ls_dir);
    ray_ls.inside = 0;

    findFirstHit(&ray_ls, &lambda, obj, &hitobj, &p_int, &n_int, &a_int, &b_int);

    lambda = lambda/distance;

    if (lambda > 0 && lambda < 1){
      // Ambient term
      local.R = obj->alb.ra;
      local.G = obj->alb.ra;
      local.B = obj->alb.ra;
    }
    else{
      // Use Phong model to determine local components
      phong(obj, p, n, &ray_ls, ray, &ambient, &diffuse, &specular, 1);
      local.R = R * (ambient + diffuse.R) + specular.R;
      local.G = G * (ambient + diffuse.G) + specular.G;
      local.B = B * (ambient + diffuse.B) + specular.B;
    }
    lightsource = lightsource->next;
  }
  
  struct object3D *curr_obj= object_list;
  int K = 1000;
  int k = 0;
  while(curr_obj != NULL){
    if(curr_obj->isLightSource){
      k = 0;
      //sample K points
      for (int j = 0; j < K; j++) {
        lambda = -1;
        struct point3D shadow;
        (*curr_obj->randomPoint)(curr_obj, &shadow.px, &shadow.py, &shadow.pz);
        shadow.pw = 1;

        ray_ls.d.px = shadow.px - p->px;
        ray_ls.d.py = shadow.py - p->py;
        ray_ls.d.pz = shadow.pz - p->pz;
        ray_ls.d.pw = 1;
        //subVectors(p, &ray_ls.p0); // switching

        distance = length(&ray_ls.d);
        normalize(&ray_ls.d);
        ray_ls.inside = 0;

        findFirstHit(&ray_ls, &lambda, curr_obj, &hitobj, &p_int, &n_int, &a_int, &b_int);

        lambda = lambda/distance;

        if(lambda < 0 || lambda > 1 || hitobj->isLightSource == 1){
          k+=1;
        }
      } 
      // SETTING LOCAL COMPONENTS
      if (lambda > 0 && lambda < 1)
      {
        //obj->textureMap(obj->texImg, a, b, &ct.R, &ct.G, &ct.B);
        // Ambient term
        local.R = obj->alb.ra;
        local.G = obj->alb.ra;
        local.B = obj->alb.ra;
      }
      else
      {
        // Use Phong model to determine local components
        phong(obj, p, n, &ray_ls, ray, &ambient, &diffuse, &specular, k/K);
        local.R += R * (ambient + diffuse.R) + specular.R;
        local.G += G * (ambient + diffuse.G) + specular.G;
        local.B += B * (ambient + diffuse.B) + specular.B;
      }
    }
    curr_obj = curr_obj->next;
  }

  // SETTING GLOBAL COMPONENTS
  if (depth < MAX_DEPTH)
  {
    // If object has specular components
    if (obj->alb.rg != 0)
    {

      // Get mirror direction
      struct ray3D mirror;

      mirror.p0 = *p;
      mirror.d.px = ray->d.px + (-2 * dot(n, &ray->d)) * n->px;
      mirror.d.py = ray->d.py + (-2 * dot(n, &ray->d)) * n->py;
      mirror.d.pz = ray->d.pz + (-2 * dot(n, &ray->d)) * n->pz;
      mirror.d.pw = 1;

      // Update global initial value
      global.R = 0;
      global.G = 0;
      global.B = 0;

      // Tracing the new mirror ray
      rayTrace(&mirror, depth + 1, &global, obj);

      // Updating the Ispec term
      global.R = obj->alb.rg * global.R;
      global.G = obj->alb.rg * global.G;
      global.B = obj->alb.rg * global.B;
    }
  }

  // Setting limit to local and global components = 1
  col->R = (local.R + global.R <= 1) ? local.R + global.R : 1;
  col->G = (local.G + global.G <= 1) ? local.G + global.G : 1;
  col->B = (local.B + global.B <= 1) ? local.B + global.B : 1;

  return;
}


void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
  // Find the closest intersection between the ray and any objects in the scene.
  // Inputs:
  //   *ray    -  A pointer to the ray being traced
  //   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
  //              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
  //              projection
  // Outputs:
  //   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
  //   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
  //              this ray (this is required so you can do the shading)
  //   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
  //   *n      -  A pointer to a 3D point structure so you can return the normal at the intersection point
  //   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point

  // SETTING VARIABLES
  *lambda = -1;
  struct point3D intersection;
  struct point3D norm;
  double temp_lambda = *lambda;
  double smallest_lambda = 10000;

  // Copy current object
  struct object3D *obj_clone = object_list;

  // While we have an object in object list
  while (obj_clone != NULL)
  {
    // If we're dealing with a new and diff object that our current
    if (obj_clone != Os)
    {
      (*obj_clone->intersect)(obj_clone, ray, &temp_lambda, p, n, a, b);

      // If temp lambda was set through the intersection
      // If current lambda has not been set or temp lambda is better
      if ((temp_lambda > 0) && (temp_lambda < smallest_lambda))
      {
        // Set values
        intersection = *p;
        norm = *n;
        smallest_lambda = temp_lambda;
        *lambda = temp_lambda;
        *obj = obj_clone;
      }
    }
    // Go to next object in linked list
    obj_clone = obj_clone->next;
  }
}
void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
  // Trace one ray through the scene.
  //
  // Parameters:
  //   *ray   -  A pointer to the ray being traced
  //   depth  -  Current recursion depth for recursive raytracing
  //   *col   - Pointer to an RGB colour structure so you can return the object colour
  //            at the intersection point of this ray with the closest scene object.
  //   *Os    - 'Object source' is a pointer to the object from which the ray
  //            originates so you can discard self-intersections due to numerical
  //            errors. NULL for rays originating from the center of projection.

  // SETTING VARIABLES
  double lambda, a, b;
  struct object3D *obj;
  struct point3D p, n;
  struct colourRGB I;

  // If max depth reached set color to -1
  if (depth > MAX_DEPTH)
  {
    col->R = -1;
    col->G = -1;
    col->B = -1;
    return; // end function
  }

  // First thing when ray tracing is to find the first hit
  findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);

  // If we hit an object get shading
  if (lambda > 0 && obj != Os)
  {
    rtShade(obj, &p, &n, ray, depth, a, b, &I);
  }
  else
  { // Otherwise our ray goes into the distance
    // set all intensities to be 0
    I.R = 0;
    I.G = 0;
    I.B = 0;
  }

  // Set colours
  col->R = I.R;
  col->G = I.G;
  col->B = I.B;
}

int main(int argc, char *argv[])
{
  // Main function for the raytracer. Parses input parameters,
  // sets up the initial blank image, and calls the functions
  // that set up the scene and do the raytracing.
  struct image *im;       // Will hold the raytraced image
  struct view *cam;       // Camera and view for this scene
  int sx;                 // Size of the raytraced image
  int antialiasing;       // Flag to determine whether antialiaing is enabled or disabled
  char output_name[1024]; // Name of the output file for the raytraced .ppm image
  struct point3D e;       // Camera view parameters 'e', 'g', and 'up'
  struct point3D g;
  struct point3D up;
  double du, dv;               // Increase along u and v directions for pixel coordinates
  struct point3D pc, d;        // Point structures to keep the coordinates of a pixel and
                               // the direction or a ray
  struct ray3D ray;            // Structure to keep the ray from e to a pixel
  struct colourRGB col;        // Return colour for raytraced pixels
  struct colourRGB background; // Background colour
  int i, j;                    // Counters for pixel coordinates
  unsigned char *rgbIm;

  if (argc < 5)
  {
    fprintf(stderr, "RayTracer: Can not parse input parameters\n");
    fprintf(stderr, "USAGE: RayTracer size rec_depth antialias output_name\n");
    fprintf(stderr, "   size = Image size (both along x and y)\n");
    fprintf(stderr, "   rec_depth = Recursion depth\n");
    fprintf(stderr, "   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
    fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
    exit(0);
  }
  sx = atoi(argv[1]);
  MAX_DEPTH = atoi(argv[2]);
  if (atoi(argv[3]) == 0)
    antialiasing = 0;
  else
    antialiasing = 1;
  strcpy(&output_name[0], argv[4]);

  fprintf(stderr, "Rendering image at %d x %d\n", sx, sx);
  fprintf(stderr, "Recursion depth = %d\n", MAX_DEPTH);
  if (!antialiasing)
    fprintf(stderr, "Antialising is off\n");
  else
    fprintf(stderr, "Antialising is on\n");
  fprintf(stderr, "Output file name: %s\n", output_name);

  object_list = NULL;
  light_list = NULL;
  texture_list = NULL;

  // Allocate memory for the new image
  im = newImage(sx, sx);
  if (!im)
  {
    fprintf(stderr, "Unable to allocate memory for raytraced image\n");
    exit(0);
  }
  else
    rgbIm = (unsigned char *)im->rgbdata;

  ///////////////////////////////////////////////////
  // TO DO: You will need to implement several of the
  //        functions below. For Assignment 2, you can use
  //        the simple scene already provided. But
  //        for Assignment 3 you need to create your own
  //        *interesting* scene.
  ///////////////////////////////////////////////////
  buildScene(); // Create a scene. This defines all the
                // objects in the world of the raytracer

  //////////////////////////////////////////
  // TO DO: For Assignment 2 you can use the setup
  //        already provided here. For Assignment 3
  //        you may want to move the camera
  //        and change the view parameters
  //        to suit your scene.
  //////////////////////////////////////////

  // Mind the homogeneous coordinate w of all vectors below. DO NOT
  // forget to set it to 1, or you'll get junk out of the
  // geometric transformations later on.

  // Camera center is at (0,0,-1)
  e.px = 0;
  e.py = 0;
  e.pz = -1;
  e.pw = 1;

  // To define the gaze vector, we choose a point 'pc' in the scene that
  // the camera is looking at, and do the vector subtraction pc-e.
  // Here we set up the camera to be looking at the origin.
  g.px = 0 - e.px;
  g.py = 0 - e.py;
  g.pz = 0 - e.pz;
  g.pw = 1;
  // In this case, the camera is looking along the world Z axis, so
  // vector w should end up being [0, 0, -1]

  // Define the 'up' vector to be the Y axis
  up.px = 0;
  up.py = 1;
  up.pz = 0;
  up.pw = 1;

  // Set up view with given the above vectors, a 4x4 window,
  // and a focal length of -1 (why? where is the image plane?)
  // Note that the top-left corner of the window is at (-2, 2)
  // in camera coordinates.
  cam = setupView(&e, &g, &up, -1, -2, 2, 4);

  if (cam == NULL)
  {
    fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
    cleanup(object_list, light_list, texture_list);
    deleteImage(im);
    exit(0);
  }

  // Set up background colour here
  background.R = 0;
  background.G = 0;
  background.B = 0;

  // Do the raytracing
  //////////////////////////////////////////////////////
  // TO DO: You will need code here to do the raytracing
  //        for each pixel in the image. Refer to the
  //        lecture notes, in particular, to the
  //        raytracing pseudocode, for details on what
  //        to do here. Make sure you undersand the
  //        overall procedure of raytracing for a single
  //        pixel.
  //////////////////////////////////////////////////////
  du = cam->wsize / (sx - 1);  // du and dv. In the notes in terms of wl and wr, wt and wb,
  dv = -cam->wsize / (sx - 1); // here we use wl, wt, and wsize. du=dv since the image is
                               // and dv is negative since y increases downward in pixel
                               // coordinates and upward in camera coordinates.

  fprintf(stderr, "View parameters:\n");
  fprintf(stderr, "Left=%f, Top=%f, Width=%f, f=%f\n", cam->wl, cam->wt, cam->wsize, cam->f);
  fprintf(stderr, "Camera to world conversion matrix (make sure it makes sense!):\n");
  printmatrix(cam->C2W);
  fprintf(stderr, "World to camera conversion matrix:\n");
  printmatrix(cam->W2C);
  fprintf(stderr, "\n");

  fprintf(stderr, "Rendering row: ");
  for (j = 0; j < sx; j++) // For each of the pixels in the image
  {
    fprintf(stderr, "%d/%d, ", j, sx);
    for (i = 0; i < sx; i++)
    {
      ///////////////////////////////////////////////////////////////////
      // TO DO - complete the code that should be in this loop to do the
      //         raytracing!
      ///////////////////////////////////////////////////////////////////
      // SETTING VARIABLES
      struct point3D u, w, v;
      struct point3D *t;

      // Finding w vector = -g/|g|
      double norm_g = sqrt(pow(g.px, 2) + pow(g.py, 2) + pow(g.pz, 2));
      w.px = (-g.px) / norm_g;
      w.py = (-g.py) / norm_g;
      w.pz = (-g.pz) / norm_g;
      w.pw = 0;

      // Finding u vector = wxup/|wxup|
      t = cross(&up, &w);
      double norm_t = sqrt(pow(t->px, 2) + pow(t->py, 2) + pow(t->pz, 2));
      u.px = t->px / norm_t;
      u.py = t->py / norm_t;
      u.pz = t->pz / norm_t;
      u.pw = 0;

      // Finding v vector = uxw
      struct point3D *wcrossu;
      wcrossu = cross(&w, &u);
      v.px = wcrossu->px;
      v.py = wcrossu->py;
      v.pz = wcrossu->pz;
      v.pw = 0;

      // Setting up converstion matrix for cam->world coordinates
      double mcw[4][4];
      mcw[0][0] = u.px;
      mcw[1][0] = u.py;
      mcw[2][0] = u.pz;
      mcw[3][0] = 0;
      mcw[0][1] = v.px;
      mcw[1][1] = v.py;
      mcw[2][1] = v.pz;
      mcw[3][1] = 0;
      mcw[0][2] = w.px;
      mcw[1][2] = w.py;
      mcw[2][2] = w.pz;
      mcw[3][2] = 0;
      mcw[0][3] = 0;
      mcw[1][3] = 0;
      mcw[2][3] = 0;
      mcw[3][3] = 1;
      nullSetupView(cam, &e, &g, &up, u, v, w);


      //Setting dimension for antialiasing
      int dimension;
      if(antialiasing){
        dimension = 3; 
      }
      else{
        dimension = 1;
      }

      struct colourRGB fullcol;
      fullcol.R = 0; fullcol.G = 0; fullcol.B = 0;

      for(int k = 0; k < dimension; k++){
        for(int n = 0; n < dimension; n++){
          // Setting camera point coordinates
          struct point3D point_coord;
          point_coord.px = cam->wl + du * i + k*du/dimension;
          point_coord.py = cam->wt + dv * j + n*dv/dimension;
          point_coord.pz = cam->f; // changed to f bc we can now change camera pos
          point_coord.pw = 1;

          // Convert local point coordinate to world coordinate using conversion matrix
          matVecMult(mcw, &point_coord);
          addVectors(&e, &point_coord);
          point_coord.pw = 1;

          // Set direction vector for pixel ij
          struct point3D dir;
          dir.px = point_coord.px;
          dir.py = point_coord.py;
          dir.pz = point_coord.pz;
          dir.pw = 1;

          // Setting direction to be d-e
          subVectors(&e, &dir);
          dir.pw = 1;
          normalize(&dir);

          // Creating new ray
          struct ray3D new_ray;
          initRay(&new_ray, &point_coord, &dir);

          // Setting initial colours
          col.R = background.R;
          col.G = background.G;
          col.B = background.B;

          // Trace the new ray
          rayTrace(&new_ray, 1, &col, NULL);

          // Adjust full colour
          fullcol.R += col.R;
          fullcol.G += col.G;
          fullcol.B += col.B;

        }
      }

      // Final full color adjustment 
      fullcol.R = fullcol.R/pow(dimension,2);
      fullcol.G = fullcol.G/pow(dimension,2);
      fullcol.B = fullcol.B/pow(dimension,2);

      // Setting img colour
      rgbIm[3 * (j * im->sx + (im->sx - i)) + 0] = fullcol.R * 255;
      rgbIm[3 * (j * im->sx + (im->sx - i)) + 1] = fullcol.G * 255;
      rgbIm[3 * (j * im->sx + (im->sx - i)) + 2] = fullcol.B * 255;
    } // end for i
  }   // end for j

  fprintf(stderr, "\nDone!\n");

  // Output rendered image
  imageOutput(im, output_name);

  // Exit section. Clean up and return.
  cleanup(object_list, light_list, texture_list); // Object, light, and texture lists
  deleteImage(im);                                // Rendered image
  free(cam);                                      // camera view
  exit(0);
}
