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

#include "utils.h"	// <-- This includes RayTracer.h

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct textureNode *texture_list;
int MAX_DEPTH;

void buildScene(void)
{
#include "buildscene.c"		// <-- Import the scene definition! 
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

 struct colourRGB tmp_col;	// Accumulator for colour components
 double R,G,B;			// Colour for the object in R G and B

 // This will hold the colour as we process all the components of
 // the Phong illumination model
 tmp_col.R=0;
 tmp_col.G=0;
 tmp_col.B=0;

 if (obj->texImg==NULL)		// Not textured, use object colour
 {
  R=obj->col.R;
  G=obj->col.G;
  B=obj->col.B;
 }
 else
 {
  // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
  // for the object. Note that we will use textures also for Photon Mapping.
  obj->textureMap(obj->texImg,a,b,&R,&G,&B);
 }

 //////////////////////////////////////////////////////////////
 // TO DO: Implement this function. Refer to the notes for
 // details about the shading model.
 //////////////////////////////////////////////////////////////

 // Be sure to update 'col' with the final colour computed here!
// shadow ray
 struct point3D mirrorReflectionLS;
  struct point3D mirrorReflectionRay;
  struct point3D vectorLS;
  struct ray3D rayLS;
  struct point3D toCamera;
  struct pointLS *currentLightSource = light_list;  
  struct ray3D reflectedRay;
  struct ray3D cameraRay;
  double lambda;


  // Compute the c vector going to camera which is the opposite direction of ray
  toCamera.px = -(ray->d.px);
  toCamera.py = -(ray->d.py);
  toCamera.pz = -(ray->d.pz);
  toCamera.pw = 1;

  normalize(&toCamera);
  initRay(&cameraRay, p, &toCamera);

  // Make iterator variable to loop through all the lightsources
  struct pointLS *lightsource_iterator = light_list;

  // Loop through all the lightsources
  while(lightsource_iterator != NULL){

    // Need to make the s vector
    vectorLS = lightsource_iterator->p0;
    subVectors(p, &vectorLS);
    normalize(&vectorLS);
    // From the newly made light source vector make the light source ray
    initRay(&rayLS, p, &vectorLS);
    findFirstHit(&rayLS, &lambda, obj, &obj, p, n, &a, &b);
    // Compute the ambient term ra*Ia from Phong
    tmp_col.R += obj->alb.ra * R;
    tmp_col.G += obj->alb.ra * G;
    tmp_col.B += obj->alb.ra * B;

    if(0 < lambda && lambda < 1){
      //If lambda is between 0 and 1 we only need the
      //ambient term

      //Compute the ambient term ra*Ia from Phong
      tmp_col.R += 1; // obj->alb.ra * R;
      tmp_col.G += 1; //obj->alb.ra * G;
      tmp_col.B += 1;// obj->alb.ra * B;

      //Also need to accomodate the previous ray traced lighting rg*Ispec
      tmp_col.R += obj->alb.rg*col->R;
      tmp_col.G += obj->alb.rg*col->G;
      tmp_col.B += obj->alb.rg*col->B;

    }
    else{
      // If lambda is greater then we compute and
      // use the entire Phong model

      // Compute the ambient term ra*Ia from Phong
      tmp_col.R += obj->alb.ra * R;
      tmp_col.G += obj->alb.ra * G;
      tmp_col.B += obj->alb.ra * B;

      // Compute the dot product for the surface normal and the lightsource
      double norm_ls_dot = dot(n, &rayLS.d);

      // Compute diffuse term rd*Id*max(0, n*s)
      tmp_col.R += obj->alb.rd * R * lightsource_iterator->col.R * max(0, norm_ls_dot);
      tmp_col.G += obj->alb.rd * G * lightsource_iterator->col.G * max(0, norm_ls_dot);
      tmp_col.B += obj->alb.rd * B * lightsource_iterator->col.B * max(0, norm_ls_dot);


      // Compute the perfect mirror reflection of the lightsource
      // m = 2(s*n)n - s
      mirrorReflectionLS.px = 2*(dot(n, &rayLS.d))*n->px - rayLS.d.px;
      mirrorReflectionLS.py = 2*(dot(n, &rayLS.d))*n->py - rayLS.d.py;
      mirrorReflectionLS.pz = 2*(dot(n, &rayLS.d))*n->pz - rayLS.d.pz;
      normalize(&mirrorReflectionLS);

      // Compute dot product for the camera and the mirror reflection ray
      double cam_mirror_dot = dot(&mirrorReflectionLS, &toCamera);

      // Compute the specular term rs*Is*max(0, c*m)^alpha (alpha is shinyness here)
      tmp_col.R += obj->alb.rs * R * lightsource_iterator->col.R * pow(max(0, cam_mirror_dot),obj->shinyness);
      tmp_col.G += obj->alb.rs * G * lightsource_iterator->col.R * pow(max(0, cam_mirror_dot),obj->shinyness);
      tmp_col.B += obj->alb.rs * B * lightsource_iterator->col.R * pow(max(0, cam_mirror_dot),obj->shinyness);

      // Also need to accomodate the previous ray traced lighting rg*Ispec
      tmp_col.R += obj->alb.rg*col->R;
      tmp_col.G += obj->alb.rg*col->G;
      tmp_col.B += obj->alb.rg*col->B;
    }
    //Compute the global components
    if(depth < MAX_DEPTH){


      //Check if object has specular components
      if(obj->alb.rs != 0){

        //Calculate the mirror direction according to the formula
        //ms = -d + 2(n*d)n
        double norm_ray_dot = dot(n, &ray->d);
        struct colourRGB specularCol;

        mirrorReflectionRay.px = -ray->d.px + 2*(norm_ray_dot)*n->px;
        mirrorReflectionRay.py = -ray->d.py + 2*(norm_ray_dot)*n->py;
        mirrorReflectionRay.pz = -ray->d.pz + 2*(norm_ray_dot)*n->pz;

        // Make the reflected ray and trace this one
        initRay(&reflectedRay, p, &mirrorReflectionRay);
        rayTrace(&reflectedRay, depth + 1, &specularCol, obj);

        // Updating the Ispec term
        tmp_col.R += obj->alb.rg * specularCol.R;
        tmp_col.G += obj->alb.rg * specularCol.G;
        tmp_col.B += obj->alb.rg * specularCol.B;

      }
      // // Check if the object is refractive
      // if (obj->alpha < 1){


      //   // Call ray trace
      //   // rayTrace()
      //   break;

      // }


   }
    lightsource_iterator = lightsource_iterator->next;
  }

 // Be sure to update 'col' with the final colour computed here!
 
 // When setting the components cap it at 1 if they are above 1

 col->R = (tmp_col.R > 1) ? 1 : tmp_col.R;
 col->G = (tmp_col.G > 1) ? 1 : tmp_col.G;
 col->B = (tmp_col.B > 1) ? 1 : tmp_col.B;

 return;
 /* struct ray3D shadow_ray;
 initRay(&shadow_ray, &ray->p0, &ray->d);
 // ray from p to lightsource 
 struct ray3D LS_ray;
 // direction vector to ls (s)
 struct point3D s;
 s.px = light_list->p0.px - p->px;
 s.py = light_list->p0.py - p->py;
 s.pz = light_list->p0.pz - p->pz;
 initRay(&LS_ray, p, &s);

 // test for intersection btwn ray and other objects in the scene
 double lambda = -1;
 struct object3D *hitobj;
 struct point3D p_int, n_int;
 double a_int;
 double b_int;

 findFirstHit(&LS_ray,&lambda,obj, &hitobj, &p_int, &n_int, &a_int, &b_int);
 if(lambda > 0 && lambda < 1){
  // ambient term
    tmp_col.R += obj->alb.ra * R + obj->alb.rg * R;
    tmp_col.G += obj->alb.ra * G + obj->alb.rg * G;
    tmp_col.B += obj->alb.ra * B + obj->alb.rg * B;
 }

 else{
  //phong model
  // ambient
    tmp_col.R += obj->alb.ra * R;
    tmp_col.G += obj->alb.ra * G;
    tmp_col.B += obj->alb.ra * B;

  //diffuse
    double ns = dot(n, &s);
    tmp_col.R += obj->alb.rd * R*max(0, ns);
    tmp_col.G += obj->alb.rd * G*max(0, ns);
    tmp_col.B += obj->alb.rd * B*max(0, ns);

  //spec
    struct point3D c;
    c.px = ray->d.px*-1;
    c.py = ray->d.py*-1;
    c.pz = ray->d.pz*-1;
    c.pw = 1;
    normalize(&c);
    struct point3D m;
    m.px = 2*dot(n, &s)*n->px - s.px;
    m.py = 2*dot(n, &s)*n->py - s.py;
    m.pz = 2*dot(n, &s)*n->pz - s.pz;  
 }
 return; */

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
 
 //Setting Variables
 //Set lambda to be -1 for no hits default
 *lambda = -1;
 struct point3D intersection;
 struct point3D norm;
 double temp_lambda;

 //Copy current object
 struct object3D *obj_clone = object_list;

 //While we have an object in object list
 while(obj_clone != NULL){
  //If we're dealing with a new and diff object that our current
  if(obj_clone != Os){ 
   obj_clone->intersect(obj_clone, ray, &temp_lambda, &intersection, &norm, a, b);

   //if temp lambda was set through the intersection
   if (temp_lambda != -1){ 
    //if current lambda has not been set or temp lambda is better
     if ((*lambda == -1) || (temp_lambda < *lambda)){ 
       //Set values
       *p = intersection;
       *n = norm;
       *lambda = temp_lambda;
       *obj = obj_clone;
      }
   }

   //Go to next object in linked list
   obj_clone = obj_clone->next;

  }
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
 
 double lambda;		// Lambda at intersection
 double a,b;		// Texture coordinates
 struct object3D *obj;	// Pointer to object at intersection
 struct point3D p;	// Intersection point
 struct point3D n;	// Normal at intersection
 struct colourRGB I;	// Colour returned by shading function

 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
 {
   col->R=-1;
   col->G=-1;
   col->B=-1;
   return;
 }

 //First thing when ray tracing is to find the first hit
 findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);

 //If we hit an object
 if(lambda > 0  && obj != Os){
   rtShade(obj, &p, &n, ray, depth, a,b, &I);
 }
 else{//Otherwise our ray goes into the distance
   //set all intensities to be 0
   I.R = 0;
   I.G = 0;
   I.B = 0;
 }

 //Set colours
 col->R = I.R;
 col->G = I.G;
 col->B = I.B;

}

int main(int argc, char *argv[])
{
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;	// Will hold the raytraced image
 struct view *cam;	// Camera and view for this scene
 int sx;		// Size of the raytraced image
 int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];	// Name of the output file for the raytraced .ppm image
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct ray3D ray;		// Structure to keep the ray from e to a pixel
 struct colourRGB col;		// Return colour for raytraced pixels
 struct colourRGB background;   // Background colour
 int i,j;			// Counters for pixel coordinates
 unsigned char *rgbIm;

 if (argc<5)
 {
  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 light_list=NULL;
 texture_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 if (!im)
 {
  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
  exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 ///////////////////////////////////////////////////
 // TO DO: You will need to implement several of the
 //        functions below. For Assignment 2, you can use
 //        the simple scene already provided. But
 //        for Assignment 3 you need to create your own
 //        *interesting* scene.
 ///////////////////////////////////////////////////
 buildScene();		// Create a scene. This defines all the
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
 e.px=0;
 e.py=0;
 e.pz=-1;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin.
 g.px=0-e.px;
 g.py=0-e.py;
 g.pz=0-e.pz;
 g.pw=1;
 // In this case, the camera is looking along the world Z axis, so
 // vector w should end up being [0, 0, -1]

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=1;

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -1, -2, 2, 4);

 if(cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list,light_list, texture_list);
  deleteImage(im);
  exit(0);
 }

 // Set up background colour here
 background.R=0;
 background.G=0;
 background.B=0;

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
 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix:\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 fprintf(stderr,"Rendering row: ");
 for (j=0;j<sx;j++)		// For each of the pixels in the image
 {
  fprintf(stderr,"%d/%d, ",j,sx);
  for (i=0;i<sx;i++)
  {
    ///////////////////////////////////////////////////////////////////
    // TO DO - complete the code that should be in this loop to do the
    //         raytracing!
    ///////////////////////////////////////////////////////////////////
   //Set up the pixels coordinates
   pc.pz = cam->wl + du *i;
   pc.py = cam->wt + dv *j;
   pc.pz = cam->f;
   pc.pw = 1;

   //Convert local pixel coordinate to world coordinate
   matVecMult(cam->C2W, &pc);

   //Assign direction vector for pixel 
   d.px = pc.px;
   d.py = pc.py;
   d.pz = pc.pz;
   d.pw = 1;

   //Subtract e and d to make d the direction based on pixel originally and the axis
   subVectors(&e, &d);
   normalize(&d);
   struct ray3D ray;
  
   //Setting up ray
   initRay(&ray,&pc, &d);
   
   //trace ray
   struct object3D obj;
   col = background;
   rayTrace(&ray, 0, &col,&obj);

  } // end for i
 } // end for j

 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list,light_list,texture_list);		// Object, light, and texture lists
 deleteImage(im);					// Rendered image
 free(cam);						// camera view
 exit(0);
}

