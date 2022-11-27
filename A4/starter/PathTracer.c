/*
  CSC D18 - Path Tracer code.

  Derived from the ray tracer starter code. Most function 
  names are identical, though in practice the implementation
  should be much simpler!

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
* 1) UtorID: tolent55
* 2) UtorID alhala16
* 
* We hereby certify that the work contained here is our own
*
* __Alexia Tolentino__             __Reem Al Halabi____
* (sign with your name)            (sign with your name)
********************************************************************************/

#include "utils_path.h"			// <-- This includes PathTracer.h
//#define __USE_IS			// Use importance sampling for diffuse materials
//#define __USE_ES			// Use explicit light sampling
//#define __DEBUG			// <-- Use this to turn on/off debugging output

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct textureNode *texture_list;
unsigned long int NUM_RAYS;
int MAX_DEPTH;

#include "buildScene.c"			// Import scene definition

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b, int hit_type)
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

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////
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


int explicitLS(struct ray3D *ray, struct point3D *old_p, struct point3D *old_n, struct object3D *obj){
  struct object3D *curr_obj = object_list;
  struct point3D d, p, n;
  double x,y,z, lambda, a, b, add_light;
  struct object3D *object;
  struct ray3D *new_ray;
  int CEL=1;
  // Looping through objects
  while(curr_obj != NULL){
    //  Check for LS
    if(curr_obj->isLightSource){
      (curr_obj->randomPoint)(curr_obj,&x,&y,&z);
      d.px = x - old_p->px;
      d.py = y - old_p->py;
      d.pz = z - old_p->pz;
      normalize(&d);
      initRay(new_ray, &old_p, &d);

      findFirstHit(new_ray, &lambda, NULL, &object, &p, &n, &a, &b);
      // Check if object is a LS
      if(object != NULL && object->isLightSource && lambda > 0){
        add_light = 2 * PI * object->LSweight * -dot(&n, &d) * dot(old_n, &d)/(pow(lambda, 2));
        if(add_light > 1){
          add_light = 1;
        }
        ray->Ir += ray->R * obj->col.R * add_light;
        ray->Ig += ray->G * obj->col.G * add_light;
        ray->Ib += ray->B * obj->col.B * add_light;
        CEL = 0;
      }
      
    }
    curr_obj = curr_obj->next;
  }
  return CEL;
}

void refraction(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, struct colourRGB *refract, struct colourRGB *reflect, int depth, struct colourRGB *tmp_col, int CEL){
  // If the object is refractive  (has index of refraction)
  if (obj->r_index != 1){
    struct ray3D refracted_ray;
    refracted_ray.p0 = *p;
    double n1, n2;

    struct point3D normal;

    //NO EMBEDDED OBJECTS W REFRACTION
    if(dot(n, &ray->d)<0){
      n1=1;
      n2 = obj->r_index; 
      normal.px = n->px;
      normal.py = n->py;
      normal.pz = n->pz;
      normal.pw = 1;
    }
    else{ 
      n1 = obj->r_index;
      n2 = 1; 
      normal.px = -n->px;
      normal.py = -n->py;
      normal.pz = -n->pz;
      normal.pw = 1;
    }
    
    struct point3D neg_norm;
    neg_norm.px = -n->px;
    neg_norm.py = -n->py;
    neg_norm.pz = -n->pz;

    // c is dot(-n, b)
    double c;
    c = dot(&neg_norm, &ray->d);

    double r;
    r = n1/n2;

    // rb
    struct point3D rb;
    rb.px = r*ray->d.px;
    rb.py = r*ray->d.py;
    rb.pz = r*ray->d.pz;

    // dt = rb + (rc - sqrt(1-r^2(1-c^2)))n
    double inner_factor = 1 - (pow(r,2) * (1 - pow(c,2)));
  

    //Check for total internal reflection
    if (inner_factor < 0 && n1 > n2){
      refracted_ray.p0 = *p;

      refracted_ray.d.px = rb.px + (r*c - sqrtl(inner_factor)) * normal.px;
      refracted_ray.d.py = rb.py + (r*c - sqrtl(inner_factor)) * normal.py;
      refracted_ray.d.pz = rb.pz + (r*c - sqrtl(inner_factor)) * normal.pz;
        
      PathTrace(&refracted_ray, depth+1, refract, obj, CEL);
    }
    else {
      refract->R = 0; 
      refract->G = 0;
      refract->B = 0;
    }
  }
}


void PathTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os, int CEL)
{
  // Trace one light path through the scene.
  //
  // Parameters:
  //   *ray   -  A pointer to the ray being traced
  //   depth  -  Current recursion dointer to an RGB epth for recursive raytracing
  //   *col   - Pcolour structure so you can return the object colour
  //            at the intersection point of this ray with the closest scene object.
  //   *Os    - 'Object source' is a pointer to the object from which the ray
  //            originates so you can discard self-intersections due to numerical
  //            errors. NULL for rays originating from the center of projection.

  double lambda; // Lambda at intersection
  double a, b;   // Texture coordinates
  int hit_type = 0;
  struct object3D *obj;   // Pointer to object at intersection
  struct point3D p;       // Intersection point
  struct point3D n;       // Normal at intersection
  double R, G, B;         // Handy in case you need to keep track of some RGB colour value
  double dice;            // Handy to keep a random value
  struct ray3D *next_ray; // For the new ray to be used in recursive calls

  if (depth > MAX_DEPTH) // Max recursion depth reached. Return black (no light coming into pixel from this path).
  {
    col->R = ray->Ir; // These are accumulators, initialized at 0. Whenever we find a source of light these
    col->G = ray->Ig; // get incremented accordingly. At the end of the recursion, we return whatever light
    col->B = ray->Ib; // we accumulated into these three values.
    // hitobj needs to account for both objects and lightsources
    struct object3D *hitobj = (struct object3D *)calloc(1, sizeof(struct object3D));
    findFirstHit(ray, &lambda, obj, &hitobj, &p, &n, &a, &b, hit_type);
    // Check if it's in the viewplane/ we hit an object
    int isrefractive = -1;
    if (lambda > 0)
    {
      if (hit_type == 1)
      {
        // hit an object - random sample
        // Generate a random number
        dice = drand48();
        ray->R *= obj->col.R;
        ray->G *= obj->col.G;
        ray->B *= obj->col.B;

        // Check the diffuse 
        struct point3D new_dir;
        if (dice < obj->diffPct)
        {
          // BRDF 
          cosWeightedSample(&n, &new_dir);
          CEL = explicitLS(ray, &p, &n, obj);
          ray->R *= dot(&n, &new_dir);
          ray->G *= dot(&n, &new_dir);
          ray->B *= dot(&n, &new_dir);
        }
        else
        {
          // Reflection
          new_dir.px = ray->d.px - (2 * dot(&ray->d, &n) * (n.px));
          new_dir.py = ray->d.py - (2 * dot(&ray->d, &n) * (n.py));
          new_dir.pz = ray->d.pz - (2 * dot(&ray->d, &n) * (n.pz));

          
          // FINE TUNEEEEEEEEE BEFORE!!! MULTIPLIED BY exp(1) * drand48() * 
          struct point3D rand_p;
          rand_p.px = cos(2 * PI * drand48()) * obj->refl_sig;
          rand_p.py = cos(2 * PI * drand48()) * obj->refl_sig;
          rand_p.pz = cos(2 * PI * drand48()) * obj->refl_sig;

          // Alligning the Ray properly
          if (dot(&rand_p, &new_dir) >= 0){
            addVectors(&rand_p, &new_dir);
          }
          else{
            subVectors(&rand_p, &new_dir);
          }

          // Check for refraction
          if (dice < obj->tranPct / (obj->reflPct + obj->tranPct))
          {
            isrefractive = 1;
            refraction(&n, &new_dir, obj, ray, &obj->col, CEL);
          }
        }

        // double a = drand48() * PI;            // get random angle from 0 to 2pi
        // double b = (drand48() - 0.5) * PI;        // get random angle from -pi to pi
        
        normalize(&new_dir);
        struct ray3D next_ray;
        initRay(&next_ray, &p, &new_dir);

        next_ray.R = ray->R;
        next_ray.G = ray->G;
        next_ray.B = ray->B;
        next_ray.Ir = ray->Ir;
        next_ray.Ig = ray->Ig;
        next_ray.Ib = ray->Ib;

        // Path trace recursively
        if(isrefractive == 1){
          obj = NULL;
        }
        NUM_RAYS++;
        PathTrace(&next_ray, depth + 1, col, obj, CEL);
      }
    
      else{
        // hit a LS - BRDF
        if (obj->isLightSource)
        {
          // Used to accumalte light
          if (CEL)
          {
            ray->Ir += ray->R * obj->col.R;
            ray->Ig += ray->G * obj->col.G;
            ray->Ib += ray->B * obj->col.B;
          }
          col->R = ray->Ir > 1 ? 1 : ray->Ir;
          col->G = ray->Ig > 1 ? 1 : ray->Ig;
          col->B = ray->Ib > 1 ? 1 : ray->Ib;
        }
      }
    }
  }
  return;
}

int main(int argc, char *argv[])
{
 // Main function for the path tracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;		// Will hold the final image
 struct view *cam;		// Camera and view for this scene
 int sx;			// Size of the  image
 int num_samples;		// Number of samples to use per pixel
 char output_name[1024];	// Name of the output file for the .ppm image file
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct ray3D ray;		// Structure to keep the ray from e to a pixel
 struct colourRGB col;		// Return colour for pixels
 int i,j,k;			// Counters for pixel coordinates and samples
 double *rgbIm;			// Image is now double precision floating point since we
				// will be accumulating brightness differences with a 
				// wide dynamic range
 struct object3D *obj;		// Will need this to process lightsource weights
 double *wght;			// Holds weights for each pixel - to provide log response
 double pct,wt;
 
 time_t t1,t2;
 FILE *f;
				
 if (argc<5)
 {
  fprintf(stderr,"PathTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: PathTracer size rec_depth num_samples output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   num_samples = Number of samples per pixel\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 num_samples=atoi(argv[3]);
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 fprintf(stderr,"NUmber of samples = %d\n",num_samples);
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 texture_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 wght=(double *)calloc(sx*sx,sizeof(double));
 if (!im||!wght)
 {
  fprintf(stderr,"Unable to allocate memory for image\n");
  exit(0);
 }
 else rgbIm=(double *)im->rgbdata;
 for (i=0;i<sx*sx;i++) *(wght+i)=1.0;
 
 buildScene();		// Create a scene. 
 
 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll get junk out of the
 // geometric transformations later on.

 // Camera center
 e.px=0;
 e.py=0;
 e.pz=-15;
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
 cam=setupView(&e, &g, &up, -3, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list, texture_list);
  deleteImage(im);
  exit(0);
 }

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

 // Update light source weights - will give you weights for each light source that add up to 1
 obj=object_list;
 pct=0;
 while (obj!=NULL)
 {
  if (obj->isLightSource)
   pct+=obj->LSweight;
  obj=obj->next;
 }
 obj=object_list;
 while (obj!=NULL)
 {
  if (obj->isLightSource)
  {
   obj->LSweight/=pct;
  }
  obj=obj->next;
 }
 fprintf(stderr,"\n");

 NUM_RAYS=0;

 t1=time(NULL);

 fprintf(stderr,"Rendering pass... ");
 for (k=0; k<num_samples; k++)
 {
  fprintf(stderr,"%d/%d, ",k,num_samples);
#pragma omp parallel for schedule(dynamic,1) private(i,j,pc,wt,ray,col,d)
  for (j=0;j<sx;j++)		// For each of the pixels in the image
  {
   for (i=0;i<sx;i++)
   {
    // Random sample within the pixel's area
    pc.px=(cam->wl+((i+(drand48()-.5))*du));
    pc.py=(cam->wt+((j+(drand48()-.5))*dv));
    pc.pz=cam->f;
    pc.pw=1;

    // Convert image plane sample coordinates to world coordinates
    matVecMult(cam->C2W,&pc);

    // Now compute the ray direction
    memcpy(&d,&pc,sizeof(struct point3D));
    subVectors(&cam->e,&d);		// Direction is d=pc-e
    normalize(&d);

    // Create a ray and do the raytracing for this pixel.
    initRay(&ray, &pc,&d);

    wt=*(wght+i+(j*sx));
    PathTrace(&ray,1, &col,NULL,1);
    (*(rgbIm+((i+(j*sx))*3)+0))+=col.R*pow(2,-log(wt));
    (*(rgbIm+((i+(j*sx))*3)+1))+=col.G*pow(2,-log(wt));
    (*(rgbIm+((i+(j*sx))*3)+2))+=col.B*pow(2,-log(wt));
    wt+=col.R;
    wt+=col.G;
    wt+=col.B;
    *(wght+i+(j*sx))=wt;
   } // end for i
  } // end for j  
  if (k%25==0)  dataOutput(rgbIm,sx,&output_name[0]);  		// Update output image every 25 passes
 } // End for k 
 t2=time(NULL);

 fprintf(stderr,"\nDone!\n");

 dataOutput(rgbIm,sx,&output_name[0]);
 
 fprintf(stderr,"Total number of rays created: %ld\n",NUM_RAYS);
 fprintf(stderr,"Rays per second: %f\n",(double)NUM_RAYS/(double)difftime(t2,t1));

 // Exit section. Clean up and return.
 cleanup(object_list,texture_list);			// Object and texture lists
 deleteImage(im);					// Rendered image
 free(cam);						// camera view
 free(wght);
 exit(0);
}

