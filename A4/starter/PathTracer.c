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
* We hereby certify that the work add_lightained here is our own
*
* __Alexia Tolentino__             __Reem Al Halabi____
* (sign with your name)            (sign with your name)
********************************************************************************/

#include "utils_path.h"			// <-- This includes PathTracer.h
//#define __USE_IS			// Use importance sampling for diffuse materials
#define __USE_ES			// Use explicit light sampling
#define __DEBUG			// <-- Use this to turn on/off debugging output

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct textureNode *texture_list;
unsigned long int NUM_RAYS;
int MAX_DEPTH;

#include "buildScene.c"			// Import scene definition

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

  *lambda = -1;
  struct point3D intersection;
  struct point3D norm;
  double smallest_lambda;

  // Copy current object
  struct object3D *obj_clone = object_list;

  // While we have an object in object list
  while (obj_clone != NULL)
  {
    // If we're dealing with a new and diff object that our current
    if (obj_clone != Os)
    {
      obj_clone->intersect(obj_clone, ray, &smallest_lambda, &intersection, &norm, a, b);

      // If temp lambda was set through the intersection
      // If current lambda has not been set or temp lambda is better
      if((*lambda < 0 || smallest_lambda < *lambda) && smallest_lambda > 0){
        // Set values
        *p = intersection;
        *n = norm;
        *lambda = smallest_lambda;
        *obj = obj_clone;
      }
    }
    // Go to next object in linked list
    obj_clone = obj_clone->next;
  }
  return;
    
}

void PathTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os, int CEL, struct ray3D *last)
{
 // Trace one light path through the scene.
 //
 // Parameters:
 //   *ray   -  A pointer to the ray being traced
 //   depth  -  Current recursion depth for recursive raytracing
 //   *col   - Pointer to an RGB colour structure so you can return the object colour
 //            at the intersection point of this ray with the closest scene object.
 //   *Os    - 'Object source' is a pointer to the object from which the ray 
 //            originates so you can discard self-intersections due to numerical
 //            errors. NULL for rays originating from the center of projection. 

  double lambda;			// Lambda at intersection
  double a,b;			// Texture coordinates
  struct object3D *obj;		// Pointer to object at intersection
  struct point3D p;		// Intersection point
  struct point3D n;		// Normal at intersection
  double R,G,B;			// Handy in case you need to keep track of some RGB colour value
  double dice;			// Handy to keep a random value
  struct ray3D *next_ray = (struct ray3D *)calloc(1, sizeof(struct ray3D));	// For the new ray to be used in recursive calls

  dice = .5*drand48();
  if (depth>MAX_DEPTH || dice > maximum(ray->R,ray->G,ray->B))	// Max recursion depth reached. Return black (no light coming into pixel from this path).
  {
    col->R=ray->Ir;	// These are accumulators, initialized at 0. Whenever we find a source of light these
    col->G=ray->Ig;	// get incremented accordingly. At the end of the recursion, we return whatever light
    col->B=ray->Ib;	// we accumulated into these three values.
    memcpy(last,ray,sizeof(struct point3D));
    return;
  }

  // Pointer to hit object
  struct object3D *hitobj;
  findFirstHit(ray,&lambda,Os,&hitobj,&p,&n,&a,&b);
  
  if(lambda > 0){

    if(!hitobj->isLightSource){ // ray hits a surface
      dice = drand48();
      // get a random new direction
      ray->R *= hitobj->col.R;
      ray->G *= hitobj->col.G;
      ray->B *= hitobj->col.B;

      // If the hit object is not diffuse 
      int hit_type;
      if(dice > hitobj->diffPct){
        // Pick your ray to either reflect or refract depending on
        // your rays components, more likely to transmissive if tranpct is larger and vise versa
          hit_type = dice < hitobj->tranPct/(hitobj->reflPct + hitobj->tranPct) ? 2 : 0;
      }
      else{
        hit_type = 1;
      }
      struct point3D new_dir;
      // BRDF color for next ray
      if(hit_type == 1){ // diffuse
        cosWeightedSample(&n,&new_dir);
        CEL = explicitLS(ray, &p,&n,hitobj);
        // add to accumulator
        ray->R*=dot(&n,&new_dir);
        ray->G*=dot(&n,&new_dir);
        ray->B*=dot(&n,&new_dir);
      }
      //reflection
      else{
        // Perfect Reflection
        new_dir.px = ray->d.px - (2*dot(&ray->d,&n)*(n.px)); 
        new_dir.py = ray->d.py - (2*dot(&ray->d,&n)*(n.py)); 
        new_dir.pz = ray->d.pz - (2*dot(&ray->d,&n)*(n.pz)); 
        double ndotd = dot(&n, &new_dir);
       
        // Blur Reflection by random factor
        struct point3D rand_p;
        rand_p.px = exp(1)*drand48() * cos(2*PI*drand48())*hitobj->refl_sig;
        rand_p.py = exp(1)*drand48() * cos(2*PI*drand48())*hitobj->refl_sig;
        rand_p.pz = exp(1)*drand48() * cos(2*PI*drand48())*hitobj->refl_sig;

        // adjust new direction
        if (dot(&rand_p, &new_dir) >= 0)
        {
            addVectors(&rand_p, &new_dir);
        }
        else
        {
            subVectors(&rand_p, &new_dir);
        }
      }
      // refraction
      if(hit_type == 2){ 
        refraction(hitobj, ray, &n, &new_dir);
      }

      // Set up ray
      normalize(&new_dir);
      initRay(next_ray, &p,&new_dir);
      next_ray->R=ray->R;
      next_ray->G=ray->G;
      next_ray->B=ray->B;
      next_ray->Ir=ray->Ir;
      next_ray->Ig=ray->Ig;
      next_ray->Ib=ray->Ib;
      // recursively path trace
      if (hit_type == 2){
        hitobj = NULL;
      }
      NUM_RAYS++;
      PathTrace(next_ray,depth+1,col,hitobj,CEL,last);
      free(next_ray);
    }
    // ray hits lightsource
    else{
      if(CEL){
        ray->Ir += ray->R * hitobj->col.R;
        ray->Ig += ray->G * hitobj->col.G;
        ray->Ib += ray->B * hitobj->col.B;
      }
      col->R=ray->Ir > 1 ? 1 : ray->Ir;
      col->G=ray->Ig > 1 ? 1 : ray->Ig;
      col->B=ray->Ib > 1 ? 1 : ray->Ib;
    }
  }else{
    col->R=ray->Ir;
    col->G=ray->Ig;
    col->B=ray->Ib;
    return;
  }
}

int explicitLS(struct ray3D *ray, struct point3D *pt, struct point3D *norm,struct object3D *obj){
  struct object3D *curr_obj = object_list;
  struct point3D d;
  double x,y,z;
  double lambda,a,b;
  struct point3D p,n;
  struct object3D *object;
  double add_light; // add_lightribution of the light source
  struct ray3D *next_ray = (struct ray3D *)calloc(1, sizeof(struct ray3D));
  int CEL=1;
  // loop through each object in the object list
  
  while(curr_obj != NULL){
    // find lightsource to cast a ray to
    if(curr_obj->isLightSource){
      (curr_obj->randomPoint)(curr_obj,&x,&y,&z);
      d.px = x - pt->px;
      d.py = y - pt->py;
      d.pz = z - pt->pz;
      normalize(&d);
      initRay(next_ray, pt,&d);
      findFirstHit(next_ray,&lambda,NULL,&object,&p,&n,&a,&b);
      if(object != NULL){
        if(object->isLightSource && lambda > 0){
          add_light = (2*PI*object->LSweight * -dot(&n,&d)*dot(norm,&d))/(lambda*lambda);
          // If light is less than 1 then we can add the factor
          if (add_light < 1){
            ray->Ir += ray->R * obj->col.R * add_light;
            ray->Ig += ray->G * obj->col.G * add_light;
            ray->Ib += ray->B * obj->col.B * add_light;
          }else{
            ray->Ir += ray->R * obj->col.R;
            ray->Ig += ray->G * obj->col.G;
            ray->Ib += ray->B * obj->col.B;
          }
          CEL = 0;
        }
      }
      free(next_ray);
    }
    curr_obj = curr_obj->next;
  }
  return CEL;
}

void refraction(struct object3D *obj, struct ray3D *ray, struct point3D *n, struct point3D *dir){
  // If the object is refractive  (has index of refraction)
  if (obj->r_index != 1){
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
    c = -1*dot(&normal, &ray->d);

    double r;
	r = n1/n2;

    // rb
    struct point3D rb;
    rb.px = r*ray->d.px;
    rb.py = r*ray->d.py;
    rb.pz = r*ray->d.pz;

    // dt = rb + (rc - sqrtl(1-r^2(1-c^2)))n
    double inner_factor = 1 - (pow(r,2) * (1 - pow(c,2)));

    dir->px = rb.px + (r*c - sqrtl(inner_factor)) * normal.px;
    dir->py = rb.py + (r*c - sqrtl(inner_factor)) * normal.py;
    dir->pz = rb.pz + (r*c - sqrtl(inner_factor)) * normal.pz;
  }
}

// return the maximum double value of a,b,c
double maximum (double a,double b,double c)
{
   if(a>b)
   {
      if(a>c)
      {
         return a;
      }
   }
   else if(b>a) 
   {
      if(b>c)
      {
         return b;
      }
   }
   else if(c>b)
   {
   	if(c>a)
      {
      	return c;
      }
   }else{
   	return a;
   }
}

void transformNormal(struct point3D *n_orig, struct point3D *n_transformed, struct object3D *obj)
{
 // Computes the normal at an affinely transformed point given the original normal and the
 // object's inverse transformation. From the notes:
 // n_transformed=A^-T*n normalized.
  ///////////////////////////////////////////
 // TO DO: Complete this function
 ///////////////////////////////////////////
 // Setting transformed ray to be original ray
  *n_transformed = *n_orig;

  //Transform matrix based T
  n_transformed->px=obj->Tinv[0][0]*n_orig->px;
  n_transformed->px+=obj->Tinv[1][0]*n_orig->py;
  n_transformed->px+=obj->Tinv[2][0]*n_orig->pz;

  n_transformed->py=obj->Tinv[0][1]*n_orig->px;
  n_transformed->py+=obj->Tinv[1][1]*n_orig->py;
  n_transformed->py+=obj->Tinv[2][1]*n_orig->pz;

  n_transformed->pz=obj->Tinv[0][2]*n_orig->px;
  n_transformed->pz+=obj->Tinv[2][1]*n_orig->py;
  n_transformed->pz+=obj->Tinv[2][2]*n_orig->pz;
  
  // Normalize normal!
  normalize(n_transformed);
  n_transformed->pw = 1;
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
 double eps = 0.00001; // epsilon to hold a small number 
 int CRUNCHY = 1; // Set CRUNCHY to 0 to do normal pathrtacing
///////////////////////////////////////////////////////////////////
 // Depth Of Field:

 // Set CRUNCHY to 1 to activate Depth of Field Effect
 // int CRUNCHY = 1;

 // Parameters for the apeture size (ap_size)  and Z1 (x_1)
 double ap_size = 1;
 double x_1 = 800;
/////////////////////////////////////////////////////////////////////
 // Bi-directional Path Tracing:

 // Set CRUNCHY to 2 to activate Bi-directional Path Tracing
 // int CRUNCHY = 2;

 // Variables for ls
 struct ray3D ray_ls;
 struct colourRGB col_ls;
 int connected_count = 0;
//////////////////////////////////////////////////////////////////////

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
 // Original Rendering Pass
 fprintf(stderr,"Rendering pass... ");
 for (k=0; k<num_samples; k++)
 {
  fprintf(stderr,"%d/%d, ",k,num_samples);
#pragma omp parallel for schedule(dynamic,1) private(i,j,pc,wt,ray,col,d)
  for (j=0;j<sx;j++)		// For each of the pixels in the image
  {
   for (i=0;i<sx;i++)
   {

    // depth of field 
    if(CRUNCHY == 1){

      // Get x_0 using thin lens eqn:
      double x_0 = x_1 * (-cam->f) / (x_1 - (-cam->f));
      x_0 = -(x_0);

      // Random sample within the pixel's area
      pc.px=(cam->wl+((i+(drand48()-.5))*du));
      pc.py=(cam->wt+((j+(drand48()-.5))*dv));
      pc.pz=x_0;
      pc.pw=1;
      // Convert image plane sample coordinates to world coordinates
      matVecMult(cam->C2W,&pc);
      // Now compute the ray direction
      memcpy(&d,&pc,sizeof(struct point3D));
      subVectors(&cam->e,&d);		// Direction is d=pc-e
      normalize(&d);


      // Create a ray and do the raytracing for this pixel.
      struct ray3D *r_0 = (struct ray3D *)calloc(1, sizeof(struct ray3D));;
      initRay(r_0, &pc,&d);

      // compute p_d of r_0 intersect with plane at distance x_1 along optical axis
      struct point3D p_d;
      struct point3D p1;
      subVectors(&pc, &p1);

      // find lambda and get the intersecion using ray eqn
      double r_lambda = dot(&p1, &cam->w)/dot(&d, &cam->w);
      p_d.px = r_0->p0.px + r_lambda*r_0->d.px;
      p_d.py = r_0->p0.py + r_lambda*r_0->d.py;
      p_d.pz = r_0->p0.pz + r_lambda*r_0->d.pz;
      
      // randomly select point p_0 on apeture
      // use circle of specified diameter
      struct point3D p_0;
      double theta = drand48() * 2* PI;
      p_0.px = (cos(theta))*ap_size* drand48();
      p_0.py = (sin(theta))*ap_size* drand48();
      p_0.pz = e.pz;
      p_0.pw = 1;
      

      // cast r_i from p_0 towards p_d
      struct point3D p_ddir = p_d;
      subVectors(&p_0, &p_ddir);
      normalize(&p_ddir);
      struct ray3D *r_i = (struct ray3D *)calloc(1, sizeof(struct ray3D));
      initRay(r_i, &p_0, &p_ddir);

      // raytrace r_i and accumulate color onto pix i,j color
      struct ray3D *last = (struct ray3D *)calloc(1, sizeof(struct ray3D));
      struct point3D pt;
      pt.px = 0;
      pt.py = 0;
      pt.pz = 0;
      initRay(last, &pt, &pt);
      PathTrace(r_i, 1, &col, NULL, 1, last);
      free(r_i);
      free(last);
      // Changing the color
      wt=*(wght+i+(j*sx));
      (*(rgbIm+((i+(j*sx))*3)+0))+=(col.R/num_samples)*pow(2,-log(wt));
      (*(rgbIm+((i+(j*sx))*3)+1))+=(col.G/num_samples)*pow(2,-log(wt));
      (*(rgbIm+((i+(j*sx))*3)+2))+=(col.B/num_samples)*pow(2,-log(wt));
      wt+=col.R;
      wt+=col.G;
      wt+=col.B;
      *(wght+i+(j*sx))=wt;
    free(r_0);
    }
    else if(CRUNCHY == 2){
      // SETTING UP INFO FOR PIXEL
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


      // SETTING UP INFO FOR LS RAY
      // Random sample within LS area for start point
      struct point3D p_ls,d_ls,og_n,ls_norm;
      struct object3D *obj_clone = object_list;
      struct ray3D *connectedRay = (struct ray3D *)calloc(1, sizeof(struct ray3D));
      struct ray3D *last = (struct ray3D *)calloc(1, sizeof(struct ray3D));
      struct ray3D *last_ls = (struct ray3D *)calloc(1, sizeof(struct ray3D));
      struct object3D *hitobj;
      struct object3D *ls_obj;
      double dice = drand48();
      double lambda,a,b;
      struct point3D p,n;

      p_ls.px = 0;
      p_ls.py = 0;
      p_ls.pz = 0;

      // Check if connection between ray and ray_ls
      initRay(connectedRay, &p_ls, &p_ls);
      initRay(last, &p_ls, &p_ls);
      initRay(last_ls, &p_ls, &p_ls);

      //Iterate through possible LS
      while(obj_clone != NULL){
        double x,y,z;
        if(obj_clone->isLightSource){
          //if starting point is unassigned
          if(p_ls.px == 0 && p_ls.py == 0 && p_ls.pz == 0){
            (obj_clone->randomPoint)(obj_clone,&x,&y,&z);
            p_ls.px = x;
            p_ls.py = y;
            p_ls.pz = z;
            normalize(&p_ls);
            p_ls.pw = 1;
            ls_obj = obj_clone;
          }  
          else if(dice > 0.5 && (p_ls.px != 0 || p_ls.py != 0 || p_ls.pz != 0)){ 
            // by random chance we change the direction wrt another ls
            (obj_clone->randomPoint)(obj_clone,&x,&y,&z);
            p_ls.px = x;
            p_ls.py = y;
            p_ls.pz = z;
            normalize(&p_ls);
            p_ls.pw = 1;
            ls_obj = obj_clone;
          }
        }
        obj_clone = obj_clone->next;
      }

      
      // Random Sample Direction
      // Normal Vector
      double theta = 2*PI*drand48();
      og_n.px = cos(theta);
      og_n.py = sin(theta);
      og_n.pz = 1;
      og_n.pw = 1;

      //Compute normal
      transformNormal(&og_n, &d_ls, ls_obj);
      normalize(&d_ls);


      // Create a ray and do the raytracing for this pixel.
      initRay(&ray, &pc,&d);
      initRay(&ray_ls, &p_ls, &d_ls);

      // Path Trace Function
      wt=*(wght+i+(j*sx));
      PathTrace(&ray,1,&col,NULL,1,last);
      PathTrace(&ray_ls,1,&col_ls,NULL,1,last_ls);

      // Check if connection between ray and ray_ls
      initRay(connectedRay, &last->p0, &last->p0);
      subVectors(&last_ls->p0,&connectedRay->d);	//Direction is from end of ray and end of rayls paths

      normalize(&connectedRay->d);
      findFirstHit(connectedRay, &lambda, NULL,&hitobj,&p,&n,&a,&b);

      // If we did not hit ray_ls end point and intersected something before
      if(abs(p.px - last_ls->p0.px) < eps && abs(p.px - last_ls->p0.px) < eps && abs(p.px - last_ls->p0.px) < eps ){
        // Give it the old color it should've been
        (*(rgbIm+((i+(j*sx))*3)+0))+=col.R*pow(2,-log(wt));
        (*(rgbIm+((i+(j*sx))*3)+1))+=col.G*pow(2,-log(wt));
        (*(rgbIm+((i+(j*sx))*3)+2))+=col.B*pow(2,-log(wt));
        
        wt+=col.R;
        wt+=col.G;
        wt+=col.B;
        *(wght+i+(j*sx))=wt;
      }
      else{ // Adjust by probability of rays hitting eachother
        //Ray colour up until it's point:
        //fprintf(stderr, "WOO CONNECT");
        double alpha = 1/(pow(PI,2)); // Probability that it went in the exact angle of ray_ls
        connected_count+=1;
        //double alpha = 1;
        // accumuate light based on the probability of ray hitting that direction
        col.R += (col_ls.R * alpha); 
        col.G += (col_ls.G * alpha);
        col.B += (col_ls.B * alpha);

        (*(rgbIm+((i+(j*sx))*3)+0))+=col.R*pow(2,-log(wt));
        (*(rgbIm+((i+(j*sx))*3)+1))+=col.G*pow(2,-log(wt));
        (*(rgbIm+((i+(j*sx))*3)+2))+=col.B*pow(2,-log(wt));

        wt+=((1-alpha)*col.R + alpha*col_ls.R);
        wt+=((1-alpha)*col.G + alpha*col_ls.G);
        wt+=((1-alpha)*col.B + alpha*col_ls.B);
        *(wght+i+(j*sx))=wt;
      }
      free(connectedRay);
      free(last);
      free(last_ls);
    }
    else{
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

      struct ray3D *last = (struct ray3D *)calloc(1, sizeof(struct ray3D));
      struct point3D pt;
      pt.px = 0;
      pt.py = 0;
      pt.pz = 0;
      initRay(last, &pt, &pt);

      wt=*(wght+i+(j*sx));
      PathTrace(&ray,1, &col,NULL,1, last);
      (*(rgbIm+((i+(j*sx))*3)+0))+=col.R*pow(2,-log(wt));
      (*(rgbIm+((i+(j*sx))*3)+1))+=col.G*pow(2,-log(wt));
      (*(rgbIm+((i+(j*sx))*3)+2))+=col.B*pow(2,-log(wt));
      wt+=col.R;
      wt+=col.G;
      wt+=col.B;
      *(wght+i+(j*sx))=wt;

      free(last);
    }

   } // end for i
  } // end for j  
  if (k%25==0)  dataOutput(rgbIm,sx,&output_name[0]);  		// Update output image every 25 passes
 } // End for k 
 t2=time(NULL);

 // Output image 
 dataOutput(rgbIm,sx,&output_name[0]);

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
