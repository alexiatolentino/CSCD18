/*
  CSC D18 - Assignment 1 - 2D light propagation

  This is the place where you will be doing most of your work for solving this
  assignment. Before you start working here, you shold have become familiar
  with the data structures and functions available to you from light2D.h, and
  with the way a scene is built in buildScene.c

  Go over each of the functions below, and implement the different components
  of the solution in the sections marked

  /************************
  / TO DO:
  ************************ /

  Do not add or modify code outside these sections.

  Details about what needs to be implemented are described in the comments, as
  well as in the assignment handout. You must read both carefully. 

  Starter by: F.J. Estrada, Aug. 2017
*/

/****************************************************************************
 * Uncomment the #define below to enable debug code, add whatever you need
 * to help you debug your program between #ifdef - #endif blocks
 * ************************************************************************/
#define __DEBUG_MODE

/*****************************************************************************
* COMPLETE THIS TEXT BOX:
*
* 1) Student Name: Reem Al Halabi
* 2) Student Name: Alexia Tolentino
*
* 1) Student number: 1006321927
* 2) Student number: 1006293967
* 
* 1) UtorID alhala16
* 2) UtorID tolent55
* 
* We hereby certify that the work contained here is our own
*
* ___Reem Al Halabi____             ___Alexia Tolentino______
* (sign with your name)            (sign with your name)
********************************************************************************/
// Helper functions here:
// This is a helper function that takes the minimum of 2 numbers
int min(double a, double b){
  return(a > b) ? b : a;
}

struct ray2D makeLightSourceRay(void)
{
 /*
   This function should return a light ray that has its origin at the light
   source, and whose direction depends on the type of light source.
   For point light sources (which emit light in all directions) the direction
    has to be chosen randomly and uniformly over a unit circle (i.e. any
    direction is equally likely)

   For a laser type light source, the direction of the ray is exactly the same
    as the direction of the lightsource.
    
   Set the colour of the ray to the same colour as the light source, and
    set the inside_outside flag to 0 (we assume light sources are 
    outside objects)

   In either case, the direction vector *must be unit length*.
*/
 
 /************************************************************************
 *  TO DO: Complete this function so that we can sample rays from the
 *         lightsource for the light propagation process.
 ************************************************************************/
 
 struct ray2D ray;
 //Set Ray's Origin
  ray.p.px=lightsource.l.p.px;		
  ray.p.py=lightsource.l.p.py;
 
 //FIRST CASE: Lightsource is a Laser
 if (lightsource.light_type == 1){
    //Normalize Lightsource direction:
    normalize(&lightsource.l.d);	

    //Set Ray's direction
    ray.d.px=lightsource.l.d.px;			
    ray.d.py=lightsource.l.d.py;	

 }
 //SECOND CASE: Lighsource is a Point
 else{
    //Generate direction:
    struct point2D rand_d;
    
    //Generate random point on unit circle:
    double theta = (double)rand()/(double)(RAND_MAX/(2*PI));//Getting random angle in 2pi
    rand_d.px = ray.p.px + cos(theta); //attain x coordinate of angle
    rand_d.py = ray.p.py + sin(theta); //attain y coordinate of angle

    //Normalize direction:
    normalize(&rand_d);	

    //Set Ray's direction
    ray.d.px=rand_d.px;			
    ray.d.py=rand_d.py;			
 }

  ray.inside_out=lightsource.l.inside_out;		// Initially 0 since the ray starts outside an object
  ray.monochromatic=lightsource.l.monochromatic;		// Initially 0 since the ray is white (from lightsource)
  ray.R=lightsource.R;			// Ray colour in RGB must be the same as the lightsource
  ray.G=lightsource.G;
  ray.B=lightsource.B;
  
  return(ray);

 // This creates a dummy ray (which won't go anywhere since the direction
 // vector d=[0 0]!. But it's here so you see which data values you have to
 // provide values for given the light source position, and type.  
 // ** REPLACE THE CODE BELOW with code that provides a valid ray that
 //    is consistent with the lightsource.

}

void propagateRay(struct ray2D *ray, int depth)
{
 /*
   This function carries out the light propagation process. It is provided with access
   to a ray data structure, and must perform the following steps (in order!):
   
   - Check if maximum recursion depth has been reached (in which case, it just returns)
   - Find the *closest* intersection between the ray and objects in the scene. This
     means you have to check against the 4 walls, and any circles added in buildScene,
     determine which intersection is closest, and obtain the intersection point, the
     normal at the intersection, and the lambda value at which the intersection happens.
   - Renders the ray onto the image from its starting point all the way up to the 
     intersection point.
   - At the intersection, use the material properties to determine how the propagation
     process proceeds:
         * For mirror materials, compute the mirror reflection direction, create a ray
           along that direction whose origin is the intersection point, and propagate
           that ray
         * For scattering materials, choose a random direction within +- 90 degrees of
           the normal, create a ray with that direction with origin at the intersection
           point, and propagate that ray
         * For refracting materials you will need to propagate two rays - one in the
           mirror reflection direction (just like for reflecting materials), and
           another in the refraction direction. Propagate both rays one after the other!
           
   NOTE: You should only care about intersections for which lambda is POSITIVE (in front
         of the ray), and greater than zero (e.g. if the ray is reflected from some
         object, you do not care about the intersection with the object itself which will
         have a lambda of *very close* to zero)
    
   In every case, make sure the ray's direction vector has unit length. You will need to
   complete other functions as part of your work here.
*/
  
 /*********************************************************************************
  * TO DO: Complete this function to implement light propagation!
  ********************************************************************************/
 
 // Define your local variables here
 
 if (depth>=max_depth) return;	 	// Leave this be, it makes sure you don't
					// recurse forever
 

 // Step 1 - Find *closest* intersection with the 4 walls (the written part of A1
 //          should help you figure out how to do that.
 int intersected = 0; // set as Fasle
 int j = 0;
 struct point2D q;
 while (intersected == 0){
    struct wall2D current_wall = walls[j];
    // finding the lambda value
    // checking if the top and bottom walls intersect
    if(j%2 == 0){
      printf(" we made it here %f", j);
      // y coord of intersection point is y value of wall
      q.py = current_wall.w.p.py;
      double lambda = (current_wall.w.p.py - ray->p.py) / ray->d.py;
      printf(" this is lambda : %f", lambda);

      // If lambda is negative, we know the ray isn't heading towards that wall
      if(lambda < 0){
        j++;
      }
      // otherwise we intersect at point (qx,qy) if lambda is positive
      else{
        q.px = ray->p.px + (lambda* ray->d.px);
        // checking if ray intersects withing bounds
        if(abs(q.px) < abs(current_wall.w.p.px)){
          intersected = 1;
          printf(" we are intersected with wall:  %f, %f", current_wall.w.p.py, current_wall.w.p.px );
        }
        j++;
      }
    }
    else{
      // x coord of intersection point is x value of wall
      q.px = current_wall.w.p.px;
      double lambda = (current_wall.w.p.px - ray->p.px) / ray->d.px;
      printf(" this is lambda : %f", lambda);
      printf(" this is q.px : %f", q.px);
      printf(" this is ray->p.px : %f", ray->p.px);
      printf(" this is ray->d.px : %f", ray->d.px);
      printf(" this is current_wall.w.p.px : %f",current_wall.w.p.px);
      // If lambda is negative, then we know it won't intersect
      if(lambda < 0){
        j++;
      }
      // otherwise we intersect at point (qx,qy) if lambda is positive
      
      else{
        q.py = ray->p.py + (lambda* ray->d.py);
        // checking if ray intersects withing bounds
        if(abs(q.py) < abs(current_wall.w.p.py)){
          intersected = 1;
          printf(" we are intersected with wall:  %f, %f", current_wall.w.p.py, current_wall.w.p.px );
        }
        j++;
      }
    }

 }

 // How many walls can the ray intersect? how many walls can the ray intersect in the
 // forward direction?

 // Step 2 - Check for intersection against objects in the object array - you must
 //          complete the intersectRay() function, call it, and obtain the closest
 //          intersection (in the forward ray direction) with objects in the scene.
 //          Note that you must provide variables for intersectRay() to return
 //          the point of intersection, normal at intersection, lambda, material type,
 //          and refraction index for the closest object hit by the ray.
 struct point2D p;
 // Assigning large values for coordinates
 //.px = 100000;
 //p.py = 100000;
 struct point2D n;
 double lambda;
 int type;
 double r_idx;
 // closest object


 printf(" this is RAY : %f", ray->p.px);
 intersectRay(ray, &p, &n, &lambda, &type, &r_idx);
 // Step 3 - Check whether the closest intersection with objects is closer than the
 //          closest intersection with a wall. Choose whichever is closer.
 
 
 // Distance between our ray and and intersected line 
 double magnitude_wall = sqrt(pow((q.px - ray->p.px),2)+pow((q.py - ray->p.py),2));
 double magnitude_obj = sqrt(pow((p.px - ray->p.px),2)+pow((p.py - ray->p.py),2));
 printf("magnitude wall is : %f\n, magnitude object is : %f\n",magnitude_wall,magnitude_obj);
 if(magnitude_wall < magnitude_obj){
  p = q;
 }
 

 // Step 4 - Render the ray onto the image. Use renderRay(). Provide renderRay() with
 //          the origin of the ray, and the intersection point (it will then draw a
 //          ray from the origin to the intersection). You also need to provide the
 //          ray's colour.

 renderRay(&ray->p, &p, ray->R, ray->G, ray->B);

 // Step 5 - Decide how to handle the ray's bounce at the intersection. You will have
 //          to provide code for 3 cases:
 //          If material type = 0, you have a mirror-reflecting object. 
 //                                Create a ray in the mirror reflection direction,
 //                                with the same colour as the incoming ray, and
 //                                with origin at the intersection point.
 //                                Then call propagateRay() recursively to trace it.
 //          if material type = 1, you have a scattering surface. 
 //                                Choose a random direction within +- 90 degrees 
 //                                from the normal at the intersection. Create a
 //                                ray in this direction, with the same colour as
 //                                the incoming ray, and origin at the intersection,
 //                                then call propagateRay() recursively to trace it.
 //          if material type = 2, you have a refracting (transparent) material.
 // 				   Here you need to process two rays:
 //                                * First, determine how much of the incoming light is
 //                                  reflected and how much is transmitted, using 
 //				     Schlick's approximation:
 // 					 R0 = ((n1-n2)/(n1+n2))^2   
 // 					 R(theta)=R0+((1-R0)*(1-cos(theta))^5)
 //				     If the ray is travelling from air to the inside
 //                                  of an object, n1=1, n2=object's index of refraction.
 //                                  If the ray is travelling from inside an object
 //                                  back onto air, n1=object's index of refraction, n2=1
 //				     And 'theta' is the angle between the normal and the
 // 				     ray direction.
 //				     R(theta) gives the amount Rs of reflected light, 
 //				     1.0-R(theta) gives the amount Rt of transmitted light.
 //                                * Now, make a ray in the mirror-reflection direction
 //				     (same as for material type 0), with the same colour
 //				     as the incoming ray, but with intensity modulated
 //				     by Rs. (e.g. if the incoming's colour is R,G,B,
 //                                  the reflected ray's colour will be R*Rs, G*Rs, B*Rs)
 //				     trace this ray.
 //				   * Make a ray in the refracted-ray direction. The 
 //				     angle for the transmitted ray is given by Snell's law
 //				     n1*sin(theta1) = n2*sin(theta2). The colour of the
 //				     transmitted ray is the same as the incoming ray but
 //			             modulated by Rt. Trace this ray.
 //	That's it! you're done!
   
}

void intersectRay(struct ray2D *ray, struct point2D *p, struct point2D *n, double *lambda, int *type, double *r_idx)
{

 /*
  This function checks for intersection between the ray and any objects in the objects 
  array. The objects are circles, so we are in fact solving for the intersection
  between a ray and a circle.
  
  For a unit circle centered at the origin, we would have the equation
  
  x^2 + y^2 = 1
  
  Using vector notation, with C=[x y]', we get
  
  ||C||^2 = 1
  
  A point on the ray is given by p + lambda*d
  
  Substituting in the equation for the circle we have 
  
  (p + lambda*d)(p + lambda*d) - 1 = 0
  
  If we expand the product above (here the product of two vectors is a DOT product), 
  we can form a quadratic equation
  
  A*lambda^2 + B*lambda + C = 0
  
  Which as you know, has a very simple solution. 
  
  Your task is to 
  * Figure out A, B, and C, considering that your circles don't necessarily have r=1
  * Figure out how to deal with the fact that circles in the scene are likely
    *not* centered at the origin
    
  Then implement the code that will find the value(s) of lambda at the intersection(s).
  
  Note that you can have no intersections, 1 intersection, or 2 intersections
  
  This function *must* find the closest intersection (if any) and update the value
  of lambda, the intersection point p, the normal n at the intersection, 
  the corresponding object's material type (needed outside here to figure out how
  to handle the light's bouncing off this object), and the index of refraction for
  the object (needed if this is a transparent object). 
  
  You must make sure that n is a unit-length vector.
 */
 
 /**********************************************************************************
  * TO DO: Complete this function to find the closest intersection between the
  *        ray and any objects in the scene, as well as the values at the
  *        intersection that will be needed to determine how to bounce/refract the
  *	   ray.
  * *******************************************************************************/
 //int distances[MAX_OBJECTS];
 double distances;
 //distances = (double*)calloc(MAX_OBJECTS,sizeof(double));
 double min_dist = 10000;
 int i = 0;

 // looping through the list of objects
 for(i; i<MAX_OBJECTS; i++){
   if (objects[i].r<=0) break;
    // x and y coordinates for intersection point
    struct point2D q;
    // declare variables for the quadratic function
    double point_diff = sqrt(pow((ray->p.px - objects[i].c.px),2)+pow((ray->p.py - objects[i].c.py),2));
    //double x = ray->p.px;
    //double c = objects[i].c.px;
    struct point2D d = ray->d; //d which is direction of ray
    struct point2D point_mult_d; //point diff mult by d
    point_mult_d.px = point_diff*d.px;
    point_mult_d.py = point_diff*d.py;
    struct point2D xminusc; //line from center of circle to start of ray
    xminusc.px = 2*(ray->p.px - objects[i].c.px);
    xminusc.py = 2*(ray->p.py - objects[i].c.py);
    
    double mag_d = sqrt(pow((d.px),2)+pow((d.py),2));
    

    double r = objects[i].r;
    double lambda_min;

    //CHECK DIRECTIONS: If the direction of the object to starting point 
    //                  is oppsite our direction vector, they will intersect
    if(((xminusc.px < 0) == (-d.px < 0)) && ((xminusc.py < 0) == (-d.py < 0))){
      // DETERMINE NUMBER OF SOLUTIONS: using discriminant we see how many solutions we have
      double disc = (pow((dot(&xminusc,&d)),2))-(4*(pow(mag_d,2))*(pow((point_diff),2)-pow(r,2)));

      //Positive Discriminant => 2 solutions
      if(disc > 0){
        printf("WE HAVE TWO SOLUTIONS %f\n", disc);

        // Finding lambda using quadratic formula
        double lambda_1 = (-dot(&xminusc,&d) + sqrt(disc))/(2*pow(mag_d,2)); // for plus case
        //calculate the first intersection point
        q.px = ray->p.px + lambda_1*ray->d.px;
        q.py = ray->p.py + lambda_1*ray->d.py;
        // calculate the distance between ray source and intersection point
        double distance_1 = sqrt(pow((q.px - ray->p.px),2)+pow((q.py - ray->p.py),2));
        
        double lambda_2 = (-dot(&xminusc,&d) - sqrt(disc))/(2*pow(mag_d,2)); // for minus case
        //calculate the second intersection point
        q.px = ray->p.px + lambda_2*ray->d.px;
        q.py = ray->p.py + lambda_2*ray->d.py;
        // calculate the distance between ray source and intersection point
        double distance_2 = sqrt(pow((q.px - ray->p.px),2)+pow((q.py - ray->p.py),2));

        // Take the minimum of the two and append to array of 
        if(distance_1 < distance_2){
          distances = distance_1;
          *lambda = lambda_1;
          lambda_min = lambda_1;
        }
        else{
          distances = distance_2;
          *lambda = lambda_2;
          lambda_min = lambda_2;
        }

      } 
      
      //Zero Discriminant => 1 solution
      else if(disc = 0){
        printf("WE HAVE one SOLUTION %f\n", disc);
        // Finding lambda using quadratic formula
        lambda_min = (-dot(&xminusc,&d) + sqrt(disc))/(2*pow(mag_d,2));

        //calculate the first intersection point
        q.px = ray->p.px + lambda_min*ray->d.px;
        q.py = ray->p.py + lambda_min*ray->d.py;
        // Take the minimum of the two and append to array of distances
        distances = sqrt(pow((q.px - ray->p.px),2)+pow((q.py - ray->p.py),2));
      } 
      
      //Negative Discriminant => No solutions 
      else{
        printf("WE HAVE no SOLUTIONS %f\n", disc);
        distances = min_dist;
      }
    }//DIRECTION CHECK: They don't go opposite ways, they won't intersect
    else{
      // we have no solutions
        printf("Direction diff! No intersection!");
        distances = min_dist;
    }


    //GETTING MINIMUM DISTANCE!
    if(distances < min_dist){
      //printf("WE FOUND MIN !!!! %f", distances);
      // updating the min value
      min_dist = distances;
      // update parameters
      *lambda = lambda_min;
      p->px = q.px;
      p->py = q.py;
      *type = objects[i].material_type;
      *r_idx = objects[i].r_idx;
      
      // finding the unit normal
      struct ray2D tangent;
      tangent.p = q;
      tangent.d.px = -2*PI*r*(sin(2*PI*lambda_min));
      tangent.d.py = 2*PI*r*(cos(2*PI*lambda_min));
      normalize(&tangent.d);
      *n = tangent.d;
    }
 }
}
