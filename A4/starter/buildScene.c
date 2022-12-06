void buildScene(void)
{
 // Sets up all objets in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 //
 // NOTE: Light sources are now EXCLUSIVELY area light sources. They are
 //       defined by regular objects whose 'isLightSource' flag is set
 //       to 1. Therefore, you can create light sources with any shape
 //       and colour using the same object primitives and transforms
 //       you're using to set up the scene.
 //
 // To create hierarchical objects:
 //    You must keep track of transformations carried out by parent objects
 //    as you move through the hierarchy. Declare and manipulate your own
 //    transformation matrices (use the provided functions in utils.c to
 //    compound transformations on these matrices). When declaring a new
 //    object within the hierarchy
 //    - Initialize the object
 //    - Apply any object-level transforms to shape/rotate/resize/move
 //      the object using regular object transformation functions
 //    - Apply the transformations passed on from the parent object
 //      by pre-multiplying the matrix containing the parent's transforms
 //      with the object's own transformation matrix.
 //    - Compute and store the object's inverse transform as usual.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

 struct object3D *o;
 struct point3D p;

//  // Cornell box
//  o=newSphere(1.0,0.0,0.0,.75,.25,.25,.05,1.4);	// Left
//  Scale(o,500,500,500);
//  Translate(o,-510,0,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  o=newSphere(1.0,0.0,0.0,.25,.25,.75,.05,1.4);		// Right
//  Scale(o,500,500,500);
//  Translate(o,510,0,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  o=newSphere(1.0,0.0,0.0,.75,.75,.75,.05,1.4);		// Back
//  Scale(o,500,500,500);
//  Translate(o,0,0,515);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  o=newSphere(1.0,0.0,0.0,.75,.75,.75,.02,1.4);	// Bottom
//  Scale(o,500,500,500);
//  Translate(o,0,-510,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  o=newSphere(1.0,0.0,0.0,.75,.75,.75,.05,1.4);		// Top
//  Scale(o,500,500,500);
//  Translate(o,0,510,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  // Two spheres scene
//  o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,1.54);		// Refract
//  Scale(o,3.75,3.75,3.75);
//  Translate(o,-5,-4.0,4.5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);


// // Two spheres scene
//  o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,1.54);		// Refract
//  Scale(o,3.75,3.75,3.75);
//  Translate(o,-5,-4.0,4.5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  o=newSphere(0.0,1.0,0.0,.99,.99,.99,.05,2.47);		// Reflect
//  Scale(o,3.75,3.75,3.75);
//  Translate(o,4,-3.75,6.5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

/*
 // Ring of refracting spheres
 for (int i=0; i<5;i++)
 {
  o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,1.45+(.1*i));
  Scale(o,1.75,1.75,1.75);
  Translate(o,3.25*cos(2*PI*i/5),-2.45,3+3.25*sin(2*PI*i/5));
  invert(&o->T[0][0],&o->Tinv[0][0]);
  insertObject(o,&object_list);
 }

 for (int i=0; i<7;i++)
 {
  o=newSphere(0.0,0.0,1.0,.99,.99,.99,.01,2.00+(.05*i));
  Scale(o,1.75,1.75,1.75);
  Translate(o,4.60*cos(2*PI*i/7),-6.35,3+4.60*sin(2*PI*i/7));
  invert(&o->T[0][0],&o->Tinv[0][0]);
  insertObject(o,&object_list);
 }
*/

//  // Two spheres scene
//  o=newSphere(0.0,1.0,0.0,.99,.99,.99,.05,2.47);		// Reflect
//  Scale(o,3,.5,3);
//  Translate(o,0,7.5,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);


//  // Planar light source at top
//  o=newPlane(1.00,0.00,0.0,1.0,1.0,1.0,0.0,1.54);
//  Scale(o,.5,2.5,.1);
//  RotateX(o,PI/2);
//  Translate(o,0,9.995,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  o->isLightSource=1;
//  insertObject(o,&object_list);





 /// TRYING NEW SCENE

 // Cornell box
 o=newSphere(0.0,1.0,0.0,0.1,.3,.5,.05,1.4);	// Left
 Scale(o,500,500,500);
 Translate(o,-512,0,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(0.0,1.0,0.0,0.1,.3,.5,.05,1.4);		// Right
 Scale(o,500,500,500);
 Translate(o,512,0,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(0.0,1.0,0.0,0.1,.3,.5,.05,1.4);		// Back
 Scale(o,500,500,500);
 Translate(o,0,0,515);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,1.0,1.0,.75,.75,.9,.02,1.4);	// Bottom
 Scale(o,500,500,500);
 Translate(o,0,-510,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

//  o=newSphere(0.0,1.0,0.0,.99,.99,.99,.05,2.47);		// Bottom
//  Scale(o,500,500,500);
//  Translate(o,0,-510,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

 o=newSphere(1.0,0.0,0.0,.99,.99,.99,.05,1.4);		// Top
 Scale(o,500,500,500);
 Translate(o,0,510,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);


//HELMET
 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);		
 Scale(o,2,2,2);
 Translate(o,0,3.75,3);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);	
 insertObject(o,&object_list);

 // Helmet Details
 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);		
 Scale(o,.3,.3,.3);
 Translate(o,2,3.8,3);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);	
 insertObject(o,&object_list);


 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);		
 Scale(o,.3,.3,.3);
 Translate(o,-2,3.8,3);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);	
 insertObject(o,&object_list);

//SHINY PART
  o=newSphere(0.0,1.0,0.0,.99,.99,.99,.05,2.47);		
 Scale(o,2,1.7,2.3);
 Translate(o,0,3.75,3);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 //Neck
//  o=newCyl(1.0,0,0.0,.99,.99,.99,.05,2.47);
//  RotateX(o,PI/2);		
//  Translate(o,0,2.5,3);
//  //RotateZ(o,PI/6);
 RotateY(o,PI/6);
//   RotateY(o,PI/6);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

 //Chest
 o=newCyl(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,2.5,1.75,4); //z is length
 RotateX(o,PI/3);		
 Translate(o,0,-1,5);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,2.5,3,2);
 RotateX(o,-PI/3);		
 Translate(o,0,-1,3.8);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,.4,.4,.4);
 Translate(o,1.3,-0.3,1.7);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 //Shoulders
 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,2,1.5,2);		
 Translate(o,2,.5,4);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);


 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,2,1.5,2);		
 Translate(o,-2,.5,4);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 //Arms
 //right
 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,1,3,2);	
 RotateX(o,-PI/4);		
 Translate(o,3,-2,4.5);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,1,2.5,1);
 RotateX(o,-2*PI/3);		
 Translate(o,3.4,-4.25,2);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 //left
 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,1,3.3,1.5);
 RotateX(o,PI/3);		
 Translate(o,-3.2,-1.5,2);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,1,2.5,1);
 RotateX(o,-PI/3);		
 Translate(o,-3.4,-2.5,-0.5);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 //Hand time D:
 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,0.75,1,0.5);	
 RotateX(o,-PI/3);		
 Translate(o,-3.4,-1.25,-3.5);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,0.75,1,0.5);	
 RotateX(o,-2*PI/3);		
 Translate(o,3.4,-5.75,-1);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 //Fingies

 double finger_space[4] = {-0.25, 0, 0.25, 0.5};

 for (int i=0; i<4;i++)
 {
    //Left hand
    o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
    Scale(o,0.25,1,0.25);	
    RotateX(o,-PI/3);		
    Translate(o,-3.4+finger_space[i],-0.8,-4.3);
    //RotateZ(o,PI/6);
 RotateY(o,PI/6);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

    o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
    Scale(o,0.25,0.5,0.25);	
    RotateX(o,PI/3);		
    Translate(o,-3.4+finger_space[i],-0.6,-5.5);
    //RotateZ(o,PI/6);
 RotateY(o,PI/6);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

    //Right Hand
    o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
    Scale(o,0.25,0.5,0.25);	
    RotateX(o,-2*PI/3);			
    Translate(o,3.4+finger_space[i],-6.25,-1.9);
    //RotateZ(o,PI/6);
 RotateY(o,PI/6);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

    o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
    Scale(o,0.25,1,0.25);	
    RotateX(o,-PI/3);		
    Translate(o,3.4+finger_space[i],-6.75,-1.4);
    //RotateZ(o,PI/6);
 RotateY(o,PI/6);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    insertObject(o,&object_list);

 }

 // THUMBS
 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,0.25,0.5,0.25);	
 RotateZ(o,-PI/3);		
 Translate(o,-2.5,-1,-4);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,0.25,0.5,0.25);	
 RotateZ(o,-PI/4);
 Translate(o,2.5,-6.5,-1);		
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // LOWER BODY BEGINS

 // Tummyyyy

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,2.5,3,2);
 RotateX(o,-PI/3);		
 Translate(o,0,-3,4.75);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);
 
 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,2.5,3,2);
 RotateX(o,-PI/3);		
 Translate(o,0,-4.5,5.5);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newCyl(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,3,3,1); //z is length
 RotateX(o,PI/3);		
 Translate(o,0,-3.5,5.5);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 // LEGS

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,1.5,3.7,1.7);
 RotateX(o,PI/3);		
 Translate(o,-1.7,-6.3,4);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,1.5,3.7,1.7);
 RotateX(o,-PI/6);		
 Translate(o,1.7,-7,6);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);


 //KNEES

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,1.3,1.3,1.3);
 RotateX(o,PI/3);		
 Translate(o,-1.7,-8,2);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newSphere(1.0,0,0.0,.99,.99,.99,.05,2.47);
 Scale(o,1.3,1.3,1.3);
 RotateX(o,-PI/6);		
 Translate(o,1.7,-9,7);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);



 // PLANETS

 o=newSphere(0.0,.5,.6,.6,.25,.75,.05,2.47); // pink
 Translate(o,0,0,-2);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);


 o=newSphere(0.0,.5,.6,1,0,0,.05,2.47); // red 
 Translate(o,-6,3,2);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newCyl(0.0,.5,.6,1,.53,.0,.05,2.47); // yellow ring
 Scale(o,1.5,1.5,.1);
 RotateX(o,PI/2);
 Translate(o,-6,3,2);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);


 o=newSphere(0.0,.5,.5,.5,.7,.9,.05,2.47); //blue
 Scale(o,1,1,1);
 Translate(o,6,-1,-4);
 //RotateZ(o,PI/6);
 RotateY(o,PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);



 // Planar light source at top
 o=newPlane(1.00,0.00,0.0,1.0,1.0,1.0,0.0,1.54);
 Scale(o,.5,2.5,.1);
 RotateX(o,PI/2);
 Translate(o,0,9.995,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 o->isLightSource=1;
 insertObject(o,&object_list);












 

}