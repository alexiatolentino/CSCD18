 struct object3D *o;
 struct pointLS *l;
 struct point3D p;

 // ORIGINAL BUILD SCENE FROM A2 STARTER:

 // Shows ALS

 // o=newSphere(.05,.95,.35,.35,1,.25,.25,1,1,6);		// Initialize a sphere
 // Scale(o,1.5,.75,.75);					// Apply a few transforms (Translate * Rotate * Scale)
 // RotateZ(o,PI/4);					
 // Translate(o,2.0,2.5,1.5);
 // invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
 
 //  insertObject(o,&object_list);			// <-- If you don't insert the object into the object list,
	// 					//     nothing happens! your object won't be rendered.

 // // That's it for defining a single sphere... let's add a couple more objects
 // o=newSphere(.05,.95,.95,.75,.75,.95,.55,1,1,6);
 // Scale(o,.95,1.65,.65);
 // RotateZ(o,-PI/1.5);
 // Translate(o,-2.2,1.75,1.35);
 // invert(&o->T[0][0],&o->Tinv[0][0]);
 // insertObject(o,&object_list);

 // o=newPlane(.05,.75,.05,.05,.55,.8,.75,1,1,2);
 // Scale(o,11,11,11);
 // RotateZ(o,PI/4);
 // RotateX(o,PI/2);
 // Translate(o,0,-4,5);
 // invert(&o->T[0][0],&o->Tinv[0][0]);
 // insertObject(o,&object_list);

 // // Insert area light source
 // addAreaLight(.2, .2, .2, 0, 0, 0, 5, .5, 1, 1, 1, 1, &object_list, &light_list, 1);


 // VERY SIMPLE SCENE (texture + alpha mapping tester):

//  o=newSphere(.05,.95,.5,.5,1,1,1,1,1,6);		// Initialize a sphere
//  Scale(o,1.1,.5,.25);				// Apply a few transforms (Translate * Rotate * Scale)			
//  Translate(o,1.7,2.05,0.9);
//  RotateZ(o,PI/1.5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
//  loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
//  insertObject(o,&object_list);

//  o=newSphere(.05,.95,.5,.5,1,1,1,1,1,6);		// Initialize a sphere
//  Scale(o,.5,.5,.5);	
//  Translate(o,2.5,-1,1.5);
//  RotateY(o,-PI/20);
//  invert(&o->T[0][0],&o->Tinv[0][0]);	
//  loadTexture(o,"./Texture/alphaMap.pgm",3,&texture_list);
//  insertObject(o,&object_list);	

 
// // Plane
//  o=newPlane(.05,.75,.05,.05,1,1,1,1,1,2);
//  Scale(o,11,11,11);
//  RotateZ(o,PI/4);
//  RotateX(o,PI/2);
//  Translate(o,0,-4,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  loadTexture(o,"./Texture/Water_texture.ppm",1,&texture_list);
//  insertObject(o,&object_list);


//  // Insert a single point light source. We set up its position as a point structure, and specify its
//  // colour in terms of RGB (in [0,1]).
//  p.px=0;
//  p.py= 25.5;
//  p.pz=-3.5;
//  p.pw=1;
//  l=newPLS(&p,0.85,1,.92);
//  insertPLS(l,&light_list);


 // SIMPLE DOG SCENE:
 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)

//  //HEAD OF DOG
//  //back ear
//  o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
//  Scale(o,1.3,.5,.5);					// Apply a few transforms (Translate * Rotate * Scale)			
//  Translate(o,-3.9,1,2);
//  //RotateZ(o,-PI/1.49);
//  RotateZ(o,-PI/1.9);
//  invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
//  insertObject(o,&object_list);
 
//  //front ear
//  o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
//  Scale(o,1.3,.5,.5);					// Apply a few transforms (Translate * Rotate * Scale)			
//  Translate(o,-3.45,1.35,1.4);
//  RotateZ(o,-PI/1.9);
//  invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
//  insertObject(o,&object_list);


//  //neck
//  o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
//  Scale(o,1.3,.5,.5);					// Apply a few transforms (Translate * Rotate * Scale)			
//  Translate(o,-1.2,1.35,1.7);
//  RotateZ(o,-PI/1.9);
//  invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
//  insertObject(o,&object_list);

 
//  //head
//  o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
//  Scale(o,1.1,.5,.25);					// Apply a few transforms (Translate * Rotate * Scale)			
//  Translate(o,-2.5,-2,1.4);
//  RotateZ(o,PI);
//  invert(&o->T[0][0],&o->Tinv[0][0]);	
//  insertObject(o,&object_list);	


//  //TAIL OF DOG
//  o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
//  Scale(o,1.1,.5,.25);				// Apply a few transforms (Translate * Rotate * Scale)			
//  Translate(o,1.7,2.05,0.9);
//  RotateZ(o,PI/1.5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
//  insertObject(o,&object_list);

//  //body
//  o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
//  //Scale(o,1.7,.6,.3);					// Apply a few transforms (Translate * Rotate * Scale)			
//  Translate(o,-.1,-.4,1.2);
//  RotateY(o,-PI/20);
//  invert(&o->T[0][0],&o->Tinv[0][0]);	
//  loadTexture(o,"./Textures/alphaMap.ppm",3,&texture_list);
//  insertObject(o,&object_list);	

//  o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
//  Scale(o,1.3,.5,.5);
//  // Apply a few transforms (Translate * Rotate * Scale)			
//  Translate(o,2.1,1,1.4);
//  RotateZ(o,-2*PI/6);
//  invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
//  insertObject(o,&object_list);

//  o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
//  Scale(o,1.3,.5,.5);
//  //-0.45 diff, 0.35 diff, 0.5 diff
//  // Apply a few transforms (Translate * Rotate * Scale)			
//  Translate(o,-2.7,1.3,0.9);
//  RotateZ(o,2*PI/6);
//  invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
//  insertObject(o,&object_list);

 
//  //PLANE

//  o=newPlane(.05,.75,.05,.05,1,1,1,1,1,2);
//  Scale(o,11,11,11);
//  RotateZ(o,PI/4);
//  RotateX(o,PI/2);
//  Translate(o,0,-4,5);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  // Insert a single point light source. We set up its position as a point structure, and specify its
//  // colour in terms of RGB (in [0,1]).
//  p.px=0;
//  p.py=25.5;
//  p.pz=-3.5;
//  p.pw=1;
//  l=newPLS(&p,.95,.95,.95);
//  insertPLS(l,&light_list);




// COOL HELLO KITTY SCENE!!
//FLOOR
o = newPlane(.3, .75, .05, .05, 0.5, 0.2, 0.7, 1, 1, 2);
Scale(o, 18, 11, 11);
RotateZ(o, PI);
RotateX(o, PI / 2);
Translate(o, 0, -4, 2);
invert(&o->T[0][0], &o->Tinv[0][0]);
loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
insertObject(o, &object_list);

// HEAD
o = newSphere(.05, .95, .95, .8, 1, 1, 1, 0.9, 1, 0.5);
Scale(o, 0.8, 0.6, 0.4);
Translate(o, 0, 0.5, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

// NOSE
o = newSphere(.05, .95, .95, .8, 0.91, 0.79, 0, 1, 1, 6);
Scale(o, 0.05, 0.05, 0.1);
Translate(o, 0,.7, -0.5);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

// // EYES
// o = newSphere(0, 0, 0, 0, 0.3, 0.3, 0.6, 1, 1, 6);
// Scale(o, 0.15, 0.2, 0.2);
// Translate(o, -0.5, 0.5, -0.1);
// invert(&o->T[0][0], &o->Tinv[0][0]);
// insertObject(o, &object_list);

// o = newSphere(0, 0, 0, 0, 0.3, 0.3, 0.6, 1, 1, 6);
// Scale(o, 0.15, 0.2, 0.2);
// Translate(o, 0.5, 0.5, -0.1);
// invert(&o->T[0][0], &o->Tinv[0][0]);
// insertObject(o, &object_list);

/// EARS
o = newSphere(.05, .95, .95, .8, 1, 1, 1, 1, 1, 1);
Scale(o, 0.15, 0.2, 0.2);
RotateY(o, PI/8);
RotateZ(o, PI/8);
Translate(o, -0.55, 1, -0.1);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newSphere(.05, .95, .95, .8, 1, 1, 1, 1, 1, 1);
Scale(o, 0.15, 0.2, 0.2);
RotateY(o, -PI/8);
RotateZ(o, -PI/8);
Translate(o, 0.55,1, -0.1);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

//BOW
o = newSphere(1, 1, 0, .8, 1, 0.01, 0.01, 0.9, 1, 0);
Scale(o, 0.1, 0.2, 0.2);
RotateY(o, PI/4);
RotateZ(o, PI/8);
Translate(o, 0.71,0.8, -0.1);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);


// BODY
o = newSphere(.05, .95, .95, .8, 1, 1, 1, 0.9, 1, 0);
Scale(o, 0.6, 0.62, 0.4);
Translate(o, 0, -0.7, 0);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);


// HANDS
o = newSphere(.05, .95, .95, .8, 1, 1, 1, 0.9, 1, 0.5);
Scale(o, 0.3, 0.2, 0.2);
RotateY(o, PI/8);
RotateZ(o, PI/8);
Translate(o, 0.75,-0.3, -0.1);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newSphere(.05, .95, .95, .8, 1, 1, 1, 0.9, 1, 0.5);
Scale(o, -0.3, 0.2, 0.2);
RotateY(o, -PI/8);
RotateZ(o, -PI/8);
Translate(o, -0.7,-0.3, -0.1);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

// FEET
o = newSphere(.05, .95, .95, .8, 1, 1, 1, 0.9, 1, 0.5);
Scale(o, 0.2, 0.1, 0.1);
RotateY(o, PI/8);
RotateZ(o, PI/8);
Translate(o, 0.75,-1.2, 0.1);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newSphere(.05, .95, .95, .8, 1, 1, 1, 0.9, 1, 0.5);
Scale(o, 0.2, 0.1, 0.1);
RotateY(o, -PI/8);
RotateZ(o, -PI/8);
Translate(o, -0.75,-1.2, 0.1);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);



/********************************************************/
// End of simple scene for Assignment 2
// Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
// or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
//
// Remember: A lot of the quality of your scene will depend on how much care you have put into defining
//           the relflectance properties of your objects, and the number and type of light sources
//           in the scene.

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// TO DO: For Assignment 3 you *MUST* define your own cool scene.
//	   We will be looking for the quality of your scene setup, the use of hierarchical or composite
//	   objects that are more interesting than the simple primitives from A2, the use of textures
//        and other maps, illumination and illumination effects such as soft shadows, reflections and
//        transparency, and the overall visual quality of your result. Put some work into thinking
//        about these elements when designing your scene.
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// WALLS
o = newPlane(.9, .4, .05, 0.2, 0.97, 0.56, 0.69, 1, 1, 0.5);
Scale(o, 10, 6.75, 7);
RotateZ(o, PI / 2);
RotateY(o, PI / 2);
RotateX(o, PI / 2);
Translate(o, 10, 3, 8);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

o = newPlane(.9, .4, .05, 0.2, 0.97, 0.56, 0.69, 1, 1, 0.5);
Scale(o, 10, 6.75, 7);
RotateZ(o, PI / 2);
RotateY(o, PI / 2);
RotateX(o, PI / 2);
Translate(o, -10, 3, 8);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

// TOP
o = newPlane(.3, .75, .05, .05, 0.69, 0.08, 0.33, 1, 1, 2);
Scale(o, 15, 11, 11);
RotateZ(o, PI);
RotateX(o, PI / 2);
Translate(o, 0, 11.5, 5);
invert(&o->T[0][0], &o->Tinv[0][0]);
insertObject(o, &object_list);

// Insert a single point light source. We set up its position as a point structure, and specify its
// colour in terms of RGB (in [0,1]).

p.px = 0;
p.py = 7;
p.pz = -3.5;
p.pw = 1;
l = newPLS(&p, .95, .95, .95);
insertPLS(l, &light_list);

