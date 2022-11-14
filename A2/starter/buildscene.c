 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
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
 struct pointLS *l;
 struct point3D p;

 // Simple scene for Assignment 3:
 // Insert a couple of objects. A plane and two spheres
 // with some transformations.

 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)
 
 //HEAD OF DOG
 //back ear
 o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
 Scale(o,1.3,.5,.5);					// Apply a few transforms (Translate * Rotate * Scale)			
 Translate(o,-3.9,1,2);
 //RotateZ(o,-PI/1.49);
 RotateZ(o,-PI/1.9);
 invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
 loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
 insertObject(o,&object_list);
 
 //front ear
 o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
 Scale(o,1.3,.5,.5);					// Apply a few transforms (Translate * Rotate * Scale)			
 Translate(o,-3.45,1.35,1.4);
 RotateZ(o,-PI/1.9);
 invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
 loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
 insertObject(o,&object_list);


 //neck
 o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
 Scale(o,1.3,.5,.5);					// Apply a few transforms (Translate * Rotate * Scale)			
 Translate(o,-1.2,1.35,1.7);
 RotateZ(o,-PI/1.9);
 invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
 loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
 insertObject(o,&object_list);

 
 //head
 o=newSphere(.05,.95,.5,.5,0.63,.12,.94,1,1,6);		// Initialize a sphere
 Scale(o,1.1,.5,.25);					// Apply a few transforms (Translate * Rotate * Scale)			
 Translate(o,-2.5,-2,1.4);
 RotateZ(o,PI);
 invert(&o->T[0][0],&o->Tinv[0][0]);	
 loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
 insertObject(o,&object_list);	


 //TAIL OF DOG
 o=newSphere(.05,.95,.5,.5,0.63,.12,.94,0.5,1.5,6);		// Initialize a sphere
 Scale(o,1.1,.5,.25);				// Apply a few transforms (Translate * Rotate * Scale)			
 Translate(o,1.7,2.05,0.9);
 RotateZ(o,PI/1.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
 loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
 insertObject(o,&object_list);

 //body
 o=newSphere(.05,.95,.5,.5,0.63,.12,.94,0.5,1.5,6);		// Initialize a sphere
 Scale(o,1.7,.6,.3);					// Apply a few transforms (Translate * Rotate * Scale)			
 Translate(o,-.1,-.4,1.2);
 RotateY(o,-PI/20);
 invert(&o->T[0][0],&o->Tinv[0][0]);	
 loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
 insertObject(o,&object_list);	

 o=newSphere(.05,.95,.5,.5,0.63,.12,.94,0.5,1.5,6);		// Initialize a sphere
 Scale(o,1.3,.5,.5);
 // Apply a few transforms (Translate * Rotate * Scale)			
 Translate(o,2.1,1,1.4);
 RotateZ(o,-2*PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
 loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
 insertObject(o,&object_list);

 o=newSphere(.05,.95,.5,.5,0.63,.12,.94,0.5,1.5,6);		// Initialize a sphere
 Scale(o,1.3,.5,.5);
 //-0.45 diff, 0.35 diff, 0.5 diff
 // Apply a few transforms (Translate * Rotate * Scale)			
 Translate(o,-2.7,1.3,0.9);
 RotateZ(o,2*PI/6);
 invert(&o->T[0][0],&o->Tinv[0][0]);			// Compute the inverse transform * DON'T FORGET TO DO THIS! *
 loadTexture(o,"./Texture/wavy.ppm",1,&texture_list);
 insertObject(o,&object_list);



 // If needed, this is how you load a texture map
 // loadTexture(o,"./Texture/mosaic2.ppm",1,&texture_list);	// This loads a texture called 'mosaic2.ppm'. The
								// texture gets added to the texture list, and a
								// pointer to it is stored within this object in the
								// corresponding place. The '1' indicates this image
								// will be used as a texture map. Use '2' to load
								// an image as a normal map, and '3' to load an
								// alpha map. Texture and normal maps are RGB .ppm
								// files, alpha maps are grayscale .pgm files.
								// * DO NOT * try to free image data loaded in this
								// way, the cleanup function already provided will do
								// this at the end.
 
		// <-- If you don't insert the object into the object list,
						//     nothing happens! your object won't be rendered.

 // That's it for defining a single sphere... let's add a couple more objects

 o=newPlane(.05,.75,.05,.05,1,1,1,1,1,2);
 Scale(o,11,11,11);
 RotateZ(o,PI/4);
 RotateX(o,PI/2);
 Translate(o,0,-4,5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);


 // Insert a single point light source. We set up its position as a point structure, and specify its
 // colour in terms of RGB (in [0,1]).
 p.px=0;
 p.py=25.5;
 p.pz=-3.5;
 p.pw=1;
 l=newPLS(&p,.95,.95,.95);
 insertPLS(l,&light_list);

 // End of simple scene for Assignment 2
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
 
//ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)

//  o=newSphere(1,.5,.95,.75,0,1,0,1,1,6);
//  //Scale(o,.95,1.65,.65);
//  //RotateZ(o,-PI/1.5);
//  Translate(o,3,3,3);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  o=newSphere(1,.5,.95,.75,0,1,0,1,1,6);
//  Scale(o,-0.5,-0.5,-0.5);
//  //RotateZ(o,-PI/1.5);
//  Translate(o,3,3,4);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);

//  o=newCyl(0.4,.5,.5,.5,0,1,0,1,1,6);
//  Scale(o,2,0,0);
//  RotateZ(o,PI/2);
//  Translate(o,-1,-1,3);
//  invert(&o->T[0][0],&o->Tinv[0][0]);
//  insertObject(o,&object_list);
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
 // TO DO: For Assignment 3 you *MUST* define your own cool scene.
 //	   We will be looking for the quality of your scene setup, the use of hierarchical or composite
 //	   objects that are more interesting than the simple primitives from A2, the use of textures
 //        and other maps, illumination and illumination effects such as soft shadows, reflections and
 //        transparency, and the overall visual quality of your result. Put some work into thinking
 //        about these elements when designing your scene.
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
