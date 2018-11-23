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
 o=newSphere(.65,.95,.85,0.5,0,.7,.7,1,1,10);		// Inner
 Scale(o,6,6,1);		
 Translate(o,0,21,30);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.55,.85,.65,0.45,0,.6,.7,.9,1.1,9);		// left 1
 Scale(o,6,6,1);		
 Translate(o,-8,18,27);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.55,.85,.65,0.45,0,.7,.6,.9,1.1,9);		// right 1
 Scale(o,6,6,1);		
 Translate(o,8,18,27);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.45,.75,.55,0.4,0,.5,.7,.8,1.2,8);		// left 2
 Scale(o,6,6,1);		
 Translate(o,-12,15,24);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.45,.75,.55,0.4,0,.7,.5,.8,1.2,8);		// right 2
 Scale(o,6,6,1);		
 Translate(o,12,15,24);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.35,.65,.45,0.35,0,.4,.7,.7,1.3,7);		// left 3
 Scale(o,6,6,1);		
 Translate(o,-15,11,19);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.35,.65,.45,0.35,0,.7,.4,.7,1.3,7);		// right 3
 Scale(o,6,6,1);		
 Translate(o,15,11,19);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.25,.55,.35,0.3,0,.3,.7,.6,1.4,6);		// left 4
 loadTexture(o,"./Texture/alpha_map.pgm",3,&texture_list);
o->alphaMapped=1;
 Scale(o,6,6,1);		
 Translate(o,-17,7,14);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.25,.55,.35,0.3,0,.7,.3,.6,1.4,6);		// right 4
 loadTexture(o,"./Texture/alpha_map.pgm",3,&texture_list);
o->alphaMapped=1;
 Scale(o,6,6,1);		
 Translate(o,17,7,14);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.15,.45,.25,0.2,0,.2,.7,.5,1.5,5);		// left 5
loadTexture(o,"./Texture/bump1.ppm",2,&texture_list);
o->normalMapped=1;
 Scale(o,6,6,1);		
 Translate(o,-19,4,10);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);

 o=newSphere(.15,.45,.25,0.2,0,.7,.2,.5,1.5,5);		// right 5
loadTexture(o,"./Texture/bump1.ppm",2,&texture_list);
o->normalMapped=1;
 Scale(o,6,6,1);		
 Translate(o,19,4,10);
 invert(&o->T[0][0],&o->Tinv[0][0]);		
insertObject(o,&object_list);


//  // If needed, this is how you load a texture map
 	
// This loads a texture called 'mosaic2.ppm'. The
// 								// texture gets added to the texture list, and a
// 								// pointer to it is stored within this object in the
// 								// corresponding place. The '1' indicates this image
// 								// will be used as a texture map. Use '2' to load
// 								// an image as a normal map, and '3' to load an
// 								// alpha map. Texture and normal maps are RGB .ppm
// 								// files, alpha maps are grayscale .pgm files.
// 								// * DO NOT * try to free image data loaded in this
// 								// way, the cleanup function already provided will do
// 								// this at the end.
 

 o=newPlane(.5,.75,.05,.05,.6,.35,.33,1,1,2); 
loadTexture(o,"./Texture/normalbrick.ppm",2,&texture_list);
o->normalMapped = 1;
RotateX(o,PI/4);
 Scale(o,37,37,37);
 Translate(o,0,0,15);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 insertObject(o,&object_list);

 o=newPlane(.75,.75,.75,.75,.9,.9,.9,1,1,5);
 Scale(o,2,2,2);
 RotateX(o,PI/4);
 Translate(o,0,30.5,13.5);
 invert(&o->T[0][0],&o->Tinv[0][0]);
 o->isLightSource=1;
 insertObject(o,&object_list);

 // Insert a single point light source. We set up its position as a point structure, and specify its
 // colour in terms of RGB (in [0,1]).
//  p.px=0;
//  p.py=25.5;
//  p.pz=-3.5;
//  p.pw=1;
//  l=newPLS(&p,.95,.95,.95);
//  insertPLS(l,&light_list);

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
