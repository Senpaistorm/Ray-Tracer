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
* 1) Student Name: Xiang Li		
* 2) Student Name:		
*
* 1) Student number: 1002486927
* 2) Student number:
* 
* 1) UtorID lixian51
* 2) UtorID
* 
* We hereby certify that the work contained here is our own
*
* ____________Xiang Li________             _____________________
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
  struct point3D pts, ptc, ptm;
  struct pointLS *ls = light_list;
  struct object3D *currentObj = object_list;
  double lda=-1;
  struct ray3D *lightray = newRay(newPoint(0,0,0), newPoint(0,0,0));
  struct object3D *objhit;
  struct point3D ptp,ptn;
  double K = 10;
  double k = 0;
  double distance;
  struct ray3D *r_ray = newRay(p,n);
  struct point3D r_b;struct point3D lsp0;
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
 // normal and alpha mapping
  if(obj->normalMapped){
    obj->textureMap(obj->normalMap,a,b,&n->px,&n->py,&n->pz);
  }
  if(obj->alphaMapped){
    obj->aMap(obj->alphaMap,a,b,&obj->alpha);
  }
  // ambient component
  tmp_col.R += obj->alb.ra *R;
  tmp_col.G += obj->alb.ra *G;
  tmp_col.B += obj->alb.ra *B;

  ptc.px = -(ray->d.px);
  ptc.py = -(ray->d.py);
  ptc.pz = -(ray->d.pz);
  ptc.pw = 1;
  // loop through each point light source and cast shadow rays
  while(ls != NULL){
    pts.px = ls->p0.px - p->px;
    pts.py = ls->p0.py - p->py;
    pts.pz = ls->p0.pz - p->pz;
    pts.pw = 1;
    distance = length(&pts);
    normalize(&pts);
    // ray trace in the direction of the lightsource
    // and scale the lambda down by the distance to the lightsource
    lightray = newRay(p, &pts);
    lightray->inside = 0;
    findFirstHit(lightray, &lda, obj, &objhit, &ptp,&ptn,&a,&b);
    lda = lda/distance;
    if(lda < 0 || lda > 1){
      phongDS(ls->col,pts,ptc,obj,R,G,B,n,&tmp_col,1.0);
    }
    ls = ls->next;
  }
  // loop through list of objects, if object is an area light source, cast K shadow rays randomly
  // and apply multiply k/K to diffuse and specular component
  while(currentObj != NULL){
    if(currentObj->isLightSource == 1){
      k=0;
      for(int i = 1; i <= K;i++){
        currentObj->randomPoint(currentObj, &lsp0.px,&lsp0.py,&lsp0.pz);
        //fprintf(stderr,"point on the area ls %f %f %f",lsp0.px,lsp0.py,lsp0.pz);
        pts.px = lsp0.px - p->px;
        pts.py = lsp0.py - p->py;
        pts.pz = lsp0.pz - p->pz;
        pts.pw = 1;
        distance = length(&pts);
        normalize(&pts);
        lightray = newRay(p, &pts);
        lightray->inside = 0;
        findFirstHit(lightray, &lda, obj, &objhit, &ptp,&ptn,&a,&b);
        lda = lda/distance;
        free(lightray);
        if(lda < 0 || lda > 1 || objhit->isLightSource == 1){
          k += 1;
        }
      }
        
      phongDS(currentObj->col,pts,ptc,obj,R,G,B,n,&tmp_col,k/K);
    }
    currentObj = currentObj->next;
  }
  
  // scale local component light by alpha
  multColor(&tmp_col,obj->alpha);
  
  if(depth < MAX_DEPTH){
    // refraction component if needed
    if(obj->alpha < 1.0){
      refract(obj,ray,r_ray,n,depth,&tmp_col);
    }
    // global component if needed
    if(!ray->inside){
      global(obj,ray,p,n,depth,&tmp_col);
    }
  }
  // update 'col' with the final colour computed
  col->R = tmp_col.R > 1.0 ? 1.0 : tmp_col.R;
  col->G = tmp_col.G > 1.0 ? 1.0 : tmp_col.G;
  col->B = tmp_col.B > 1.0 ? 1.0 : tmp_col.B;
  return;

}

void multColor(struct colourRGB *col, double r){
  col->R = col->R*r;
  col->G = col->G*r;
  col->B = col->B*r;
}

void phongDS(struct colourRGB col, struct point3D pts,struct point3D ptc,struct object3D *obj, double R, double G, double B, struct point3D *n, struct colourRGB *tmp_col,double k){
  // Calculates the phong components diffuse and specular and store it in col
  double diffuse,specular;
  // diffuse component
  struct point3D ptm;
  diffuse = dot(n,&pts);
  if(diffuse < 0){
    diffuse = 0;
  }
  tmp_col->R +=  obj->alb.rd * R * col.R * diffuse *k;
  tmp_col->G +=  obj->alb.rd * G * col.G * diffuse*k;
  tmp_col->B +=  obj->alb.rd * B * col.B * diffuse*k;
  
  //specular component
  ptm.px = (2*(dot(&pts,n))*(n->px)) - pts.px;
  ptm.py = (2*(dot(&pts,n))*(n->py)) - pts.py;
  ptm.pz = (2*(dot(&pts,n))*(n->pz)) - pts.pz;
  ptm.pw = 1;
  normalize(&ptm);
  specular = dot(&ptm,&ptc);
  if(specular < 0){
    specular = 0;
  }else{
    specular = pow(specular, obj->shinyness);
  }
  tmp_col->R +=  obj->alb.rs * col.R * specular*k;
  tmp_col->G +=  obj->alb.rs * col.G * specular*k;
  tmp_col->B +=  obj->alb.rs * col.B * specular*k;
}

void refract(struct object3D *obj, struct ray3D *ray, struct ray3D *r_ray, struct point3D *n, int depth, struct colourRGB *tmp_col){
  double n1,n2,c,r,rfactor;
  struct point3D r_b;
  struct colourRGB refractRGB;
  struct point3D normal;
  double thetac;
  // scale refracted light by 1-alpha
  // calculate r=n1/n2
  if(ray->inside){
    r_ray->inside = 0;
    n1 = obj->r_index;
    n2 = 1.0;
    normal.px = -n->px;
    normal.py = -n->py;
    normal.pz = -n->pz;
    normal.pw = 1;
  }else{
    r_ray->inside = 1;
    n1 = 1.0;
    n2 = obj->r_index;
    normal.px = n->px;
    normal.py = n->py;
    normal.pz = n->pz;
    normal.pw = 1;
  }
  r = n1/n2;
  // calculate c = -n.b
  c = -1*dot(&normal,&ray->d);
  // d = rb+(rc-sqrt(1-r^2(1-c^2)))
  rfactor = r*c - sqrtl((1-r*r*(1 - c*c)));
  r_ray->d.px = normal.px * rfactor + ray->d.px * r;
  r_ray->d.py = normal.py * rfactor + ray->d.py * r;
  r_ray->d.pz = normal.pz * rfactor + ray->d.pz * r;

  if(ray->inside){
    r_ray->inside = 0;
  }else{
    r_ray->inside = 1;
  }
  if(ray->inside && n1>n2){
    // total internal reflection
    thetac = asin(n2/n1);
    if(acos(dot(&ray->d,n)) >= thetac){
      refractRGB.R=0;
      refractRGB.G=0;
      refractRGB.B=0;
    }else{
      rayTrace(r_ray, depth + 1, &refractRGB, obj);
    }
  }else{
    rayTrace(r_ray, depth + 1, &refractRGB, obj); 
  }
  tmp_col->R += refractRGB.R * (1-obj->alpha);
  tmp_col->G += refractRGB.G * (1-obj->alpha);
  tmp_col->B += refractRGB.B * (1-obj->alpha);
}

void global(struct object3D *obj, struct ray3D *ray, struct point3D *p, struct point3D *n, int depth, struct colourRGB *tmp_col){
  struct colourRGB globalRGB;
  struct point3D ptr;
  // global component
  // r = -2(d.n)n + d
  ptr.px = ray->d.px - (2*dot(&ray->d,n)*(n->px)) ;
  ptr.py = ray->d.py - (2*dot(&ray->d,n)*(n->py)) ;
  ptr.pz = ray->d.pz - (2*dot(&ray->d,n)*(n->pz)) ;
  ptr.pw = 1;
  normalize(&ptr);
  struct ray3D *new_ray = newRay(p,&ptr);
  new_ray->inside = 0;
  // recursive call to ray trace from intersection point
  rayTrace(new_ray, depth + 1, &globalRGB, obj);
  tmp_col->R += obj->alb.rg * globalRGB.R;
  tmp_col->G += obj->alb.rg * globalRGB.G;
  tmp_col->B += obj->alb.rg * globalRGB.B;
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

  *lambda = -1;
  struct point3D point;
  struct point3D normal;
  double ld;

  // store the current object we are looking at
  struct object3D *cur_obj = object_list;
  // loop through each object in the object list
  while(cur_obj != NULL){
    // center of projection, check against every object
    if(cur_obj != Os){
      cur_obj->intersect(cur_obj, ray, &ld, &point, &normal, a, b);
      if((*lambda < 0 || ld < *lambda) && ld > 0){
        *p = point;
        *n = normal;
        *lambda = ld;
        *obj = cur_obj;
      }
    }
    cur_obj = cur_obj->next;
  }
  return;
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Trace one ray through the scene.
 //
 // Parameters:
 //   *ray   -  A pointer to the rNULLay being traced
 //   depth  -  Current recursion NULLdepth for recursive raytracing
 //   *col   - Pointer to an RGB cNULLolour structure so you can return the object colour
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

  findFirstHit(ray,&lambda,Os,&obj,&p,&n,&a,&b);

  if(lambda > 0 && obj != Os){
    rtShade(obj,&p,&n,ray,depth,a,b,&I);
  }else {
    // background
    I.R = 0;
    I.G = 0;
    I.B = 0;
  }
  col->R = I.R;
  col->G = I.G;
  col->B = I.B;
  return;
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
 struct ray3D *ray;		// Structure to keep the ray from e to a pixel
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

 if (cam==NULL)
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
    // get camera coord of the pixel
    int dim;
    if(antialiasing){
      dim = 3; // dimension of antialiasing
    }
    else{
      dim = 1;
    }
    struct colourRGB totalcol;
    totalcol.R=0;    
    totalcol.G=0;   
    totalcol.B=0;
    for(int k = 0;k<dim;k++){
      for(int l = 0;l<dim;l++){
        pc.px = cam->wl + du * i + k*du/dim;
        pc.py = cam->wt + dv * j + l*dv/dim;
        pc.pz = cam-> f;
        pc.pw = 1;
        // convert to world coord
        matVecMult(cam->C2W, &pc);
        // get direction vector
        d.px = pc.px;
        d.py = pc.py;
        d.pz = pc.pz;
        d.pw = 1;
        subVectors(&e, &d);
        normalize(&d);
        // get ray
        ray = newRay(&pc, &d);
        ray->inside = 0;
        // cast ray
        struct object3D *Os = NULL;

        rayTrace(ray, 1, &col, Os);
        totalcol.R += col.R;
        totalcol.G+= col.G;
        totalcol.B += col.B;
        free(ray);
      }
    }
    totalcol.R = totalcol.R/(dim*dim);
    totalcol.G = totalcol.G/(dim*dim);
    totalcol.B = totalcol.B/(dim*dim);
    rgbIm[j * sx * 3 + i * 3 ] = (unsigned char)(totalcol.R * 255);
    rgbIm[j * sx * 3 + i * 3 + 1] = (unsigned char)(totalcol.G * 255);
    rgbIm[j * sx * 3 + i * 3 + 2] = (unsigned char)(totalcol.B * 255);
    
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

