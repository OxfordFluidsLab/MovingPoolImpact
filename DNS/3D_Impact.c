// Main developer: Radu Cimpeanu
// Contributor: Thomas Sykes
// Date: 24/10/2025

#include <stdlib.h>

#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include <sys/stat.h>                // required for mkdir
#include "grid/octree.h"
#include "navier-stokes/centered.h"  // solve NS equations
#include "two-phase.h"               // two fluid phases
#include "tension.h"                 // include surf tension between phases
#include "vof.h"                     // solve using VoF method
#include "fractions.h"               // initially define fraction of fluid phases
#include "view.h"                    // need to make the animations
#include "tag.h"                     // helps track droplet properties
#include "draw.h"                    // visualisation helper
#include "contact.h"                 // contact angle specification at boundaries
#include "tracer.h"                  // passive scalar tracer
#include "maxruntime.h"              // kill simulation after a specified time

/* Physical Parameters */
// density: kg m^-3
double rhol = 1089.0;
double rhoa = 1.2;

// dynamic viscosity: kg m^-1 s^-1
double mul = 3.0e-3;
double mua = 1.8e-5;

// surface tension: N m^-1
double sig = 70.3e-3;

// drop radius: m
double R0;

// pool depth: m
double Pool_depth;

// drop velocity: m s^-1
double Udrop;

// drop velocity: m s^-1
double Upool;

// acceleration due to gravity: m s^-2
double Gacc = 9.81;

/* Dimensionless Constants */
double Re; // Reynolds number
double Fr; // Froude number
double We; // Weber number

int maxlevel; //  = max refinement level
double impact_angle;

double t_end; // max simulation time
double domainSize;
double southPoleHeight = 0.1;

scalar omega[]; // vorticity
scalar velfield[]; // velocity field norm

vector h[];
double theta0 = 90; // specified (default) contact angle

face vector av[];

scalar drop_tracer[], pool_tracer[];
scalar * tracers = {drop_tracer, pool_tracer};

// output files
FILE * fp_interface;
FILE * fp_stats;

/* ====== */
/* Set Up */
/* ====== */
int main(int argc, char * argv[]) {

  // Function to call for maximum runtime, which is an optional arguement (-m)
  maxruntime(&argc, argv);

  maxlevel = atoi(argv[1]);
  impact_angle = atof(argv[2]);
  Udrop = atof(argv[3]);
  Upool = atof(argv[4]);
  R0 = atof(argv[5]);
  Pool_depth = atof(argv[6]);
  domainSize = atof(argv[7]);
  t_end = atof(argv[8]);
  
  Re = (rhol * Udrop * R0 / mul);  
  Fr = (Udrop / sqrt(Gacc * R0));
  We = (rhol * Udrop * Udrop * R0 / sig);

  size(domainSize); 
  origin(0.0, 0.0, -(domainSize/2.0));
 
   // Make folders to store the outputs
   mkdir("Slices", 0700);
   mkdir("Animations", 0700);
   mkdir("Interfaces", 0700);
 
  //periodic (front); 

  a = av;
  
  f.sigma = 1./We; // surface tension
  f.height = h;
  
  mu1 = 1./Re;         // fluid viscosity
  mu2 = mu1/(mul/mua); // gas viscosity
  
  rho1 = 1.;              // fluid density
  rho2 = 1./(rhol/rhoa);   // gas density

  init_grid(128);

    // Pointers to files
    {
      char name[200];
      sprintf(name, "loginterface.dat");
      fp_interface = fopen(name, "a");
    }

    {
      char name[200];
      sprintf(name, "logstats.dat");
      fp_stats = fopen(name, "a");
    }

    DT = 1.0e-3;
    NITERMIN = 1;     // default 1
    NITERMAX = 200;   // default 100
    TOLERANCE = 1e-4; // default 1e-3
 
    run(); // Runs the simulation

    // Close files
    fclose(fp_interface);
    fclose(fp_stats);
}


/* =================== */
/* Boundary Conditions */
/* =================== */

// outflow at top
u.n[top] = neumann(0);
p[top] = dirichlet(0);
pf[top] = dirichlet(0);

// inflow at front
u.n[front] = dirichlet(f[]*((Upool/Udrop)));

// outflow at back
u.n[back] = neumann(0);
p[back] = neumann(0);
pf[back] = neumann(0);

// gravity inclusion via Froude number
event acceleration (i++) {
  foreach_face(x)  
    av.x[] += 0.0;
  foreach_face(y)  
    av.y[] -= 1./(Fr*Fr);
  foreach_face(z)  
    av.z[] += 0.0;
}

/* ================== */
/* Initial Conditions */
/* ================== */
event init(t=0.) {

if (!restore(file = "restart"))
    {

  // Double-check dimensionless groupings
  fprintf(stdout, "Reynolds number Re = %0.6f \n", Re);
  fflush(stdout);
  fprintf(stdout, "Weber number We = %0.6f \n", We);
  fflush(stdout);
  fprintf(stdout, "Froude number Fr = %0.6f \n", Fr);
  fflush(stdout);
  fprintf(stdout, "Density ratio = %0.6f \n", (rhol/rhoa));
  fflush(stdout);
  fprintf(stdout, "Viscosity ratio = %0.6f \n", (mul/mua));
  fflush(stdout);

  double dropr = 1.;

  double dropx = 0.;
  double dropy = (Pool_depth/R0) + dropr + southPoleHeight;
  double dropz = 0.;
  if (impact_angle < 90.0) // check for normal impact
  	dropz = -0.016666*impact_angle + 1.5; // general formula to maintain the south pole at the given height for non-normal impact, linearly from z = 1.25 for angle=15 to z = 0 for angle=90

  refine((((sq(x - dropx) + sq(y - dropy) + sq(z - dropz) < sq(dropr*1.025)) && \
          ((sq(x - dropx) + sq(y - dropy) + sq(z - dropz) > sq(dropr*0.975)))))
          && level < 9);
  refine(((y > ((Pool_depth/R0) - 0.025)) && (y < ((Pool_depth/R0) + 0.025))) && level < 9);
  fraction(f, union((Pool_depth/R0) - y, sq(dropr) - sq(x - dropx) - sq(y - dropy) - sq(z - dropz))); // The drop and pool

  fraction(drop_tracer, sq(dropr) - sq(x - dropx) - sq(y - dropy) - sq(z - dropz)); // Drop tracer
  fraction(pool_tracer, (Pool_depth/R0) - y); // Pool tracer

  DT = 1.0e-3;
  NITERMIN = 1;     // default 1
  NITERMAX = 200;   // default 100
  TOLERANCE = 1e-4; // default 1e-3

  /* Set initial velocity */
  foreach() {
    // if within the drop (and slightly outside for smoothness) set relevant velocities
    if (sq(x - dropx) + sq(y - dropy) + sq(z - dropz) < 1.05*sq(dropr)) {
	u.x[] =  0.0;      
	u.y[] = -sin(M_PI*(impact_angle/180.0));
	u.z[] = -cos(M_PI*(impact_angle/180.0)); 
    } else {
	u.x[] =  0.0;      
	u.y[] =  0.0;
	u.z[] =  (Upool/Udrop)*f[]; 
    }
  }
}
}


/* ========== */
/* Adapt Grid */
/* ========== */
event adapt (i++) {
  /* Compute vorticity */
  vorticity(u, omega);

  foreach()
    velfield[] = pow(u.x[]*u.x[]+u.y[]*u.y[]+u.z[]*u.z[], 0.5);

  /* ======== */
  /* Refining */
  /* ======== */
  // Refine based on interface location, velocity and vorticity
  adapt_wavelet((scalar*){f,drop_tracer,u}, (double[]){1e-4, 1e-2, 1e-2, 1e-2, 1e-2}, maxlevel, maxlevel-4);

  /* Unrefine far away from the impact site */
  unrefine( (sq(x) + sq(z) > sq(2.0)) && level > (maxlevel-2));

}

// Output interface statistics (for mass conservation check as well)
event loginterface (t = 0.0; t += 0.01) {

    scalar posX[],posY[],posZ[];
    position (f, posX, {1,0,0}); // (1,0,0) indicates the unit vector in the x-direction
    position (f, posY, {0,1,0}); // (0,1,0) indicates the unit vector in the y-direction
    position (f, posZ, {0,0,1}); // (0,1,0) indicates the unit vector in the y-direction

    fprintf(fp_interface, "%i %g %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f \n", i, t, statsf(f).sum, statsf(posX).min, statsf(posX).max, statsf(posY).min, statsf(posY).max, statsf(posZ).min, statsf(posZ).max);
    fflush(fp_interface);
}


// Output stats
event logstats (t = 0.0; t += 0.001; t <= t_end) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}

// Save regular simulation slices
event gfsview (t += 0.1) {
    char name_gfs[200];
    sprintf(name_gfs,"Slices/DropImpact-%0.2f.gfs",t);

    FILE* fp_gfs = fopen (name_gfs, "w");
    output_gfs(fp_gfs);
    fclose(fp_gfs);  
}

event small_droplet_removal (i++) {
// Removes any small droplets that have formed, that are smaller than a specific size
    remove_droplets(f, 8);       // Removes droplets of diameter 8 cells or less
    remove_droplets(f, 8, true); // Removes bubbles of diameter 8 cells or less
}

// Raw output of fluid-fluid interfaces
event interfaceShapeFinder (t += 0.01) {

  char nameInterfacesLiquidRaw[200];
  char nameInterfacesDropRaw[200];
  char nameInterfacesPoolRaw[200];

  sprintf(nameInterfacesLiquidRaw,"Interfaces/interfacesLiquidRaw-%0.3f.dat",t);
  sprintf(nameInterfacesDropRaw,"Interfaces/interfacesDropRaw-%0.3f.dat",t);
  sprintf(nameInterfacesPoolRaw,"Interfaces/interfacesPoolRaw-%0.3f.dat",t);

  FILE *fpLiquidRaw = fopen(nameInterfacesLiquidRaw, "w");
  FILE *fpDropRaw = fopen(nameInterfacesDropRaw, "w");
  FILE *fpPoolRaw = fopen(nameInterfacesPoolRaw, "w");

  scalar ffLiquid[], ffDrop[], ffPool[];

  foreach() {
    if (x < 1e-2){ // 1e-2 for level 9, 5e-3 for level 10
      ffLiquid[] = f[] < 1.0e-6 ? 0 : f[] > 1. - 1.0e-6 ? 1. : f[];
      ffDrop[] = drop_tracer[] < 0.4 ? 0 : drop_tracer[] > 0.6 ? 1. : drop_tracer[];
      ffPool[] = pool_tracer[] < 0.4 ? 0 : pool_tracer[] > 0.6 ? 1. : pool_tracer[];
    } else {
      ffDrop[] = 1.0;
      ffPool[] = 1.0;
    }
  }

  output_facets (ffLiquid, fpLiquidRaw);
  output_facets (ffDrop, fpDropRaw);
  output_facets (ffPool, fpPoolRaw);

  fclose(fpLiquidRaw);
  fclose(fpDropRaw);
  fclose(fpPoolRaw);
}

// Output fluid-fluid interfaces
event saveInterfaces (t += 0.01) {

    char nameInterfaces1[200];
    char nameInterfacesDrop[200];
    char nameInterfacesPool[200];

    sprintf(nameInterfaces1,"Interfaces/interfacesLiquid-%0.1f.dat",t);
    sprintf(nameInterfacesDrop,"Interfaces/interfacesDrop-%0.1f.dat",t);
    sprintf(nameInterfacesPool,"Interfaces/interfacesPool-%0.1f.dat",t);

    FILE * fp1 = fopen(nameInterfaces1, "w");
    FILE * fpDrop = fopen(nameInterfacesDrop, "w");
    FILE * fpPool = fopen(nameInterfacesPool, "w");

    output_facets (f, fp1);	
    output_facets (drop_tracer, fpDrop);
    output_facets (pool_tracer, fpPool);

    fclose(fp1);
    fclose(fpDrop);
    fclose(fpPool);
}


// Create flow animations (using a few helper variables)

scalar velnorm[], liquids[];

event movies (t += 0.001)
{

  scalar omega[];
  vorticity (u, omega);
  
    foreach(){
  	velnorm[] = sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
  	liquids[] = 1-f[]+drop_tracer[]/2.;
  }

  view (fov = 20.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  //cells(n = {1,0,0});
  squares ("u.x", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  save ("Animations/Vel_Ux.mp4");

  view (fov = 20.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  //cells(n = {1,0,0});
  squares ("u.y", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  save ("Animations/Vel_Uy.mp4");
  
  view (fov = 20.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  //cells(n = {1,0,0});
  squares ("u.z", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  save ("Animations/Vel_Uz.mp4");  
  
  view (fov = 20.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  cells(n = {1,0,0});
  squares ("liquids", n = {1,0,0}, map = cool_warm, min = 0.0, max = 2.0);
  save ("Animations/LiquidsGrid.mp4");  
  
  view (fov = 20.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  //cells(n = {1,0,0});
  squares ("liquids", n = {1,0,0}, map = cool_warm, min = 0.0, max = 2.0);
  save ("Animations/Liquids.mp4");  

  view (fov = 20.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  //cells(n = {1,0,0});
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  save ("Animations/Velocity.mp4");  
  
  view (fov = 20.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  //cells(n = {1,0,0});
  squares ("omega", n = {1,0,0}, map = cool_warm, min = -2.5, max = 2.5);
  save ("Animations/Vorticity.mp4");  

  view (fov = 20.0, camera = "left", tx = 0.0, ty = -0.5, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  //cells(n = {1,0,0});
  squares ("p", n = {1,0,0}, map = cool_warm, min = -0.3, max = 0.6);
  save ("Animations/Pressure.mp4");
  
  view (fov = 20.0, camera = "right", tx = 0.0, ty = 0.0, bg = {1,1,1},
	  width = 900, height = 900);
  clear();
  squares ("velnorm", spread = -1, linear = true, map = cool_warm, n = {1,0,0});
  cells(n = {1,0,0});
  isosurface("f", 0.5, fc = {1.0,1.0,1.0});
  save ("Animations/Velocity_Front_All.mp4");

}

