#include "simpleRayTracer.h"
#include <mpi.h>

// to compile:
// gcc -O3 -o simpleRayTracer *.c -I.  -fopenmp -lm

// to run:
//  ./simpleRayTracer

// to compile animation:
//   ffmpeg -y -i image_%05d.ppm -pix_fmt yuv420p foo.mp4

int main(int argc, char *argv[]){

  int rank, size;

  // initialize MPI backend infrastructure
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("MPI process from %d of %d ranks\n", rank, size);
    
  initTimer();
  
  // initialize triangles and spheres
  scene_t *scene = sceneSetup();

  grid_t     *grid      = scene->grid;
  shape_t    *shapes    = scene->shapes;
  material_t *materials = scene->materials;
  light_t    *lights    = scene->lights;
  
  /* Will contain the raw image */
  unsigned char *img = (unsigned char*) calloc(3*WIDTH*HEIGHT, sizeof(char));

  // 1. location of observer eye (before rotation)
  sensor_t sensor;

  // background color
  sensor.bg.red   = 126./256;
  sensor.bg.green = 192./256;
  sensor.bg.blue  = 238./256;

  dfloat br = 3.75f*BOXSIZE;

  // angle elevation to y-z plane
  dfloat eyeAngle = M_PI/4.f; // 0 is above, pi/2 is from side.  M_PI/3; 0; M_PI/2.;

  // target view
  vector_t targetX = vectorCreate(BOXSIZE/2., BOXSIZE, BOXSIZE/2.); // this I do not understand why target -B/2
  sensor.eyeX = vectorAdd(targetX, vectorCreate(0, -br*cos(eyeAngle), -br*sin(eyeAngle))); 
  dfloat sensorAngle = eyeAngle; +15.*M_PI/180.;
  sensor.Idir   = vectorCreate(1.f, 0.f, 0.f);
  sensor.Jdir   = vectorCreate(0.f, sin(sensorAngle), -cos(sensorAngle));
  vector_t sensorNormal = vectorCrossProduct(sensor.Idir, sensor.Jdir);
  
  // 2.4 length of sensor in axis 1 & 2
  sensor.Ilength = 20.f;
  sensor.Jlength = HEIGHT*20.f/WIDTH;
  sensor.offset  = 0.f;

  // 2.5 normal distance from sensor to focal plane
  dfloat lensOffset = 50;
  sensor.lensC = vectorAdd(sensor.eyeX, vectorScale(lensOffset, vectorCrossProduct(sensor.Idir, sensor.Jdir)));

  // why 0.25 ?
  sensor.focalPlaneOffset = 0.22f*fabs(vectorTripleProduct(sensor.Idir, sensor.Jdir, vectorSub(targetX,sensor.eyeX))); // triple product
  
  //  sensor.focalOffset = 0.8*BOXSIZE - sensor.lensC.z; // needs to be distance to plane from sensor

  printf("lensOffset = %g, sensor.focalPlaneOffset = %g\n", lensOffset, sensor.focalPlaneOffset);
  
  dfloat *randomNumbers = (dfloat*) calloc(2*NRANDOM, sizeof(dfloat));
  for(int i=0;i<NRANDOM;++i){
    dfloat r1 = 2*drand48()-1;
    dfloat r2 = 2*drand48()-1;

    randomNumbers[2*i+0] = r1/sqrt(r1*r1+r2*r2);
    randomNumbers[2*i+1] = r2/sqrt(r1*r1+r2*r2);
  }
  
  // number of angles to render at
  int Ntheta = 20; // instead of 240;
  
  // loop over scene angles
  for(int thetaId=0;thetaId<Ntheta;++thetaId){
    
    /* rotation angle in y-z */
    dfloat theta = thetaId*M_PI*2./(dfloat)(Ntheta-1);

    /* sort objects into grid */
    gridPopulate(grid, scene->Nshapes, shapes);

    /* for Question 7, Start timer now replaced with MPI_Wtime() */
    // ticTimer();
    double tic = MPI_Wtime();
    
    /* render scene */
    renderKernel(WIDTH,
		 HEIGHT,
		 scene[0],
		 sensor,
		 cos(theta), 
		 sin(theta),
		 randomNumbers,
		 img);
 
    // MPI modification for Q6 starts here:
    // number of elements in send and recv buffers,
    // shared among the size of ranks:
    int messageN = 3*WIDTH*HEIGHT/size;
    int messageTag = 999; // message tag(integer value)
    int messageDest = 0; // rank of destination(integer)
    
    // pointers to the initial address of send and recv buffers:
    unsigned char *sendBuf = (unsigned char*) calloc(3*WIDTH*HEIGHT, sizeof(char));
    unsigned char *recvBuf = (unsigned char*) calloc(3*WIDTH*HEIGHT, sizeof(char));

    if (rank != 0){
      sendBuf[rank] = *img;
      // sends all pixels to the img array in rank 0:
      MPI_Send(sendBuf,
	       messageN,
	       MPI_UNSIGNED_CHAR,
	       messageDest,
	       messageTag,
	       MPI_COMM_WORLD);

      // debugging assistance:
      printf("rank %d sends pixels to rank %d\n", rank, messageDest);
    }

    // rank 0 receives
    if (rank == 0){
      // rank 0 will now properly align all image sections to their respectively position
      for (int i = 1; i< size; ++i){
	MPI_Status status;
	int messageSource = i;
	recvBuf[i] = *img;
	MPI_Recv(recvBuf+i*messageN,
		 messageN,
		 MPI_UNSIGNED_CHAR,
		 i,
		 messageTag,
		 MPI_COMM_WORLD,
		 &status);

	// debugging assistance:
	printf("rank %d receives pixels from rank %d\n", rank, i);
      }
    }
    /*
    // -------------------------------------------------------------
    // for Q6 Extra credit: using MPI_Gather() in place of MPI_Send+MPI_Recv
    // allocate enough space to contain final rendered image - from image fragments on each rank
    unsigned char *img_complete = (unsigned char*) calloc(3*WIDTH*HEIGHT, sizeof(char));
    // number of elements shared among the size of ranks:
    int messageN = 3*WIDTH*HEIGHT/size;
    int root = 0;
    
    // Gathers together values from a group of processes
    // create a complete image by adding segments behind receiving process
    MPI_Gather(img+(3*WIDTH*HEIGHT*(rank)/size),
	       messageN,
	       MPI_UNSIGNED_CHAR,
	       img_complete,
	       messageN,
	       MPI_UNSIGNED_CHAR,
	       root,
	       MPI_COMM_WORLD);
    
    printf("image makeup complete.\n");
    
    // -------------------------------------------------------------
    */
    
    /* for Q7, report elapsed time now replaced with MPI_Wtime() */
    // tocTimer("recursiveRenderKernel");
    double toc = MPI_Wtime();

    /* to calculate total elapsed time spent by rank 0 generating each frame: */
    // will result ending on rank = messageDest = 0
    double time_elapsed = toc - tic;
    int messageLength = 1;

    double *messageOut = (double*) calloc(messageLength, sizeof(double));
    double *messageIn = (double*) calloc(messageLength, sizeof(double));

    messageOut[0] = time_elapsed;

    MPI_Reduce(messageOut,
	       messageIn,
	       messageLength,
	       MPI_DOUBLE,
	       MPI_SUM,
	       messageDest,
	       MPI_COMM_WORLD);
    
    // only rank == messageDest finalizes the sum
    if (rank==messageDest){
      printf("Total elapsed time by rank 0 generating each frame = %lg\n", messageIn[0]);
      printf("Average time to generate each fram = %lg\n", messageIn[0]/size);
    }

    
    dfloat dt = .025, g = 1;
    int NsubSteps= 40;

    ticTimer();
    
    // collide and move spheres in time and update grid
    for(int subStep=0;subStep<NsubSteps;++subStep){
      
      sphereCollisions(grid, dt, g, scene->Nshapes, shapes);

      sphereUpdates(grid, dt, g, scene->Nshapes, shapes);

      gridPopulate(grid, scene->Nshapes, shapes);
    }

    // report time taken to move and collide spheres
    tocTimer("move and collide Spheres");

    /* save scene as ppm file */
    char fileName[BUFSIZ];

    // make sure images directory exists
    mkdir("images", S_IRUSR | S_IREAD | S_IWUSR | S_IWRITE | S_IXUSR | S_IEXEC);

    // write image as ppm format file
    if (rank == 0){
      sprintf(fileName, "images/image_%05d.ppm", thetaId);
      saveppm(fileName, img, WIDTH, HEIGHT);
    }
  }
  
  free(img);

  MPI_Finalize();
  
  return 0;
}
