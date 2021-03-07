/* =============== */
/* Header Includes */
/* =============== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

/* ===================== */
/* Constants Declaration */
/* ===================== */

/* Integration Method */

//#define EULER
#define VERLET
//#define LEAPFROG
//#define RUNGEKUTTA

/* Physical Constants */

#define NBODY 3
#define NDIM 3

#define G 1.0
#define dt 0.002

#define RMIN 0.05

/* Graphical Configuration */

#define WINDOWWIDTH 750
#define WINDOWHEIGHT 750

#define FPS 60.0
#define PATHSIZE 1000000

#define RADIANSPERDEGREE 0.0174533
#define SPHEREDENSITY 2000.0
#define ARROWSCALE 0.5
#define ARROWDIAMETER 0.01

#define NGRID 10
#define BOXLEN 4.0
#define BOXEDGE (float)BOXLEN/2.0

#define TOGGLEROTATION 1
#define TOGGLEPATH 2
#define TOGGLEMOMENTUM 3

#define SHOWFPS 1
#define SHOWTIME 1
#define SHOWENERGY 1
#define DRAWCENTEROFMASS 0

/* ========== */
/* Structures */
/* ========== */

struct body {
    float m;
    float x[NDIM];
    float p[NDIM];
};

/* ================ */
/* Global Variables */
/* ================ */

/* Initial Conditions */

/*
const float mass[NBODY] = {1.0,2.0};
const float x0[NBODY][NDIM] = {{1.0,0.0,0.0},{-1.0,0.0,0.0}};
const float p0[NBODY][NDIM] = {{0.0,0.5,0.0},{0.0,-0.4,0.0}};
*/

/*
const float mass[NBODY] = {1.0,1.0};
const float x0[NBODY][NDIM] = {{-1.0,0.2,0.0},{0.8,0.2,0.0}};
const float p0[NBODY][NDIM] = {{0.20,0.20,0.0},{-0.25,0.60,0.0}};
*/

/*
const float mass[NBODY] = {1.0,1.0};
const float x0[NBODY][NDIM] = {{0.0,0.0,0.0},{-1.0,0.0,0.0}};
const float p0[NBODY][NDIM] = {{0.0,0.0,0.0},{0.0,0.5,0.0}};
*/

/*
const float mass[NBODY] = {1.0,1.0};
const float x0[NBODY][NDIM] = {{0.0,0.0,0.0},{1.0,0.0,0.0}};
const float p0[NBODY][NDIM] = {{0.0,0.0,0.0},{0.2,1.0,0.0}};
*/

/*
const float mass[NBODY] = {3.0,1.0};
const float x0[NBODY][NDIM] = {{0.0,0.0,0.0},{2.0,0.5,0.0}};
const float p0[NBODY][NDIM] = {{-0.25,-0.25,0.0},{-0.75,0.60,0.0}};
*/

/* --------------- */
/* 3-Body Periodic */
/* --------------- */

/*
const float p1=0.464445;
const float p2=0.396060;
*/

const float p1=0.347111;
const float p2=0.532728;

const float mass[NBODY] = {1.0,1.0,1.0};
const float x0[NBODY][NDIM] = {{-1.0,0.0,0.0},{1.0,0.0,0.0},{0.0,0.0,0.0}};
const float p0[NBODY][NDIM] = {{p1,p2,0.0},{p1,p2,0.0},{-2.0*p1,-2.0*p2,0.0}};

/* pthread Variables */

pthread_t integrateThreadID;

/* Camera Configuration */

float rCam=1.0;
float thetaCam=M_PI_2*3.0/4.0;
float fiCam=-M_PI_4/2.0;

float rStepCam=0.2;
float angleStepCamActive=0.02f;
float angleStepCamPassive=0.01f;

/* Control Variables */

int rotationEnabled=0;
int pathEnabled=1;
int momentumEnabled=1;
int integrate=1;

float deltaRCam=0.0;
float deltaThetaCam=0.0;
float deltaFiCam=0.0;

int pathCurrent=0;

float stepsPerSecond=1000.0;

float deltaStepsPerSecondStep=1.0;

float deltaStepsPerSecond=0.0;

/* Data Variables */

float t;

struct body m[NBODY];

float path[NBODY][NDIM][PATHSIZE];

float initialEnergy;

/* Graphical Variables */

GLfloat light_position[] = {1.0,1.0,1.0,0.0};
GLfloat mat_specular[] = {1.0,1.0,1.0,1.0};
GLfloat mat_shininess[] = {50.0};
GLfloat qaWhite[] = {1.0, 1.0, 1.0, 1.0};
GLfloat qaBrown[] = {0.82, 0.70, 0.55, 1.0};
GLfloat qaBlack[] = {0.0, 0.0, 0.0, 1.0};
GLfloat qaYellow[] = {1.0, 1.0, 0.0, 1.0};
GLfloat qaOffWhite[] = {0.25, 0.25, 0.25, 1.0};
GLfloat qaOffRed[] = {0.25, 0.0, 0.0, 1.0};
GLfloat qaOffGreen[] = {0.0, 0.25, 0.0, 1.0};
GLfloat qaOffBlue[] = {0.0, 0.0, 0.25, 1.0};
GLfloat qaOffYellow[] = {0.25, 0.25, 0.0, 1.0};

/* =================== */
/* Function Prototypes */
/* =================== */

float r(float x, float y, float z);
float computeEnergy(void);

void printUsage(void);
void printStatus(int status);

void setInitialConditions(void);
void *integrateThread(void *arg);
void takeStep(void);
void redisplayTimer(int value);

void normalKey(unsigned char key, int x, int y);
void normalKeyRelease(unsigned char key, int x, int y);
void specialKey(int key, int x, int y);
void specialKeyRelease(int key, int x, int y);
void processMenuEvent(int option);

void glPrintString(char *string, float x, float y, float z);
void glArrow(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat D, GLfloat *qaColor);
void resize(int w, int h);
void display(void);

int main(int argc, char **argv);