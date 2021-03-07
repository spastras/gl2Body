#include "gl2Body.h"

/* r */

float r(float x, float y, float z) {
    return sqrt(x*x+y*y+z*z);
}

/* computeEnergy */

float computeEnergy(void) {
    int i, j;
    float energy=0.0;
    for(i=0;i<NBODY;i++) {
        for(j=0;j<NDIM;j++) {
            energy=energy+0.5*m[i].p[j]*m[i].p[j]/m[i].m;
        }
    }
    for(i=0;i<(NBODY-1);i++) {
        for(j=i+1;j<NBODY;j++) {
            energy=energy-G*m[i].m*m[j].m/r(m[i].x[0]-m[j].x[0],m[i].x[1]-m[j].x[1],m[i].x[2]-m[j].x[2]);
        }
    }
    return energy;
}

/* printUsage */

void printUsage(void) {
    printf("\n");
    printf("[i] -<=============>-\n");
    printf("[i] -<=[ gl2Body ]=>-\n");
    printf("[i] -<=============>-\n");
    printf("\n");
    printf("[i] Controls:\n");
    printf("\n");
    printf("[i] (arrows)\t->\trotate\n");
    printf("[i] r\t\t->\tenable\\disable rotation\n");
    printf("[i] ALT + r\t->\tswitch direction of rotation\n");
    printf("[i] p\t\t->\tplay\\pause\n");
    printf("[i] ALT + p\t->\treset\n");
    printf("[i] + \\ -\t->\tincrease\\decrease steps per second\n");
    printf("[i] q\t\t->\tquit\n");
    printf("\n");
    return;
}

/* printStatus */

void printStatus(int status) {
    printf("\e[2K\r[i] Status: ");
    switch (status) {
        case 0 :
            printf("paused");
            break;
		case 1 :
            printf("running");
            break;
		case -1 :
            printf("initializing");
            break;
		case -2 :
            printf("exiting\n");
            break;
	}
    fflush(stdout);
    return;
}

/* setInitialConditions */

void setInitialConditions(void) {
    int i, j;
    for(i=0;i<NBODY;i++) {
        m[i].m=mass[i];
        for(j=0;j<NDIM;j++) {
            m[i].x[j]=x0[i][j];
            m[i].p[j]=p0[i][j];
        }
    }
    t=0.0;
    return;
}

/* integrateThread */

void *integrateThread(void *arg) {
    int i, j, rlimitReached=0;
    printStatus(1);
    while(integrate) {
        for(i=0;i<(NBODY-1);i++) {
            for(j=i+1;j<NBODY;j++) {
                if(r(m[i].x[0]-m[j].x[0],m[i].x[1]-m[j].x[1],m[i].x[2]-m[j].x[2])<RMIN) {
                    rlimitReached=1;
                }
            }
        }
        if(rlimitReached) {
            break;
        }
        stepsPerSecond+=deltaStepsPerSecond;
        usleep((useconds_t)1000000.0/stepsPerSecond);
        takeStep();
    }
    if(integrate) {
        integrate=!integrate;
    }
    printStatus(0);
    return NULL;
}

/* takeStep */

void takeStep(void) {

    int i, j, k;
    float x[NBODY][NDIM], p[NBODY][NDIM];
    float rr, dxx, h;

    #ifdef RUNGEKUTTA

    float xx[NBODY][NDIM], pp[NBODY][NDIM];
    float adx[NBODY][NDIM], adp[NBODY][NDIM];
    float bdx[NBODY][NDIM], bdp[NBODY][NDIM];
    float cdx[NBODY][NDIM], cdp[NBODY][NDIM];
    float ddx[NBODY][NDIM], ddp[NBODY][NDIM];
    float dxdt[NBODY][NDIM], dpdt[NBODY][NDIM];

    #endif

    #ifdef LEAPFROG

    float xmid[NBODY][NDIM], pdot[NBODY][NDIM];

    #endif

    #ifdef RUNGEKUTTA

    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            adx[i][k]=m[i].p[k]/m[i].m;
        }
    }
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            adp[i][k]=0.0;
        }
        for(j=0;j<NBODY;j++) {
            if(i==j) continue;
            h=G*m[i].m*m[j].m;
            rr=r(m[i].x[0]-m[j].x[0],m[i].x[1]-m[j].x[1],m[i].x[2]-m[j].x[2]);
            for(k=0;k<NDIM;k++) {
                dxx=m[i].x[k]-m[j].x[k];
                adp[i][k]=adp[i][k]-(h*dxx/(rr*rr*rr));
            }
        }
    }

    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            xx[i][k]=m[i].x[k]+adx[i][k]*0.5*dt;
            pp[i][k]=m[i].p[k]+adp[i][k]*0.5*dt;
        }
    }
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            bdx[i][k]=pp[i][k]/m[i].m;
        }
    }
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            bdp[i][k]=0.0;
        }
        for(j=0;j<NBODY;j++) {
            if(i==j) continue;
            h=G*m[i].m*m[j].m;
            rr=r(xx[i][0]-xx[j][0],xx[i][1]-xx[j][1],xx[i][2]-xx[j][2]);
            for(k=0;k<NDIM;k++) {
                dxx=xx[i][k]-xx[j][k];
                bdp[i][k]=bdp[i][k]-(h*dxx/(rr*rr*rr));
            }
        }
    }
    
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            xx[i][k]=m[i].x[k]+bdx[i][k]*0.5*dt;
            pp[i][k]=m[i].p[k]+bdp[i][k]*0.5*dt;
        }
    }
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            cdx[i][k]=pp[i][k]/m[i].m;
        }
    }
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            cdp[i][k]=0.0;
        }
        for(j=0;j<NBODY;j++) {
            if(i==j) continue;
            h=G*m[i].m*m[j].m;
            rr=r(xx[i][0]-xx[j][0],xx[i][1]-xx[j][1],xx[i][2]-xx[j][2]);
            for(k=0;k<NDIM;k++) {
                dxx=xx[i][k]-xx[j][k];
                cdp[i][k]=cdp[i][k]-(h*dxx/(rr*rr*rr));
            }
        }
    }
    
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            xx[i][k]=m[i].x[k]+cdx[i][k]*dt;
            pp[i][k]=m[i].p[k]+cdp[i][k]*dt;
        }
    }
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            ddx[i][k]=pp[i][k]/m[i].m;
        }
    }
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            ddp[i][k]=0.0;
        }
        for(j=0;j<NBODY;j++) {
            if(i==j) continue;
            h=G*m[i].m*m[j].m;
            rr=r(xx[i][0]-xx[j][0],xx[i][1]-xx[j][1],xx[i][2]-xx[j][2]);
            for(k=0;k<NDIM;k++) {
                dxx=xx[i][k]-xx[j][k];
                ddp[i][k]=ddp[i][k]-(h*dxx/(rr*rr*rr));
            }
        }
    }
    
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            dxdt[i][k]=(1.0/6.0)*(adx[i][k]+2.0*(bdx[i][k]+cdx[i][k])+ddx[i][k]);
            dpdt[i][k]=(1.0/6.0)*(adp[i][k]+2.0*(bdp[i][k]+cdp[i][k])+ddp[i][k]);
        }
    }
    
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            x[i][k]=m[i].x[k]+dxdt[i][k]*dt;
            p[i][k]=m[i].p[k]+dpdt[i][k]*dt;
        }
    }

    #endif

    #ifdef LEAPFROG

    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            xmid[i][k]=m[i].x[k]+m[i].p[k]*0.5*dt/m[i].m;
        }
    }
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            pdot[i][k]=0.0;
        }
        for(j=0;j<NBODY;j++) {
            if(i==j) continue;
            h=G*m[i].m*m[j].m;
            rr=r(xmid[i][0]-xmid[j][0],xmid[i][1]-xmid[j][1],xmid[i][2]-xmid[j][2]);
            for(k=0;k<NDIM;k++) {
                dxx=xmid[i][k]-xmid[j][k];
                pdot[i][k]=pdot[i][k]-(h*dxx/(rr*rr*rr));
            }
        }
    }
    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            p[i][k]=m[i].p[k]+pdot[i][k]*dt;
            x[i][k]=xmid[i][k]+p[i][k]*0.5*dt/m[i].m;
        }
    }

    #endif

    #ifdef VERLET

    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            p[i][k]=m[i].p[k];
        }
        for(j=0;j<NBODY;j++) {
            if(i==j) continue;
            h=G*m[i].m*m[j].m;
            rr=r(m[i].x[0]-m[j].x[0],m[i].x[1]-m[j].x[1],m[i].x[2]-m[j].x[2]);
            for(k=0;k<NDIM;k++) {
                dxx=m[i].x[k]-m[j].x[k];
                p[i][k]=p[i][k]-dt*(h*dxx/(rr*rr*rr));
            }
        }
        for(k=0;k<NDIM;k++) {
            x[i][k]=m[i].x[k]+dt*p[i][k]/m[i].m;
        }
    }

    #endif

    #ifdef EULER

    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            p[i][k]=m[i].p[k];
        }
        for(j=0;j<NBODY;j++) {
            if(i==j) continue;
            h=G*m[i].m*m[j].m;
            rr=r(m[i].x[0]-m[j].x[0],m[i].x[1]-m[j].x[1],m[i].x[2]-m[j].x[2]);
            for(k=0;k<NDIM;k++) {
                dxx=m[i].x[k]-m[j].x[k];
                p[i][k]=p[i][k]-dt*(h*dxx/(rr*rr*rr));
            }
        }
        for(k=0;k<NDIM;k++) {
            x[i][k]=m[i].x[k]+dt*m[i].p[k]/m[i].m;
        }
    }
    
    #endif

    for(i=0;i<NBODY;i++) {
        for(k=0;k<NDIM;k++) {
            m[i].x[k]=x[i][k];
            m[i].p[k]=p[i][k];
        }
    }

    t=t+dt;

    if(pathEnabled) {

        /* Save path point */

        if(pathCurrent==PATHSIZE) {
            pathCurrent=0;
        }

        for(i=0;i<NBODY;i++) {
            for(k=0;k<NDIM;k++) {
                path[i][k][pathCurrent]=m[i].x[k];
            }
        }
        pathCurrent++;

    }

    //printf("%lf %lf %lf\n",m[0].x[0],m[0].x[1],m[0].x[2]);

    return;

}

/* redisplayTimer */

void redisplayTimer(int value) {
    glutPostRedisplay();
    glutTimerFunc(1000.0/FPS, redisplayTimer, 0);
    return;
}

/* normalKey */

void normalKey(unsigned char key, int x, int y) {
    int modifiers;
    switch (key) {
		case 27 :
		case 113 :
            printStatus(-2);
            exit(0);
            break;
		case 114 :
            modifiers=glutGetModifiers();
            if(modifiers==GLUT_ACTIVE_ALT) {
                angleStepCamPassive=-angleStepCamPassive;
            }
            else {
                rotationEnabled=!rotationEnabled;
            }
            break;
		case 112 :
            modifiers=glutGetModifiers();
            if(modifiers==GLUT_ACTIVE_ALT) {
                pathCurrent=0;
                setInitialConditions();
            }
            else {
                integrate=!integrate;
                if(integrate) {
                    pthread_create(&integrateThreadID, NULL, integrateThread, NULL);
                }
                else {
                    pthread_join(integrateThreadID, NULL);
                }
            }
            break;
        case 43 :
            deltaStepsPerSecond=deltaStepsPerSecondStep;
            break;
        case 45 :
            deltaStepsPerSecond=-deltaStepsPerSecondStep;
            break;
        case 119 :
            //deltaRCam+=rStepCam;
            break;
        case 115 :
            //deltaRCam-=rStepCam;
            break;
	}
    return;
}

/* normalKeyRelease */

void normalKeyRelease(unsigned char key, int x, int y) {
	switch (key) {
        case 43 :
		case 45 :
            deltaStepsPerSecond=0.0;
            break;
		case 119 :
		case 115 :
            //deltaRCam = 0.0;
            break;
	}
    glutPostRedisplay();
    return;
}

/* specialKey */

void specialKey(int key, int x, int y) {
	switch (key) {
		case GLUT_KEY_LEFT :
			deltaFiCam-=angleStepCamActive;
			break;
		case GLUT_KEY_RIGHT :
			deltaFiCam+=angleStepCamActive;
			break;
		case GLUT_KEY_UP :
			deltaThetaCam-=angleStepCamActive;
			break;
		case GLUT_KEY_DOWN :
			deltaThetaCam+=angleStepCamActive;
			break;
	}
    glutPostRedisplay();
    return;
}

/* specialKeyRelease */

void specialKeyRelease(int key, int x, int y) {
	switch (key) {
		case GLUT_KEY_LEFT :
		case GLUT_KEY_RIGHT :
            deltaFiCam = 0.0;
            break;
		case GLUT_KEY_UP :
		case GLUT_KEY_DOWN :
            deltaThetaCam = 0.0;
            break;
	}
    glutPostRedisplay();
    return;
}

/* processMenuEvent */

void processMenuEvent(int option) {
	switch (option) {
		case TOGGLEROTATION :
			rotationEnabled=!rotationEnabled;
            break;
        case TOGGLEPATH :
			pathEnabled=!pathEnabled;
            pathCurrent=0;
            break;
        case TOGGLEMOMENTUM :
			momentumEnabled=!momentumEnabled;
            break;
	}
    return;
}

/* glPrintString */

void glPrintString(char *string, float x, float y, float z) {

    int i;

    glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT), 0);
	glMatrixMode(GL_MODELVIEW);

	glPushMatrix();
	glLoadIdentity();
    glRasterPos3f(x, y, z);
	for(i=0;string[i]!='\0';i++) {
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, string[i]);
	}
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

    return;

}

/* glArrow */

void glArrow(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat D, GLfloat *qaColor)
{

    float x=x2-x1;
    float y=y2-y1;
    float z=z2-z1;
    float L=sqrt(x*x+y*y+z*z);
    
    GLUquadricObj *quadObj;
    
    glPushMatrix();
    
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, qaColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, qaColor);
    
    glTranslated(x1,y1,z1);
    
    if((x!=0.0)||(y!=0.0)) {
        glRotatef(atan2(y,x)/RADIANSPERDEGREE,0.0,0.0,1.0);
        glRotatef(atan2(sqrt(x*x+y*y),z)/RADIANSPERDEGREE,0.0,1.0,0.0);
    } else if (z<0) {
        glRotatef(180.0,1.0,0.0,0.0);
    }

    glTranslatef(0.0,0.0,L-4.0*D);

    quadObj = gluNewQuadric();
    gluQuadricDrawStyle(quadObj, GLU_FILL);
    gluQuadricNormals(quadObj, GLU_SMOOTH);
    gluCylinder(quadObj, 2.0*D, 0.0, 4.0*D, 32, 1);
    gluDeleteQuadric(quadObj);

    quadObj = gluNewQuadric();
    gluQuadricDrawStyle(quadObj, GLU_FILL);
    gluQuadricNormals(quadObj, GLU_SMOOTH);
    gluDisk(quadObj, 0.0, 2.0*D, 32, 1);
    gluDeleteQuadric(quadObj);
    
    glTranslatef(0.0,0.0,-L+4.0*D);

    quadObj = gluNewQuadric();
    gluQuadricDrawStyle(quadObj, GLU_FILL);
    gluQuadricNormals(quadObj, GLU_SMOOTH);
    gluCylinder(quadObj, D, D, L-4.0*D, 32, 1);
    gluDeleteQuadric(quadObj);

    quadObj = gluNewQuadric();
    gluQuadricDrawStyle(quadObj, GLU_FILL);
    gluQuadricNormals(quadObj, GLU_SMOOTH);
    gluDisk(quadObj, 0.0, D, 32, 1);
    gluDeleteQuadric(quadObj);

    glPopMatrix();

    return;

}

/* resize */

void resize(int w, int h)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0,0,w,h);
    if(w<=h)
        glOrtho(-2.0,2.0,-2.0*h/w,2.0*h/w,-10.0,100.0);
    else
        glOrtho(-2.0*w/h,2.0*w/h,-2.0,2.0,-10.0,100.0);
    //gluPerspective(120,(float)h/(float)w,1,100);
    glMatrixMode(GL_MODELVIEW);
    return;
}

/* display */

void display(void)
{

    int i, j;
    #if SHOWFPS||SHOWTIME||SHOWENERGY
    char status[128]="";
    #endif
    #if SHOWFPS
    static int frameNumber;
    static long time, timebase;
    static float currentFPS;
    #endif
    #if SHOWENERGY
    float energy;
    #endif
    #if DRAWCENTEROFMASS
    float mtotal, xcm[NDIM];
    #endif

    /* Rotate the camera */

    if(rotationEnabled) {
        fiCam+=angleStepCamPassive;
    }
    /*
    if(deltaRCam) {
        rCam+=deltaRCam;
        glutPostRedisplay();
    }
    */
    if(deltaThetaCam) {
        thetaCam+=deltaThetaCam;
        glutPostRedisplay();
    }
    if(deltaFiCam) {
        fiCam+=deltaFiCam;
        glutPostRedisplay();
    }

    /* Clear buffer */

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

    /* Set the camera */

	gluLookAt(rCam*sin(thetaCam)*cos(fiCam),rCam*sin(thetaCam)*sin(fiCam),rCam*cos(thetaCam),0.0,0.0,0.0,0.0,0.0,1.0);

    /* Set shade model */

    glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);

    /* Set the lighting */

    glLightfv(GL_LIGHT0,GL_POSITION,light_position);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

    /* Set material specular, shininess, ambient and diffuse */

    glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
	glMaterialfv(GL_FRONT,GL_SHININESS,mat_shininess);
    glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,qaWhite);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, qaOffWhite);

    #if SHOWFPS

    /* Compute FPS */

	frameNumber++;
	time=glutGet(GLUT_ELAPSED_TIME);
	if(time-timebase>1000.0) {
        currentFPS=frameNumber*1000.0/(time-timebase);
		timebase=time;
		frameNumber=0;
	}

    #endif

    #if SHOWENERGY

    /* Compute energy */

    energy=computeEnergy();

    #endif

    /* Create status string */

    #if SHOWFPS
        sprintf(status+strlen(status), "FPS:%4.2f ", currentFPS);
    #endif
    #if SHOWTIME
        sprintf(status+strlen(status), "[ t=%f ] ", t);
    #endif
    #if SHOWENERGY
        sprintf(status+strlen(status), "[ E=%f dE/E0=%1.4e ]", energy, (energy-initialEnergy)/initialEnergy);
    #endif

    /* Print status */
    
    #if SHOWFPS||SHOWTIME||SHOWENERGY
        glPrintString(status, 5.0, 15.0, 0.0);
    #endif

    /* Draw the bodies */

    for(i=0;i<NBODY;i++) {
        glPushMatrix();
        glTranslatef(m[i].x[0],m[i].x[1],m[i].x[2]);
        glutSolidSphere(pow(m[i].m/(4.0*M_PI*SPHEREDENSITY/3.0),1.0/3.0),50,50);
        glPopMatrix();
    }

    #if DRAWCENTEROFMASS

    /* Draw center of mass */
    
    mtotal=0.0;
    for(i=0;i<NBODY;i++) {
        mtotal=mtotal+m[i].m;
    }
    for(i=0;i<NDIM;i++) {
        xcm[i]=0.0;
        for(j=0;j<NBODY;j++) {
            xcm[i]=xcm[i]+m[j].m*m[j].x[i];
        }
        xcm[i]=xcm[i]/mtotal;
    }

    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, qaBrown);
    glPushMatrix();
    glTranslatef(xcm[0],xcm[1],xcm[2]);
    glutSolidSphere(pow(mtotal/(4.0*M_PI*SPHEREDENSITY/3.0),1.0/3.0),50,50);
    glPopMatrix();
    
    #endif

    /* Draw grid */
    
    glBegin(GL_LINES);
        for(i=0;i<=NGRID;i++) {
            if((i-(float)NGRID/2.0)==0) {
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, qaOffWhite);
            }
            else {
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, qaBlack);
            }
            glVertex3f((i-(float)NGRID/2.0)*BOXLEN/(float)NGRID,-BOXEDGE,0.0);
            glVertex3f((i-(float)NGRID/2.0)*BOXLEN/(float)NGRID,BOXEDGE,0.0);
            glVertex3f(-BOXEDGE,(i-(float)NGRID/2.0)*BOXLEN/(float)NGRID,0.0);
            glVertex3f(BOXEDGE,(i-(float)NGRID/2.0)*BOXLEN/(float)NGRID,0.0);
        };
    glEnd();

    /* Draw axes names */

    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, qaOffWhite);
    glRasterPos3f((BOXEDGE+0.2*BOXEDGE), 0.0, 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'x');
    glRasterPos3f(0.0, (BOXEDGE+0.2*BOXEDGE), 0.0);
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, 'y');
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, qaBlack);

    if(pathEnabled) {

        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, qaOffWhite);

        /* Draw paths */

        for(i=0;i<NBODY;i++) {
            glBegin(GL_LINE_STRIP);
            for(j=0;j<pathCurrent;j++) {
                glVertex3f(path[i][0][j],path[i][1][j],path[i][2][j]);
            }
            glEnd();
        }

    }

    if(momentumEnabled) {

        /* Draw momentum arrows */

        for(i=0;i<NBODY;i++) {

            glPushMatrix();
            glTranslatef(m[i].x[0],m[i].x[1],m[i].x[2]);
            glArrow(0.0,0.0,0.0,ARROWSCALE*m[i].p[0],ARROWSCALE*m[i].p[1],ARROWSCALE*m[i].p[2],ARROWDIAMETER,qaOffWhite);
            glArrow(0.0,0.0,0.0,ARROWSCALE*m[i].p[0],0.0,0.0,ARROWDIAMETER,qaOffRed);
            glArrow(0.0,0.0,0.0,0.0,ARROWSCALE*m[i].p[1],0.0,ARROWDIAMETER,qaOffGreen);
            glPopMatrix();

        }
        
    }

    /* Swap buffers */

	glFlush();
	glutSwapBuffers();

    return;
    
}

/* main */

int main(int argc, char **argv)
{

    int menu;

    /* Initialize */

    printUsage();
    printStatus(-1);
    setInitialConditions();
    initialEnergy=computeEnergy();

    /* Initialize GLUT */

    glutInit(&argc, argv);
    glutInitWindowPosition(-1,-1);
    glutInitWindowSize(WINDOWWIDTH, WINDOWHEIGHT);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    /* Create window */

    glutCreateWindow("gl2Body");

    /* Create menu */

	menu = glutCreateMenu(processMenuEvent);
	glutAddMenuEntry("Toggle Rotation",TOGGLEROTATION);
    glutAddMenuEntry("Toggle Path",TOGGLEPATH);
    glutAddMenuEntry("Toggle Momentum",TOGGLEMOMENTUM);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

    /* Set callbacks */

    glutDisplayFunc(display);
	glutReshapeFunc(resize);

    glutIgnoreKeyRepeat(1);
    glutKeyboardFunc(normalKey);
    glutKeyboardUpFunc(normalKeyRelease);
    glutSpecialFunc(specialKey);
	glutSpecialUpFunc(specialKeyRelease);

    /* Set display timer */

    glutTimerFunc(1000.0/FPS, redisplayTimer, 0);

    /* Integrate */

    printStatus(0);
    if(integrate) {
        pthread_create(&integrateThreadID, NULL, integrateThread, NULL);
    }

    /* Enable depth test */
    
    glEnable(GL_DEPTH_TEST);

    /* Main GLUT loop */

    glutMainLoop();

    return 1;

}