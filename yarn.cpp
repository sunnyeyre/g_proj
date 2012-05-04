#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>
#include "algebra3.h"
#include <Eigen/Dense>
#include "Energy_Func.h"

#define PI 3.14159265
#define g 9.80665

using namespace std;
//using namespace Eigen;

//******************************************************************
//  Some Classes
//******************************************************************

nodeSystem* sys;

class Viewport {
public:
	int w,h;
};

float deltaT = 1;
float rho = 1.0;

Viewport viewport;
int counter = 1;
int display_counter = 1;
float radius = 0.02;
float win_x, win_y, win_z;
int rotx = 0;
int roty = 0;
float amt_x=0.0; float amt_y=0.0;
int frame=0,my_time,timebase=0;

void myReshape(int w, int h) {
	glViewport(0,0, (GLsizei) w, (GLsizei) h);// sets the rectangle that will be the window
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();				// loading the identity matrix for the screen
	gluPerspective(65.0, (GLfloat) w/(GLfloat) h, 1.0, 20.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
//	glTranslatef(0.0, 0.0, -5.0);
}

void buildFrame() {
	counter++;
	glNewList(counter, GL_COMPILE);
		
	glutSolidSphere(0.1,20,20);

	glRotatef(sys->n[1].x[0],sys->n[1].x[1],sys->n[1].x[2],1);
//	glRotatef(30,0,0,1);
	glPushMatrix();
	glBegin(GL_LINES);
//	glVertex3f(0,0,0);
//	glVertex3f(sys->n[1].x[0]-sys->n[2].x[0],sys->n[1].x[1]-sys->n[2].x[1],sys->n[1].x[2]-sys->n[2].x[2]);
	glVertex3f(0,0,0);
	glVertex3f(0,-2,0);
	glEnd();
	
	glTranslatef(0,-2,0);
//	glTranslatef(sys->n[1].x[0]-sys->n[2].x[0],sys->n[1].x[1]-sys->n[2].x[1],sys->n[1].x[2]-sys->n[2].x[2]);
	glutSolidSphere(0.1,20,20);
	glPopMatrix();
	
	glEndList();
} 

void animate() {

	takeAStep();
	buildFrame();
	//save to jpeg
  //	glutTimerFunc(33.33333, animate, 0);    //to save to buffer ~30frames/sec
}

void initScene(){
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f); // Clear to black, fully transparent
	glClearDepth(0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glTranslatef(0.0, 0.0, -5.0);
	glPushMatrix();
	
	glNewList(counter, GL_COMPILE);
	glColor3f(0.0,0.0,1.0);

	
	glutSolidSphere(0.1,20,20);

//	glRotatef(30,0,0,1);
	glRotatef(sys->n[1].x[0],sys->n[1].x[1],sys->n[1].x[2],1);
	glPushMatrix();
	glBegin(GL_LINES);
	glVertex3f(0,0,0);
	glVertex3f(0,-2,0);
	glEnd();
	
	glTranslatef(0,-2,0);
	glutSolidSphere(0.1,20,20);
	glPopMatrix();
	
	glEndList();
	glShadeModel (GL_FLAT);
	
//	glutTimerFunc(0,animate,0);
}

void initLights() {
	glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	GLfloat global_ambient[] = {.1f, .1f, .1f};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	
	GLfloat ambient[] = {0.1f, 0.1f, 0.1f};
	GLfloat diffuse[] = {1.0f, 1.0f, 1.0f};
	GLfloat specular[] = {1.0f, 1.0f, 1.0f};
	GLfloat pos[] = {0, 0, 1, 0};
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT0, GL_POSITION, pos);
	glEnable(GL_LIGHT0);
	
	GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_diffuse[] = {1, 0, 1, 1.0};
	GLfloat mat_ambient[] = {.1, .1, .1, 1.0};
	GLfloat mat_shininess[] = {50.0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
}

void mouseTracker(int button, int state, int x, int y){
	if(button == 0 && state == GLUT_DOWN){
		win_x = (x/viewport.w - 0.5)* 10.0;
		win_y = (viewport.h/2.0 - y) / viewport.h * 10.0;
		win_z = 0.0;
		printf("%d %d, %f %f\n", x, y, win_x, win_y);
		glutPostRedisplay();
	}
}

//void EulerStep() {
//	float old_t = sys->t;
//	sys->t = glutGet(GLUT_ELAPSED_TIME);
//	float deltaT = (sys->t - old_t)/100.0;
//	float deltaThetaDot = g*sin(sys->n[1].x[0])*deltaT/2.0;
//	sys->n[1].xdot[0] += deltaThetaDot;
//	sys->n[1].x[0] += sys->n[1].xdot[0]*deltaT;
//	sys->t = glutGet(GLUT_ELAPSED_TIME);
////	sys->n[1].x[0] = 30*cos(sqrt(g/2.0)*sys->t - t0);
//}
//

void myDisplay() {
	
	animate();
	
	glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer (sets everything to black)
	glViewport(0,0, (GLsizei) viewport.w, (GLsizei) viewport.h);// sets the rectangle that will be the window
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();				// loading the identity matrix for the screen
	gluPerspective(65.0, (GLfloat) viewport.w/(GLfloat) viewport.h, 1.0, 20.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glTranslatef(0.0, 0.0, -5.0);
	glTranslatef(amt_x,amt_y,0);
	glRotatef(rotx,1,0,0);
	glRotatef(roty,0,1,0);
	glPushMatrix();
	
	glColor3f(0,1,0);
	display_counter++;
	glCallList(display_counter);
	
//	glutSolidSphere(0.1,20,20);
//	
//	glRotatef(sys->n[1].x[0],sys->n[1].x[1],sys->n[1].x[2],1);
//	//	glRotatef(30,0,0,1);
//	glPushMatrix();
//	glBegin(GL_LINES);
//	//	glVertex3f(0,0,0);
//	//	glVertex3f(sys->n[1].x[0]-sys->n[2].x[0],sys->n[1].x[1]-sys->n[2].x[1],sys->n[1].x[2]-sys->n[2].x[2]);
//	glVertex3f(0,0,0);
//	glVertex3f(0,-2,0);
//	glEnd();
//	
//	glTranslatef(0,-2,0);
//	//	glTranslatef(sys->n[1].x[0]-sys->n[2].x[0],sys->n[1].x[1]-sys->n[2].x[1],sys->n[1].x[2]-sys->n[2].x[2]);
//	glutSolidSphere(0.1,20,20);
//	glPopMatrix();
	
	
	
	glShadeModel(GL_FLAT);
	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

void processSpecialKeys(int key, int x, int y) {
	int mod = glutGetModifiers();
	switch (key) { //left and right key for moving the frame rate
		case GLUT_KEY_LEFT:
			if(mod == GLUT_ACTIVE_SHIFT) {
				roty -= 10;
				roty %= 360;
			}
			else {
				if(amt_x+0.1 <=3)
					amt_x -= 0.1;
			}
			break;
		case GLUT_KEY_F1:
			if(display_counter -20 >= 1)
				display_counter -= 20;
			glutPostRedisplay();
			break;
		case GLUT_KEY_F2:
			if(display_counter + 20 <= counter)
				display_counter += 20;
			glutPostRedisplay();
			break;
		case GLUT_KEY_RIGHT:
			if(mod == GLUT_ACTIVE_SHIFT) {
				roty += 10;
				roty %= 360;
			}
			else {
				if(amt_x - 0.1 >= -3)
					amt_x += 0.1;
				glutPostRedisplay();
			}
			break;
		case GLUT_KEY_UP:
			if(mod == GLUT_ACTIVE_SHIFT){
				rotx += 10;
				rotx %= 360;
			}
			else {
				if(amt_y - 0.1 >= -3)
					amt_y += 0.1;
			}
			break;
		case GLUT_KEY_DOWN:
			if(mod == GLUT_ACTIVE_SHIFT) {
				rotx -= 10;
				rotx %= 360;
			}
			else {
				if(amt_y + 0.1 <= 3)
					amt_y -= 0.1;
			}
			break;
	}
}
	
void processNormalKeys(unsigned char key, int x, int y) {

	if (key == 32) //the space key
		glutPostRedisplay();

	if (key == 27) {
// when receiving the escape key the function escapes
		exit(0);
	}
}
 
void myFrameMove() {

	glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}

// initialization of nodes, can change into a read from an initialization text file
void node_system_initialize() {
	sys = new nodeSystem;
	sys->num = 2;
	sys->n = new struct node[2];

	sys->n[0].rho = sys->n[1].rho = rho;
	sys->n[0].x = new float[3];
	sys->n[0].x[0] = 0.0;
	sys->n[0].x[1] = 0.0;
	sys->n[0].x[2] = 0.0;
	sys->n[1].x = new float[3];
	sys->n[1].x[0] = 45.0;
	sys->n[1].x[1] = 0.0;
	sys->n[1].x[2] = 0.0;
	sys->n[0].xdot = new float[3];
	sys->n[0].xdot[0] = 0.0;
	sys->n[0].xdot[1] = 0.0;
	sys->n[0].xdot[2] = 0.0;
	sys->n[1].xdot = new float[3];
	sys->n[1].xdot[0] = 0.1;
	sys->n[1].xdot[1] = 0.1;
	sys->n[1].xdot[2] = 0.1;
	sys->n[0].s = 0.0;
	sys->n[1].s = 0.0;
	sys->n[0].sdot = 0.0;
	sys->n[1].sdot = 0.0;
}

int main(int argc, char *argv[]) {

	
	srand((unsigned)time(0));
	node_system_initialize();
	
	vec3 deltax = vec3(sys->n[1].x[0]-sys->n[0].x[0],sys->n[1].x[1]-sys->n[0].x[1],sys->n[1].x[2]-sys->n[0].x[2]);
	input_M_matrix(sys->n[1].s - sys->n[0].s, deltax, rho, radius*radius);
	input_G_matrix(1, deltax, sys->n[1].s - sys->n[0].s, radius);
	
	
  	glutInit(&argc, argv);
  	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  	viewport.w = 500;
  	viewport.h = 500;
		
  	glutInitWindowSize(viewport.w, viewport.h);
  	glutInitWindowPosition(0, 0);
  	glutCreateWindow("Pendulum");

  	initScene();							// quick function to set up scene
	glutSwapBuffers();
	glutDisplayFunc(myDisplay);

	glutKeyboardFunc(processNormalKeys);
	glutSpecialFunc(processSpecialKeys);
  	glutReshapeFunc(myReshape);				// function to run when the window gets resized
  	initLights();
	glEnable(GL_DEPTH_TEST);
	glutMouseFunc(mouseTracker);
	glutIdleFunc(myFrameMove);				// function to run when not handling any other task
	glutMainLoop();							// infinite loop that will keep drawing and resizing and whatever else

  	return 0;
}

