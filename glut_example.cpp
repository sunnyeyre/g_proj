// Simple OpenGL example for CS184 sp08 by Trevor Standley, modified from sample code for CS184 on Sp06
// modified by Sunling Yang on Fa10
//Note: I got the rotating and keystroke codes from
// www.lighthouse3d.com/opengl/glut/index.php
// Note: the circle drawing I got from 
// cboard.cprogramming.com "drawing a circle using glut"

//********************************************************************
// Note: <space> switches between scenes and <esc> escapes the program
//
//********************************************************************

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#ifdef _WIN32
#	include <windows.h>
#else
#	include <sys/time.h>
#endif

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <time.h>
#include <math.h>

#ifdef _WIN32
static DWORD lastTime;
#else
static struct timeval lastTime;
#endif

#define PI 3.14159265

using namespace std;

//****************************************************
// Some Classes
//****************************************************
class Viewport {
public:
	int w, h; // width and height
};


//****************************************************
// Global Variables
//****************************************************
Viewport	viewport;
float angle1=0.0;
int SceneNum = 0;
float x,y;
float ybase = 0.37;
float radiusx = 0.04;
float radiusy = 0.0625;
float theoCounter = 0.0;
float angle_theo = 0.0;
static float denom = sqrt(17);                // want the last triangle to be biggest
static float deltaColor = 0.1;
float hypo;
float x_theo = 1.0;
float y_theo = 0.0;

//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
	viewport.w = w;
	viewport.h = h;

	glViewport(0,0,viewport.w,viewport.h);// sets the rectangle that will be the window
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();				// loading the identity matrix for the screen

	//----------- setting the projection -------------------------
	// glOrtho sets left, right, bottom, top, zNear, zFar of the chord system

	// glOrtho(-1, 1 + (w-400)/200.0 , -1 -(h-400)/200.0, 1, 1, -1); // resize type = add
	// glOrtho(-w/400.0, w/400.0, -h/400.0, h/400.0, 1, -1); // resize type = center

	glOrtho(-1, 1, -1, 1, 1, -1);	// resize type = stretch

	//------------------------------------------------------------
}


//****************************************************
// sets the window up
//****************************************************
void initScene(){
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f); // Clear to black, fully transparent
	
	myReshape(viewport.w,viewport.h);
}


//***************************************************
// function that does the actual drawing
//***************************************************
void myDisplay() {
	
	glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer (sets everything to black)
	
	glMatrixMode(GL_MODELVIEW);					// indicate we are specifying camera transformations
	glLoadIdentity();							// make sure transformation is "zero'd"
	
	//----------------------- code to draw objects --------------------------
							// Rectangle Code
	//glColor3f(red component, green component, blue component);
	glColor3f(1.0f,0.0f,0.0f);					// setting the color to pure red 90% for the rect
	
	glBegin(GL_POLYGON);						// draw rectangle 
		//glVertex3f(x val, y val, z val (won't change the point because of the projection type));
		glVertex3f(-0.8f, 0.0f, 0.0f);			// bottom left corner of rectangle
		glVertex3f(-0.8f, 0.5f, 0.0f);			// top left corner of rectangle
		glVertex3f( 0.0f, 0.5f, 0.0f);			// top right corner of rectangle
		glVertex3f( 0.0f, 0.0f, 0.0f);			// bottom right corner of rectangle
	glEnd();
							// Triangle Code
	glColor3f(1.0f,0.5f,0.0f);					// setting the color to orange for the triangle
	
	float basey = -sqrt(0.48f);					// height of triangle = sqrt(.8^2-.4^2)
	glBegin(GL_POLYGON);
		glVertex3f(0.5f,  0.0f, 0.0f);			// top tip of triangle
		glVertex3f(0.1f, basey, 0.0f);			// lower left corner of triangle
		glVertex3f(0.9f, basey, 0.0f);			// lower right corner of triangle
	glEnd();
	//-----------------------------------------------------------------------
	
	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

//**********************************************************
// function that takes a space and moves to the next function
//**********************************************************

void processNormalKeys(unsigned char key, int x, int y) {
	if (key == 32) //the space key
		{ 
			SceneNum = (SceneNum + 1)%2; //oscillates between scenes
		}
	if (key == 27) {
// when receiving the escape key the function escapes
		exit(0);
	}
}
//****************************************************
// function that starts animation
//****************************************************

void renderScene(void) {

glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

if(SceneNum == 0) {// attempt at Spiral of Theodorus

 glPushMatrix();
 for (int i=0; i<18; i++) {
	 theoCounter = i;
  hypo = sqrt(1.0 + theoCounter);
  angle_theo += atan(1/hypo);
	glColor3f(0.1+deltaColor,0.1f+deltaColor,0.8f);					
	glBegin(GL_POLYGON);
		glVertex3f(0.0f,0.0f,0.0f);			
		glVertex3f(hypo*cos(angle_theo)/denom, hypo*sin(angle_theo)/denom, 0.0f);		
		glVertex3f(x_theo/denom, y_theo/denom, 0.0f);			
	glEnd();
	
 glPopMatrix();
	glFlush();
	glutSwapBuffers();

  deltaColor += 0.1;
 // theoCounter++;
  x_theo=hypo*cos(angle_theo);
  y_theo=hypo*sin(angle_theo);	
}
  deltaColor = 0;
}
if (SceneNum == 1) {// attempt at vitruvian man
glPushMatrix();

glRotatef(angle1, 0.0,1.0,0.0);
	glColor3f(0.0f,0.0f,0.9f);		
glBegin(GL_TRIANGLES);
	glVertex3f(1.0,0.0,0.0);
	glVertex3f(0.95,0.025,0.0);
	glVertex3f(0.9,-0.025,0.0);
glEnd();

glRotatef(angle1, 0.0,0.0,1.0); //rotates everything by angle degrees
	glColor3f(0.0f,0.0f,0.9f);		
glBegin(GL_TRIANGLES);
	glVertex3f(0.0,1.0,0.0);
	glVertex3f(0.025,0.95,0.0);
	glVertex3f(-0.025,0.9,0.0);
glEnd();

glRotatef(angle1, 0.0,1.0,0.0); // trunk
	glColor3f(0.7f,0.5f,0.0f);		
glBegin(GL_POLYGON);
	glVertex3f(-0.11,0.0,0.0);
	glVertex3f(-0.11,0.32,0.0);
	glVertex3f(0.0,0.375,0.0);
	glVertex3f(0.11,0.32,0.0);
	glVertex3f(0.11,0.0,0.0);
glEnd();

glRotatef(angle1, 0.0,1.0,0.0); //left thigh
	glColor3f(0.7f,0.5f,0.0f);		
glBegin(GL_POLYGON);
	glVertex3f(-0.08,-0.25,0.0);
	glVertex3f(-0.11,0.0,0.0);
	glVertex3f(-0.02,0.0,0.0);
	glVertex3f(-0.02,-0.25,0.0);
glEnd();

glRotatef(angle1, 0.0,1.0,0.0); //right thigh
	glColor3f(0.7f,0.5f,0.0f);		
glBegin(GL_POLYGON);
	glVertex3f(0.02,-0.25,0.0);
	glVertex3f(0.02,0.0,0.0);
	glVertex3f(0.11,0.0,0.0);
	glVertex3f(0.08,-0.25,0.0);
glEnd();

glRotatef(angle1, 0.0,1.0,0.0); //left leg 
	glColor3f(0.7f,0.5f,0.0f);		
glBegin(GL_POLYGON);
	glVertex3f(-0.08,-0.5,0.0);
	glVertex3f(-0.08,-0.25,0.0);
	glVertex3f(-0.02,-0.25,0.0);
	glVertex3f(-0.02,-0.5,0.0);
glEnd();

glRotatef(angle1, 0.0,1.0,0.0); //right leg
	glColor3f(0.7f,0.5f,0.0f);		
glBegin(GL_POLYGON);
	glVertex3f(0.02,-0.5,0.0);
	glVertex3f(0.02,-0.25,0.0);
	glVertex3f(0.08,-0.25,0.0);
	glVertex3f(0.08,-0.5,0.0);
glEnd();

glRotatef(angle1, 0.0,1.0,0.0); // arm span
	glColor3f(0.7f,0.5f,0.0f);		
glBegin(GL_POLYGON);
	glVertex3f(-0.5,0.28,0.0);
	glVertex3f(-0.5,0.32,0.0);
	glVertex3f(0.5,0.32,0.0);
	glVertex3f(0.5,0.28,0.0);
glEnd();

glRotatef(angle1, 0.0,1.0,0.0); // head 
	glColor3f(0.7f,0.5f,0.0f);		
 glBegin(GL_TRIANGLE_FAN);
 	 glVertex2f(0.0, ybase);
 for (int j=0; j<360; j++) {
	 glVertex2f( (float) radiusx * sin(j * PI/ 180.0), ybase + (float) radiusy * cos(j*PI/180.0));
}
glEnd();

	glColor3f(0.7f,0.0f,0.3f);		
glBegin(GL_LINES);         // a ring
	x = (float) 0.95 * cos(359 * PI/180.0);
	y = (float) 0.95 * sin(359 * PI/180.0);
	for (int j = 0; j<360; j=j++) {
		glVertex2f(x,y);
		x = (float) 0.95 * cos(j * PI / 180.0);
		y = (float) 0.95 * sin(j * PI / 180.0);
		glVertex2f(x,y);
		}
glEnd();

glPopMatrix();
glutSwapBuffers();

angle1++;
}

}
//****************************************************
// called by glut when there are no messages to handle
//****************************************************
void myFrameMove() {
	//nothing here for now
#ifdef _WIN32
	Sleep(1000);						//give ~10ms back to OS (so as not to waste the CPU)
#endif
	glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}


//****************************************************
// the usual stuff, nothing exciting here
//****************************************************
int main(int argc, char *argv[]) {
  	//This initializes glut
  	glutInit(&argc, argv);
  
  	//This tells glut to use a double-buffered window with red, green, and blue channels 
  	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

  	// Initalize theviewport size
  	viewport.w = 400;
  	viewport.h = 400;

  	//The size and position of the window
  	glutInitWindowSize(viewport.w, viewport.h);
  	glutInitWindowPosition(0, 0);
  	glutCreateWindow("CS184!");

//  	initScene();							// quick function to set up scene

  
//  	glutDisplayFunc(myDisplay);				// function to run when its time to draw something
   glutSwapBuffers();
	glutDisplayFunc(renderScene);


	glutKeyboardFunc(processNormalKeys);

  	glutReshapeFunc(myReshape);				// function to run when the window gets resized
  	glutIdleFunc(myFrameMove);				// function to run when not handling any other task
	glEnable(GL_DEPTH_TEST);
  	glutMainLoop();							// infinite loop that will keep drawing and resizing and whatever else
  
  	return 0;
}








