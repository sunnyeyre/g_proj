/*
 *  yarn.cpp
 *  FTPstats.yml
 *
 *  Created by Sunling Yang on 9/20/11.
 *  Copyright 2011 Cornell University. All rights reserved.
 *
 */
#include "Particle.h"

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

#ifdef _WIN32
static DWORD lastTime;
#else
static struct timeval lastTime;
#endif

#define PI 3.14159265

using namespace std;


Viewport	viewport;
float radius = 0.02;
float speed_factor = 50000.0;
ParticleSystem* sys;

void myReshape(int w, int h) {
	viewport.w = w;
	viewport.h = h;
	glViewport(0,0,viewport.w,viewport.h);// sets the rectangle that will be the window
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();				// loading the identity matrix for the screen
	glOrtho(0, 1, 0, 1, 1, -1);	// resize type = stretch
}

void initScene(){
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f); // Clear to black, fully transparent
	
	myReshape(viewport.w,viewport.h);
}

void myDisplay() {
	glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer (sets everything to black)
	glMatrixMode(GL_MODELVIEW);					// indicate we are specifying camera transformations
	glLoadIdentity();							// make sure transformation is "zero'd"
	glColor3f(0.0,0.0,1.0);
	glPushMatrix();
	
	// center of right yarn is at (my_x,my_y)
	float my_x = 0.5, my_y=0.5, my_r = 0.3,my_t=0.05;
	float my_d = my_r +2.0*my_t;
	glBegin(GL_POINTS);
	for(float t_inc = 0.0; t_inc < my_t; t_inc+=0.001){
		for(int ang=0; ang<=90; ang++) {
			glVertex3f((my_r+t_inc)*cos((180-ang)*PI/180.0)+my_x,(my_r+t_inc)*sin((180-ang)*PI/180.0)+my_y,0.0f);
			glVertex3f((my_r+t_inc)*cos(ang*PI/180.0)+my_x-my_d,(my_r+t_inc)*sin(ang*PI/180.0)+my_y,0.0f);
		}}
	glEnd();	
	glTranslatef((my_r+my_t/2.0)*cos(49.0*PI/180.0)+my_x-my_d,(my_r+my_t/2.0)*sin(49.0*PI/180.0)+my_y,0.0f);
	glutSolidSphere(0.02,20,20);
	glPopMatrix();
	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

void processNormalKeys(unsigned char key, int x, int y) {
	
	if (key == 32) //the space key
	{ 
		float temp,deltat,cos_delta,sin_delta,deltax,deltay,abs_delta,u1,u2;
		float temp_x_n[sys->n], temp_y_n[sys->n], n[2], nt[2],v1[2], v2[2];
		float temp_t = sys->t;
		sys->t = glutGet(GLUT_ELAPSED_TIME);
		temp_t = (sys->t - temp_t)/speed_factor;
		for(int i=0; i<sys->n; i++){
			deltax = sys->p[i]->x[0]+sys->p[i]->v[0]*temp_t;
			deltay = sys->p[i]->x[1]+sys->p[i]->v[1]*temp_t;
			temp_x_n[i] = deltax;
			temp_y_n[i] = deltay;
			/* collision checking                  */
			for(int j=0; j<i; j++){
				deltax = sys->p[i]->x[0]+sys->p[i]->v[0]*temp_t;
				deltay = sys->p[i]->x[1]+sys->p[i]->v[1]*temp_t;
				deltax = deltax - temp_x_n[j];
				deltay = deltay - temp_y_n[j];
				abs_delta = deltax*deltax + deltay*deltay;
				// in case of a collision, return particles to the moment they collide, swap their velocities and evolve to equal time w/
				// new velocities
				if(abs_delta <= 4*radius*radius){
					abs_delta = sqrt(abs_delta);
					n[0] = deltax / abs_delta;
					n[1] = deltay / abs_delta;
					nt[0] = -n[1];
					nt[1] = n[0];
					v1[0] = sys->p[i]->v[0]*n[0] + sys->p[i]->v[1]*n[1]; //v1n
					v1[1] = sys->p[i]->v[0]*nt[0] + sys->p[i]->v[1]*nt[1]; //v1nt
					v2[0] = sys->p[j]->v[0]*n[0] + sys->p[j]->v[1]*n[1];  //v2n
					v2[1] = sys->p[j]->v[0]*nt[0] + sys->p[j]->v[1]*nt[1]; // v2nt
					u1 = sys->p[i]->x[0]*n[0] + sys->p[i]->x[1]*n[1];
					u2 = sys->p[j]->x[0]*n[0] + sys->p[j]->x[1]*n[1];
					if(abs(v1[0] - v2[0]) < 0.00001) {
						v1[0] *= -1;	
					}
					else {
						/************* find time of collision and save in deltat*******/
						deltat = (2*radius-u1+u2)/(v1[0]-v2[0]); //time of collision
						deltat = deltat<0.0? -deltat:deltat;
						temp_x_n[j] = sys->p[j]->x[0]+sys->p[j]->v[0]*deltat;
						temp_y_n[j] = sys->p[j]->x[1]+sys->p[j]->v[1]*deltat;
						temp_x_n[i] = sys->p[i]->x[0]+sys->p[i]->v[0]*deltat;
						temp_y_n[i] = sys->p[i]->x[1]+sys->p[i]->v[1]*deltat;
						
						/*********** swap velocities here ************************/
						temp = v1[0]; //normal to axis is swapped, tangent stays the same
						v1[0] = v2[0];
						v2[0] = temp;
						cos_delta = deltax/abs_delta;
						sin_delta = deltay/abs_delta;
					}
					sys->p[i]->v[0] = v1[0]*cos_delta - v1[1]*sin_delta; //rotated back to xy axis
					sys->p[i]->v[1] = +v1[0]*sin_delta + v1[1]*cos_delta;
					sys->p[j]->v[0] = v2[0]*cos_delta - v2[1]*sin_delta; //rotated back to xy axis
					sys->p[j]->v[1] = +v2[0]*sin_delta + v2[1]*cos_delta;
					
					/************ evolve to equal time w/ new velocities here*****/
					temp_x_n[j] += sys->p[j]->v[0]*(temp_t - deltat);
					temp_x_n[j] += sys->p[j]->v[1]*(temp_t - deltat);
					temp_x_n[i] += sys->p[i]->v[0]*(temp_t - deltat);
					temp_x_n[i] += sys->p[i]->v[1]*(temp_t - deltat);
				}}
		}
		for(int i=0; i<sys->n; i++){
			/* boundary checking */
			if(temp_x_n[i] - radius < 0.0 || temp_x_n[i] + radius > 1.0) {
				sys->p[i]->v[0] = -sys->p[i]->v[0];
			}
			else
				sys->p[i]->x[0] = temp_x_n[i];
			
			if(temp_y_n[i] - radius < 0.0 || temp_y_n[i] + radius > 1.0) {
				sys->p[i]->v[1] = -sys->p[i]->v[1];
			}
			else
				sys->p[i]->x[1] = temp_y_n[i];
		}
		glutPostRedisplay();
		
	}
	if (key == 27) {
		// when receiving the escape key the function escapes
		exit(0);
	}
}

void myFrameMove() {
	
	glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}
int main(int argc, char *argv[]) {
	srand((unsigned)time(0));
	
  	glutInit(&argc, argv);
  	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	
  	viewport.w = 500;
  	viewport.h = 500;
	
	sys = new ParticleSystem;
	sys->n = 350;
	sys->t = glutGet(GLUT_ELAPSED_TIME);
	sys->p = new Particle_t[sys->n];
	int coin = 0;
	for(int i=0; i< sys->n; i++){
		sys->p[i] = new Particle;
		sys->p[i]->x = new float[2];
		sys->p[i]->v = new float[2];
		sys->p[i]->x[0] = (float) (i%10)/10.0 + 0.05;
		sys->p[i]->x[1] = (float) (i/10)/10.0 + 0.05;
		sys->p[i]->v[0] = 0.7;
		sys->p[i]->v[1] = 0.7;
		
		coin = rand()%2;
		if(coin)
			sys->p[i]->v[0] *= -1;
		coin = rand()%2;
		if(coin)
			sys->p[i]->v[1] *= -1;
	}	
  	glutInitWindowSize(viewport.w, viewport.h);
  	glutInitWindowPosition(0, 0);
  	glutCreateWindow("Toy Model");
	
  	initScene();							// quick function to set up scene
	glutSwapBuffers();
	glutDisplayFunc(myDisplay);
	
	glutKeyboardFunc(processNormalKeys);
  	glutReshapeFunc(myReshape);				// function to run when the window gets resized
  	glutIdleFunc(myFrameMove);				// function to run when not handling any other task
  	glutMainLoop();							// infinite loop that will keep drawing and resizing and whatever else
	
	for(int i=0; i<sys->n; i++){
		delete(sys->p[i]->x);
		delete(sys->p[i]->v);
		delete(sys->p[i]);
	}
	delete(sys->p);
	delete(sys);
	
  	return 0;
}


