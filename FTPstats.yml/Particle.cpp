/*
 *  Particle.cpp
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

void move_forward(Particle_System * sys) {
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
		}
	}