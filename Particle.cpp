/*
 *  Particle.cpp
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

#define PI 3.14159265

using namespace std;
extern float speed_factor;
extern float radius;

void ParticleGetState(ParticleSystem_t sys, float *dst) {
	int i;
	for(i=0; i< sys->n; i++){
		*(dst++) = sys->p[i]->x[0];
		*(dst++) = sys->p[i]->x[1];
		*(dst++) = sys->p[i]->x[2];
		*(dst++) = sys->p[i]->v[0];
		*(dst++) = sys->p[i]->v[1];
		*(dst++) = sys->p[i]->v[2];
	}
}

int ParticleDims(ParticleSystem_t sys) {
	return (6 * sys->n);
}

void ParticleSetState(ParticleSystem_t sys, float *src) {
	int i;
	for (i=0; i< sys->n; i++){
		sys->p[i]->x[0] = *(src++);
		sys->p[i]->x[1] = *(src++);
		sys->p[i]->x[2] = *(src++);
		sys->p[i]->v[0] = *(src++);
		sys->p[i]->v[1] = *(src++);
		sys->p[i]->v[2] = *(src++);
	}
}

void ParticleDerivative(ParticleSystem_t sys, float *dst) {
	int i;
	Clear_Forces(sys);
	Compute_Forces(sys);
	for (i=0; i< sys->n; i++) {
		*(dst++) = sys->p[i]->v[0];
		*(dst++) = sys->p[i]->v[1];
		*(dst++) = sys->p[i]->v[2];
		*(dst++) = sys->p[i]->f[0]/sys->p[i]->m;
		*(dst++) = sys->p[i]->f[1]/sys->p[i]->m;
		*(dst++) = sys->p[i]->f[2]/sys->p[i]->m;
	}
}

void EulerStep(ParticleSystem_t sys, float DeltaT) { //DeltaT = glutGet(GLUT_ELAPSED_TIME) - sys->t
	float * temp1 = new float(ParticleDims(sys));
	float * temp2 = new float(ParticleDims(sys));
	ParticleDerivative(sys, temp1);
	ScaleVector(temp1, DeltaT);
	ParticleGetState(sys, temp2);
	AddVectors(temp1, temp2, temp2);
	ParticleSetState(sys, temp2);
	sys->t += DeltaT;
}

void AddVectors(float *num1, float * num2, float * result) {
	int i;
	int lengthOfArray = sizeof(num1) / sizeof(float);	
	for (i=0; i< lengthOfArray; i++) {
		 *(result++) = *(num1++) + *(num2++);
	}
}

void ScaleVector(float * nums, float DeltaT) {
	int i;
	int lengthOfArray = sizeof(nums) / sizeof(float);
	for (i=0; i< lengthOfArray; i++) {
		nums[i] *= DeltaT;
	}
}

void Clear_Forces(ParticleSystem_t sys) {
	int i;
	for (i=0; i<sys->n; i++) {
		sys->p[i]->f[0] = 0;
		sys->p[i]->f[1] = 0;
		sys->p[i]->f[2] = 0;
	}
}

void Compute_Forces(ParticleSystem_t sys) {
}

void move_forward(ParticleSystem_t sys) {
	float temp_t = glutGet(GLUT_ELAPSED_TIME);
	temp_t = (sys->t - temp_t)/speed_factor;
	EulerStep(sys, temp_t);
}

void move_backward(ParticleSystem_t sys) {
	float temp_t = glutGet(GLUT_ELAPSED_TIME);
	temp_t = (sys->t - temp_t)/speed_factor;
	EulerStep(sys, -temp_t);
}