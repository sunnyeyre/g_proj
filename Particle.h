/*
 *  Particle.h
 *
 */
#ifndef PARTICLE_H
#define PARTICLE_H
#endif
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

class Viewport {
public:
	int w, h; // width and height
};

/****** declaring particle and particle system structs ******/
typedef struct Particle {
	float m;
	float *x;
	float *v;
	float *f;
};

typedef struct Particle * Particle_t;

typedef struct ParticleSystem{
	Particle_t * p;
	int n;
	float t;
};
typedef struct ParticleSystem * ParticleSystem_t;

/****** performs calculations to 'evolve' the particle system ******/

void EulerStep(ParticleSystem_t sys, float DeltaT);
typedef struct ParticleSystem * ParticleSystem_t;
int ParticleDims(ParticleSystem_t sys);
void ParticleGetState(ParticleSystem_t sys, float *dst);
void ParticleSetState(ParticleSystem_t sys, float *src);
void ParticleDerivative(ParticleSystem_t sys, float *dst);
void move_forward(ParticleSystem_t sys);

/**** helper functions to make my life easier *****/
void AddVectors(float * num1, float * num2, float * result);
void ScaleVector(float * nums, float DeltaT);
void Clear_Forces(ParticleSystem_t sys);
void Compute_Forces(ParticleSystem_t sys);
