/*
 *  Particle.h
 *  FTPstats.yml
 *
 *  Created by Sunling Yang on 9/20/11.
 *  Copyright 2011 Cornell University. All rights reserved.
 *
 */
#ifndef "Particle.h"
#define "Particle.h"
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

Class Particle {
public:
	typedef struct {
		float m;
		float *x;
		float *v;
		float *f;
	} *Particle;

	typedef struct {
		Particle_t * p;
		int n;
		float t;
	} *ParticleSystem;

	void move_forward();
};