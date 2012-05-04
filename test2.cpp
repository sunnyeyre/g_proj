/*
 *  Energy_Func.cpp
 *  
 *
 */

#include "Energy_Func.h"
#include "linalg.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_min.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>

//#ifdef OSX
//#include <GLUT/glut.h>
//#include <OpenGL/glu.h>
//#else
//#include <GL/glut.h>
//#include <GL/glu.h>
//#endif

#include <time.h>
#include "algebra3.h"
#include <Eigen/Dense>

#define PI 3.14159265
#define g 9.80665

using namespace std;
using namespace Eigen;

mat_t M;
mat_t G;
int cons_row;
float radius = 1.0;
float rho = 1.0;
float deltaT = 1.0;
nodeSystem* sys;

void input_M_matrix(float deltas, vec3& deltax, float rho, float X) {
//inputs and updates M matrix given X, deltas, deltax, and rho
	
	float mxx = deltas * X;
	vec3 mxs = vec3(deltax[0], deltax[1], deltax[2]) * X;
	float mss = deltax * deltax * X / deltas;
	
	M = new mat;
	M->rows = new matrix_row[4];
	M->rows[0].s0 = M->rows[1].s1 = 2*mxx;
	M->rows[0].s1 = M->rows[1].s0 = mxx;
	M->rows[0].x0 = -2.0*mxs;
	M->rows[0].x1 = -mxs;
	
	M->rows[1].s0 = mxx;
	M->rows[1].s1 = 2*mxx;
	M->rows[1].x0 = -mxs;
	M->rows[1].x1 = -2*mxs;
	
	M->rows[2].s0 = 2*mss;
	M->rows[2].s1 = mss;
	M->rows[2].x0 = -2*mxs;
	M->rows[2].x1 = -mxs;
	
	M->rows[3].s0 = mss;
	M->rows[3].s1 = 2*mss;
	M->rows[3].x0 = -mxs;
	M->rows[3].x1 = -2*mxs;
	
	for(int i=0; i<4; i++) {
		M->rows[i].s0 *= (6.0);
		M->rows[i].s1 *= (6.0);
		M->rows[i].x0 *= (6.0);
		M->rows[i].x1 *= (6.0);
	}
}
	
void input_G_matrix(int cons_row, vec3& deltax, float deltas, float radius) {
//inputs and updates the lagrangian matrix given parameters
	
	G = new mat;
	G->rows = new matrix_row[cons_row]; //not using num_params because it is 1 in this case
	
	G->rows[0].s0 = radius*radius/2.0;
	G->rows[0].s1 = radius*radius/2.0;
	G->rows[0].x0 = -radius*radius*deltax/(2.0*deltas);
	G->rows[0].x1 = -radius*radius*deltax/(2.0*deltas);
	
}

matrix_row & matrix_row::operator *= (vec3& b) {
	matrix_row* temp = new matrix_row;
	temp->s0 = this->x0 * b;
	temp->s1 = this->x1 * b;
	temp->x0 = this->s0 * b;
	temp->x1 = this->s1 * b;
	
	this->x0 = temp->x0;
	this->x1 = temp->x1;
	this->s0 = temp->s0;
	this->s1 = temp->s1;
	free(temp);
	return *this;
}

matrix_row & matrix_row::operator = (matrix_row & b) {
	
	this->s0 = b.s0;
	this->s1 = b.s1;
	this->x0 = b.x0;
	this->x1 = b.x1;
	return *this;
}

matrix_row & matrix_row::operator += (matrix_row& b) {
	
	this->s0 += b.s0;
	this->s1 += b.s1;
	this->x0 += b.x0;
	this->x1 += b.x1;
	return *this;
}

matrix_row & matrix_row::operator *= (float b) {
	this->s0 *= b;
	this->s1 *= b;
	this->x0 *= b;
	this->x1 *= b;
	return *this;
}

vec3 get_interpolated_x(float s0, float s1, float s, vec3& x0, vec3& x1) {
//computes x tilda using material coordinate s
	double alpha = (s-s0)/(s1-s0);
	return vec3((1-alpha)*x0 + alpha*x1);
}

vec3 get_interpolated_xdot(float s, float s0, float s1, float sdot0, float sdot1, vec3& xdot0, vec3& xdot1, vec3& deltax, vec3& derivx) {
//computes actual worlds coordinates using the constraint function, for pendulum

	double alpha = (s-s0)/(s1-s0);
	vec3 temp = vec3(deltax);
	temp *= (1/(s1-s0));
	return vec3(derivx * ((1-alpha) * xdot0 + alpha * xdot1 - (temp)*( (1-alpha)*sdot0 + alpha*sdot1)));
}

matrix_row* get_force(float deltas, vec3& x0, vec3& x1, float rho) {
//computes dV/dq and multiplies it by inverse of M and output to acceleration
	
	matrix_row* force = new matrix_row;
	
	double constant = rho*g/2.0;
	force->s0 = -deltas*constant;
	force->s1 = -deltas*constant;
	force->x0 = (x0 + x1)*constant;
	force->x1 = -(x0 + x1)*constant;
	
	return force;
}

matrix_row* M_multiply(matrix_row & q) {
	
	matrix_row* result = new matrix_row;
	result->x0 = M->rows[0].s0 * q.x0 + M->rows[0].s1 * q.x1 
	+ M->rows[0].x0 * q.s0 + M->rows[0].x1 * q.s1;
	
	result->x1 = M->rows[1].s0 * q.x0 + M->rows[1].s1 * q.x1
	+ M->rows[1].x0 * q.s0 + M->rows[1].x1 * q.s1;
	
	result->s0 = M->rows[2].x0 * q.x0 + M->rows[2].x1 * q.x1
	+ M->rows[2].s0 * q.s0 + M->rows[2].s1 * q.s1;
	
	result->s1 = M->rows[3].x0 * q.x0 + M->rows[3].x1 * q.x1
	+ M->rows[3].s0 * q.s0 + M->rows[3].s1 * q.s1;
	
	return result;
}
	
double integrate(const gsl_vector *v, void *params) {
// for pendulum, s0 and x0 are constant
	matrix_row q;
	q.s0 = sys->n[0].s;
	q.s1 = gsl_vector_get(v,0);
	q.x0[0] = sys->n[0].x[0];
	q.x0[1] = sys->n[0].x[1];
	q.x0[2] = sys->n[0].x[2];
	q.x1[0] = gsl_vector_get( v,1);
	q.x1[1] = gsl_vector_get( v,2);
	q.x1[2] = gsl_vector_get( v,3);
	double lambda = gsl_vector_get(v,4);
	
	matrix_row q0;
	q0.s0 = sys->n[0].s;
	q0.s1 = sys->n[1].s;
	q0.x0[0] = sys->n[0].x[0];
	q0.x0[1] = sys->n[0].x[1];
	q0.x0[2] = sys->n[0].x[2];
	q0.x1[0] = sys->n[1].x[0];
	q0.x1[1] = sys->n[1].x[1];
	q0.x1[2] = sys->n[1].x[2];
	
	matrix_row * force = get_force(q.s1-q.s0, q.x0, q.x1, rho);
	(*force) *= deltaT;
	matrix_row * result = M_multiply(q);
	matrix_row * temp = M_multiply(q0);
	(*temp) += (*force);
	(*temp) *= -1.0; //right hand side Mqdot + hf
	(*result) += (*temp); //Mq - Mq0 - hf
	
	(*temp) = G->rows[0]; //1 row for now, otherwise Gxqdot is an array

	float Gxqdot = temp->x0 * q.x0 + temp->x1 * q.x1 + temp->s0 * q.s0 + temp->s1 * q.s1;
	// G dot qdot = 0
	
	(*temp) *= lambda;
	(*result) += (*temp); //Mq - Mq0 - hf + lambda*Gt
	
	float res = Gxqdot*Gxqdot + result->x0*result->x0 + result->x1*result->x1 + result->s0*result->s0 + result->s1*result->s1;
	free(temp); free(result); free(force);
	
	return res;
}

double post_stabilize(const gsl_vector *v, void *params) {
	
	matrix_row diff;
	diff.s0 = sys->n[0].s;
	diff.s1 = gsl_vector_get(v, 0);
	diff.x0[0] = sys->n[0].x[0];
	diff.x0[1] = sys->n[0].x[1];
	diff.x0[2] = sys->n[0].x[2];
	diff.x1[0] = gsl_vector_get(v, 1);
	diff.x1[1] = gsl_vector_get(v, 2);
	diff.x1[2] = gsl_vector_get(v, 3);
	double lambda = gsl_vector_get(v, 4);
	
	matrix_row q0;
	q0.s0 = sys->n[0].s;
	q0.s1 = sys->n[1].s;
	q0.x0[0] = sys->n[0].x[0];
	q0.x0[1] = sys->n[0].x[1];
	q0.x0[2] = sys->n[0].x[2];
	q0.x1[0] = sys->n[1].x[0];
	q0.x1[1] = sys->n[1].x[1];
	q0.x1[2] = sys->n[1].x[2];
	
	q0 *= -1.0;
	diff += q0;
	
	matrix_row * result = M_multiply(diff); // M * deltaq
	matrix_row temp = G->rows[0];
	//1 row for now, otherwise Gxqdot is an array
	
	float Gxdq = temp.x0 * diff.x0 + temp.x1 * diff.x1 + temp.s0 * diff.s0 + temp.s1 * diff.s1;
	// G dot dq = -g
	
	temp *= lambda;
	(*result) += temp; //M*deltaq + lambda*Gt + G*deltaq
	
	// add g = l^2 - deltas^2
	vec3 temp1 = vec3(sys->n[0].x[0], sys->n[0].x[1], sys->n[0].x[2]);
	double s = gsl_vector_get(v, 0);
	vec3 temp2 = vec3(gsl_vector_get(v,1), gsl_vector_get(v,2), gsl_vector_get(v,3));
	temp2 -= temp1;
	double g_var = radius * radius * temp2*temp2 - (s - sys->n[0].s)*(s - sys->n[0].s);
	
	double res = Gxdq*Gxdq + g_var + result->x0*result->x0 + result->x1*result->x1 + result->s0*result->s0 
	+ result->s1*result->s1;

	free(result);
	return res;		
}

float* Simplex(int func_enum, int num_cons) {
// pass in the function to be minimized and number of constraints
	
	double par[2] = {0.05, 0.5};
	const gsl_multimin_fminimizer_type *T =	gsl_multimin_fminimizer_nmsimplex;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	
	size_t iter = 0;
	int status;
	double size;
	
	/* Starting point */
	x = gsl_vector_alloc (4 + num_cons);
	for (int i=0; i< (4+num_cons); i++) {
	gsl_vector_set (x, i, 0.5);
	}
		
	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (4+num_cons);
	gsl_vector_set_all (ss, 1.0);
	
	/* Initialize method and iterate */
	minex_func.n = 4+num_cons;
	minex_func.f = func_enum==0? integrate:post_stabilize;
	minex_func.params = par;
	
	s = gsl_multimin_fminimizer_alloc (T, 4+num_cons);
	gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
	
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		
		if (status) 
			break;
		
		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-2);
		
		if (status == GSL_SUCCESS)
		{
			printf ("converged to minimum at\n");
		}
		
		printf ("%5d %10.3e %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n", 
				iter,
				gsl_vector_get (s->x, 0), 
				gsl_vector_get (s->x, 1), 
				gsl_vector_get (s->x, 2), 
				gsl_vector_get (s->x, 3), 
				s->fval, size);
	}
	while (status == GSL_CONTINUE && iter < 100);
	
	float * values = new float[4];
	values[0] = gsl_vector_get (s->x, 0); 
	values[1] = gsl_vector_get (s->x, 1); 
	values[2] = gsl_vector_get (s->x, 2); 
	values[3] = gsl_vector_get (s->x, 3); 
	
	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);
	
	return values;	
}
//	
//void takeAStep() { //calls all the functions and evolve one time step
//
//	float* s = Simplex(0,1);
//	sys->n[1].sdot = s[0];
//	sys->n[1].xdot[0] = s[1];
//	sys->n[1].xdot[1] = s[2];
//	sys->n[1].xdot[2] = s[3];
//	
//	sys->n[1].s += sys->n[1].sdot*deltaT;
//	sys->n[1].x[0] += sys->n[1].xdot[0]*deltaT;
//	sys->n[1].x[1] += sys->n[1].xdot[1]*deltaT;
//	sys->n[1].x[2] += sys->n[1].xdot[2]*deltaT;
//	free(s);
//	
//	s = Simplex(1, 1);
//	sys->n[1].s += s[0];
//	sys->n[1].x[0] += s[1];
//	sys->n[1].x[1] += s[2];
//	sys->n[1].x[2] += s[3];
//	free(s);
//}
//
//void Eulerstep(nodeSystem_t sys) {
//	//time_integrate
//	
//	vec3 x0 = vec3(sys->n[0].x[0], sys->n[0].x[1], sys->n[0].x[2]);
//	vec3 x1 = vec3(sys->n[1].x[0], sys->n[1].x[1], sys->n[1].x[2]);
//	matrix_row_t f = get_force(sys->deltas, x0, x1, sys->rho);
//	vec3 qdot = vec3(sys->n[0].xdot[0], sys->n[0].xdot[1], sys->n[0].xdot[2]);
//	f *= sys->deltat;
//	f += M * qdot; // 1 x 8
//
//	MatrixXd result = MatrixXd::Zero(4 + cons_row, cons_row);
//	result.block(0,0,1,4) = f;
//	
//	MatrixXd eqns = MatrixXd::Zero(4 + cons_row,8 + cons_row); 
//	eqns.block(0,0,4,8) = M;
//	eqns.block(4,0,cons_row,8) = G;
//	eqns.block(0,8,8,cons_row) = G.transpose();
//	
//	JacobiSVD<MatrixXd> svd(eqns, ComputeThinU | ComputeThinV);
//	VectorXd newqdot() = svd.solve(result);
//	
//	sys->n[1].xdot[0] = newqdot(0);
//	sys->n[1].xdot[1] = newqdot(1);
//	sys->n[1].xdot[2] = newqdot(2);
//	sys->n[1].sdot = newqdot(3);
//	
//	Vector4d newq(sys->n[1].x[0], sys->n[1].x[1], sys->n[1].x[2], sys->n[1].s);
//	newq += sys->deltat * newqdot;
//	sys->n[1].x[0] = newq(0);
//	sys->n[1].x[1] = newq(1);
//	sys->n[1].x[2] = newq(2);
//	sys->n[1].s = newq(3);
//	
//	//post_stabilize 
//	MatrixXd post_result = MatrixXd::Zero(4+cons_row, cons_row);
//	post_result(4,0) = -(radius*radius*(deltax[0]*deltax[0] + deltax[1]*deltax[1] + deltax[2]*deltax[2]) - deltas*deltas);
//	
//	VectorXd deltaq = svd.solve(post_result);
//	sys->n[1].x[0] += deltaq(0,0);
//	sys->n[1].x[1] += deltaq(1,0);
//	sys->n[1].x[2] += deltaq(2,0);
//	sys->n[1].s += deltaq(3,0);
//}
//

//void animate() {
// outputs jpegs or pngs to be compiled by mpeg animator


//}


int main() {
   sys = new nodeSystem;
   sys->num = 2;
   sys->n = new node[sys->num];
   sys->n[0].x = new float[3];
   sys->n[0].x[0] = 0.0;
   sys->n[0].x[1] = 0.0;
   sys->n[0].x[2] = 0.0;
   sys->n[1].x = new float[3];
   sys->n[1].x[0] = 0.0;
   sys->n[1].x[1] = 0.0;
   sys->n[1].x[2] = 0.0;
   sys->n[0].s = 0;
   sys->n[1].s = 1;
   sys->n[0].sdot = 0;
   sys->n[1].sdot = 1;
   sys->n[0].xdot = new float[3];
   sys->n[0].xdot[0] = 0.0;
   sys->n[0].xdot[1] = 0.0;
   sys->n[0].xdot[2] = 0.0;
   sys->n[1].xdot = new float[3];
   sys->n[1].xdot[0] = 0.0;
   sys->n[1].xdot[1] = 0.0;
   sys->n[1].xdot[2] = 0.0;
	float* s = Simplex(0,1);
	printf("%f %f %f %f\n", s[0], s[1], s[2], s[3]);
		free(s);
		
		s = Simplex(1, 1);
	printf("%f %f %f %f\n", s[0], s[1], s[2], s[3]);
		free(s);
	return 0;
}
