/*
 *  Energy_Func.h
 *  
 *
 */

#ifndef ENERGY_FUNC_H
#define ENERGY_FUNC_H

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

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#endif

#include <OpenGL/gl.h>
#include "glui.h"
#include <time.h>
#include "algebra3.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>

#define PI 3.14159265
#define g 9.80665

using namespace std;
using namespace Eigen;

using Eigen::Vector3f;
using Eigen::Vector2f;

class node; //forward declaration
class L1node;
typedef L1node * (*fptr)(L1node *);

    // check proper enum syntax
enum type_of_constraint {CIRCLE, RIGID, FIXED, CYLINDER, NORMAL};
enum type_of_display {BEAD};

union type_of_parameter {
    
    double d;
    Vector2f* v2;
    Vector3f* v3;
};

class constraint { // for converting from constrained coords to world coords, no state copies

public:
    enum type_of_constraint type;
    virtual void get_world_coordinates(Vector3f& tmp)=0;
    node * nodes;
};

class circle_constraint : public constraint {
        //draws a circle about the z plane
public:
    circle_constraint(float a, float b, float c, node* child);
    void get_world_coordinates(Vector3f& tmp);
    float l;
    float theta;
    float z0;
};

class node {

public :
    
    void display(); // for drawing a line from 0 to L
    virtual node * update() =0;

        //    constraint * cons_node;
    enum type_of_constraint cons;
    double material_coordinate; // this is s, fixed for L nodes and variable to E nodes
                                //    double rho; // density of segment
    Vector3f * world_x;
    Vector3f * world_x_dot;
    node * left;
    node * right;
    union type_of_parameter * x;
    union type_of_parameter * x_dot;
	
    float rho;
    fptr func; // function for force accumulation
    bool isEnode; //for checking if the node is a special Enode
    node * linear_search(float s, node* rootNode); //finds node to right of material coordinate
};

class display_node { //only for displaying, no state copies
    
public:
    
    virtual void display()=0; //for displaying beads, tubes, etc.
    virtual display_node* update(Vector3f& x) = 0;
    constraint * cons;
    Vector3f * world_x;
    enum type_of_display type;
    display_node * next;
};

class bead : public display_node {

public:
    bead(float r, constraint* c,  Vector3f& x) ;
    void display();
    bead * update(Vector3f& x);
    float radius;
};
    
class L0node : public node {
public:
    L0node(double m,  Vector3f& w_x) ;
    ~L0node();
    
    L0node* update();

};

class L1node : public node {
public:
    L1node(double m, double x_tilda, double x_tilda_dot,  Vector3f& w_x) ;
    ~L1node();
    
    L1node* update();
    
};

class L2node : public node {
public:
    L2node(double m,  Vector2f & x_tilda,  Vector2f & x_tilda_dot,  Vector3f & w_x);
    ~L2node();
    
    L2node* update();
    
};
class L3node : public node {
public:

    L3node(double m, Vector3f& w_x, Vector3f& w_x_dot);
    ~L3node();
    
    Vector3f * force; // uses force accumulate to add forces
    
    L3node* update();
    L3node* updateLeft(Vector3f& x, Vector3f& xdot, float alpha); //updates and links
    L3node* updateRight(Vector3f& x, Vector3f& xdot, float alpha);
	void left_force_accumulate(Vector3f& x, Vector3f& xdot, float alpha);
    void right_force_accumulate(Vector3f& x, Vector3f& xdot, float alpha);
};
class E0node : public node {
public:
    E0node(double m, Vector3f& w_x);
    ~E0node();
    
    E0node* update();
    
    double sdot;
};
class E1node : public node {
public:
    E1node(double m, double x_tilda, double x_tilda_dot, Vector3f& w_x);
    ~E1node();
    
    E1node* update();
    E1node* update(node * q0, node * q1);
    
	double sdot;
};
class E2node : public node {
public:
    E2node(double m, Vector2f& x_tilda, Vector2f& x_tilda_dot, Vector3f& w_x);
	double sdot;
    ~E2node();
    
    E2node* update();
    
};

class E3node : public node {
public:

    E3node();
    E3node(double m, double mdot, Vector3f& w_x, Vector3f& w_x_dot, node* root);
    ~E3node();
    
    node * rootNode;
    
    float sdot;
        //    Vector3f pos; //4 element vector  
        //    Vector3f velocity; // 4 element vector
    VectorXf force; // 8 element vector
    float alpha; //interpolation parameter
    
    E3node* update();
    E3node* update(node* left);
    void force_accumulate(Vector4f& q0, Vector4f& q1);
	
};

//class matrix_row {
//public:
//	matrix_row & operator *= (const vec3& b);
//	matrix_row & operator *= (float b);
//	matrix_row & operator += (const matrix_row & b);
//	matrix_row & operator = (const matrix_row & b);
//	matrix_row();
//	~matrix_row();
//	matrix_row(const matrix_row & copy);
//	matrix_row(float a, float b, const vec3& c, const vec3& d);
//	
//	float s0;
//	float s1;
//	vec3 x0;
//	vec3 x1;
//private:
//};
//
//typedef struct mat {
//	matrix_row* rows;
//} * mat_t;
//
//extern int cons_row; //rows and columns of constraint by G
//extern mat_t M; 
//extern mat_t G;
//
//typedef struct node{
//	float rho;
//	float* x;
//	float s;
//	float* xdot;
//	float sdot;
//} * node_t;
//
//typedef struct {
//	node_t n;
//	int num;
//} nodeSystem;
//
//extern nodeSystem* sys;
//extern float radius;
//extern float deltaT;
//extern float rho;
//
//void takeAStep();
//void input_M_matrix(float deltas, vec3& deltax, float rho, float X);
//	//inputs and updates M matrix given X, deltas, deltax, and rho
//void input_G_matrix(int cons_row, vec3& deltax, float deltas, float radius);
//
//vec3 get_interpolated_x(float s0, float s1, float s, vec3& x0, vec3& x1);
//vec3 get_interpolated_xdot(float sdot0, float sdot1, vec3& xdot0, vec3& xdot1, vec3& deltax, vec3& deltas, vec3& derivx);
//
#endif

