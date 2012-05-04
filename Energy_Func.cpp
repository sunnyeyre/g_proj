/*
 *  Energy_Func.cpp
 *  
 *
 */

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
#include "Energy_Func.h"
#include "linalg.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_min.h"
#include <Eigen/Dense>
#include <Eigen/Core>

#define PI 3.14159265
#define g 9.80665
#define k 0.2 //spring constant and length of string
#define mass 0.2
#define length 0.1

using namespace Eigen;
using namespace std;
extern float dt;


circle_constraint::circle_constraint(float a, float b, float c, node* child){
    
    l = a;
    theta = b;
    z0 = c;
    type = CIRCLE;
    nodes = child;
}

void circle_constraint::get_world_coordinates(Vector3f& tmp) {

    tmp[0] = l*sin(theta-PI/2.0);
    tmp[1] = l*cos(theta-PI/2.0);
    tmp[2] = z0;
    
}

bead::bead(float r, constraint* c, Vector3f& x) {
    
    radius = r;
    cons = c;
    type = BEAD;
    world_x = new Vector3f(x);
    next = NULL;
}

void bead::display() {
    
    
    glTranslatef(world_x->operator[](0), world_x->operator[](1), world_x->operator[](2));
    
    glPushMatrix();

    if (next != NULL )
        next->display();
}

bead * bead::update(Vector3f& x) {
    
    if(next != NULL) {
        next = next->update(x);
    }
    
    if (cons == NULL)
        return this;
    else {
        // calls constraint node to update world coordinates
    Vector3f temp = Vector3f(x);
 
    bead * new_bead = new bead(radius, cons, temp);
    return new_bead;
    }
}

L0node::L0node(double m, Vector3f& w_x) {

    cons_node = NULL;
    material_coordinate = m;
    world_x = new Vector3f(w_x);
    left = NULL;
    right = NULL;
    x = NULL;
    
}

void node::display() {
    
    if (cons_node == NULL || cons_node->type == CIRCLE) {
        
        if(left == NULL)
            glVertex3f(world_x->operator[](0), world_x->operator[](1), world_x->operator[](2));
        
        else {
            glVertex3f(world_x->operator[](0)-left->world_x->operator[](0), world_x->operator[](1)
                     - left->world_x->operator[](1), world_x->operator[](2)-left->world_x->operator[](2));
        }
        
        if (right != NULL)
            right->display();

        return;
    }
    printf("sorry! haven't implemented circular constraints yet...\n");
}

L0node* L0node::update() {
    
    L0node * new_node = new L0node(material_coordinate, *world_x);
    new_node->func = this->func;
    
    if(right != NULL) {
        new_node->right = right->update();
        new_node->right->left = new_node;
    }
    
    return new_node;
}

L1node::L1node(double m, double x_tilda, double x_tilda_dot, Vector3f& w_x) {
    
    cons_node = NULL;
    material_coordinate = m;
    left = NULL;
    right = NULL;
    x = new union type_of_parameter;
    x->d = x_tilda;
    x_dot = new union type_of_parameter;
    x_dot->d = x_tilda_dot;
    world_x = new Vector3f(w_x);
    isEnode = false;
}

L1node::~L1node() {
    
    delete x;
    delete x_dot;
    delete world_x;
    
    if(cons_node != NULL)
        delete cons_node;
}
/*
L1node* pendulum_update(L1node * old_node) {
    
    L1node * new_node = new L1node(old_node->material_coordinate, old_node->x->d, old_node->x_dot->d, *(old_node->world_x));
    new_node->x_dot->d -= (g/dynamic_cast<circle_constraint *> (old_node->cons_node)->l) * sin(old_node->x->d) * dt;
    new_node->x->d += dt * new_node->x_dot->d;
    
    new_node->cons_node = new circle_constraint(dynamic_cast<circle_constraint *> (old_node->cons_node)->l, new_node->x->d, old_node->world_x->operator[](2), new_node);
    dynamic_cast<circle_constraint *> (new_node->cons_node)->get_world_coordinates(*(new_node->world_x)); //updates world coordinates
    
    return new_node;
}
*/
node * linear_search(float s, node * old_right) {
// given an s and its old rightmost node, find interval [s0, s1] and return Lnode with s0
   
   if( s <= old_right->material_coordinate) {
      if( s >= old_right->left->material_coordinate)
         return old_right->left;
      
      node * ptr = old_right->left;
      while (ptr->left != NULL) {
         if(ptr->isEnode) 
            ptr = ptr->left;
         else {

            if( s >= ptr->left->material_coordinate)
               return ptr->left;
            ptr = ptr->left;
         }
      }
   }

   else {
      node * ptr = old_right->right;
      while (ptr != NULL) {
         if(ptr->isEnode)
            ptr = ptr->right;
         else {

            if( s <= ptr->material_coordinate)
               return ptr->left;
            ptr = ptr->right;
         }
      }
   }
}

L1node* L1node::update() { //make updates to constrained variables
    
    L1node* (*ptrFunc)(L1node*) = (L1node*(*)(L1node*))(this->func);
    
    L1node * new_node = (*ptrFunc)(this);
    new_node->func = this->func;
    
    if(right != NULL) {
        new_node->right = right->update();
        new_node->right->left = this;
    }
        
    return new_node;
    
}

L2node::L2node(double m, Vector2f& x_tilda, Vector2f& x_tilda_dot, Vector3f& w_x) {
    
    cons_node = NULL;
    material_coordinate = m;
    left = NULL;
    right = NULL;
    x = new union type_of_parameter;
    x->v2 = new Vector2f(x_tilda[0], x_tilda[1]);
    x_dot = new union type_of_parameter;
    x_dot->v2 = new Vector2f(x_tilda_dot[0], x_tilda_dot[1]);
    world_x = new Vector3f(w_x);
    isEnode = false;

}

L2node::~L2node() {
    
    delete x->v2;
    delete x;
    delete world_x;
    
    if ( cons_node != NULL)
        delete cons_node;

}

L2node * L2node::update() {
    
    return this;
}

L3node::L3node(double m, Vector3f& x_tilda, Vector3f& x_tilda_dot, Vector3f& w_x) {
    cons_node = NULL;
    material_coordinate = m;
    left = NULL;
    right = NULL;
    x = new union type_of_parameter;
    x->v3 = new Vector3f(x_tilda[0], x_tilda[1], x_tilda[2]);
    x_dot = new union type_of_parameter;
    x_dot->v3 = new Vector3f(x_tilda_dot[0], x_tilda_dot[1], x_tilda_dot[2]);
    world_x = new Vector3f(w_x);
    isEnode = false;
    rho = 1.0;

}

L3node::~L3node() {
    
    delete x->v3;
    delete x;
    delete world_x;
    if ( cons_node != NULL)
        delete cons_node;

}

void L3node::left_force_accumulate(Vector4f& force_right) {
    
        //TODO : check if this is according to the spring equation in Lec23's spring system and gravity
	
	if(left != NULL) {
        Vector3f temp(*(left->world_x) - *world_x);
        Vector4f x_left(temp[0], temp[1], temp[2], left->material_coordinate - material_coordinate);
        x_left.normalize();
			
        Vector4f force_left = x_left * k * sqrt( (*world_x - *(left->world_x)).dot(*world_x - *(left->world_x)) - length*length);
        Vector4f result = force_left - force_right;
	 
		*(this->force) = result;
		
		left_force_accumulate(force_left);
	}
}

void L3node::right_force_accumulate(Vector4f& force_left) {

	if(right != NULL) {
        Vector3f temp(*(right->world_x) - *world_x);

		Vector4f x_right(temp[0], temp[1], temp[2], right->material_coordinate - material_coordinate);
		x_right.normalize();
	 
		Vector4f force_right = x_right * k * sqrt( (*world_x - *(right->world_x)).dot(*world_x - *(right->world_x)) - length*length);
		
        *(this->force) = force_right - force_left;
		
		right_force_accumulate(force_right);
	}	
}

void E3node::force_accumulate() {
	
	node * s1_node = linear_search(material_coordinate, right_s);
	float alpha = (s1_node->left->material_coordinate - material_coordinate)
	/(s1_node->left->material_coordinate - s1_node->material_coordinate);
	
	// s can also have an acceleration and force, even though its just a mesh
    Vector3f tmp_left(*(s1_node->left->world_x) - *world_x);
	Vector4f x_left(tmp_left[0], tmp_left[1], tmp_left[2],s1_node->left->material_coordinate - material_coordinate);
	x_left.normalize();
	
	Vector4f force_left  = x_left * (1-alpha) * k * sqrt( (*world_x - *(s1_node->left->world_x)).dot(*world_x - *(s1_node->left->world_x)) - length*length*(1-alpha)*(1-alpha));
	
    Vector3f tmp_right(*(s1_node->world_x) - *world_x);
    Vector4f x_right(tmp_right[0], tmp_right[1], tmp_right[2], s1_node->material_coordinate - material_coordinate);
	 x_right.normalize();
	 
	 Vector4f force_right = x_right * alpha * k * sqrt( (*world_x - *(s1_node->world_x)).dot(*world_x - *(s1_node->world_x)) - length*length*alpha*alpha);
	 
//	 Vector4f gravity = (-1.0) * x_left * mass * g * abs((-1.0 * x_left).dot(Vector4f(0,-1.0,0))) - x_right * mass * g *
//	 abs((-1.0 *x_right).dot(Vector4f(0,-1,0,0))); ***** TODO : what is this supposed to do?
    Vector4f gravity = x_left - x_right;
	 
	 Vector4f result = gravity + force_left + force_right;
    *(this->force) = result;
	   
	 dynamic_cast<L3node *>(s1_node->left)->left_force_accumulate(force_left);
	 dynamic_cast<L3node *>(s1_node)->right_force_accumulate(force_right);
}

E3node* E3node::update() {
    return NULL;
}

E3node* E3node::update(node * left) {
        // 1 find left and right Lagrangian nodes using material_coordinate
        // 2 call force_accumulate, plug in correct force equation
        // 3 solve for qdot for both left and right q, then use alpha to update E3node
        // 4 call Lagrangian nodes' update
    this->force_accumulate();
    
 //time_integrate step
 
 // input M matrix
	MatrixXf M(8,8);
	float deltas = left->right->material_coordinate - left->material_coordinate;
	Vector3f deltax = *(left->right->world_x) - *(left->world_x);
	float mss = deltax.dot(deltax) / deltas;
	
	M << 2 * deltas, deltas, -2 * deltax[0], -2*deltax[1], -2*deltax[2], -deltax[0], -deltax[1], -deltax[2],  
		deltas, 2*deltas, -deltax[0], -deltax[1], -deltax[2], -2*deltax[0], -2*deltax[1], -2*deltax[2],
		-2*deltax[0], -2*deltax[1], -2*deltax[2], -deltax[0], -deltax[1], -deltax[2], 2*mss, mss,
		-deltax[0], -deltax[1], -deltax[2], -2*deltax[0], -2*deltax[1], -2*deltax[2], mss, 2*mss;
		
// input right hand matrix, M*q0 + dt*force
	Vector4f f, qdot0, rhs;
	qdot0 << world_x_dot[0], world_x_dot[1], world_x_dot[2], sdot;
	f = *force;
	
	rhs = M * qdot0 + dt * f;
	VectorXf rhs_full(8);
	rhs_full << rhs(0), rhs(1), rhs(2), rhs(3), 0, 0, 0, 0;
	
	VectorXf qdot = M.inverse() * rhs_full;

	E3node * new_node = new E3node(*this);
	*(new_node->world_x_dot) = qdot;
	
	Vector3f myqdot(qdot(0), qdot(1), qdot(2));
	*(new_node->world_x) = *world_x + dt*myqdot;
	new_node->material_coordinate = material_coordinate + dt*qdot(3);
	
    return new_node;
}

L3node* L3node::update() {
    return NULL;
}

L3node* L3node::update(node * left) {
  
 //time_integrate step
 
 // input M matrix
	MatrixXf M(8,8);
	float deltas = left->right->material_coordinate - left->material_coordinate;
	Vector3f deltax = *(left->right->world_x) - *(left->world_x);
	float mss = deltax.dot(deltax) / deltas;
	
	M << 2 * deltas, deltas, -2 * deltax[0], -2*deltax[1], -2*deltax[2], -deltax[0], -deltax[1], -deltax[2],  
		deltas, 2*deltas, -deltax[0], -deltax[1], -deltax[2], -2*deltax[0], -2*deltax[1], -2*deltax[2],
		-2*deltax[0], -2*deltax[1], -2*deltax[2], -deltax[0], -deltax[1], -deltax[2], 2*mss, mss,
		-deltax[0], -deltax[1], -deltax[2], -2*deltax[0], -2*deltax[1], -2*deltax[2], mss, 2*mss;
		
// input right hand matrix, M*q0 + dt*force
	Vector4f force, qdot0, rhs;
	qdot0 << world_x_dot[0], world_x_dot[1], world_x_dot[2], 0;
	force << force[0], force[1], force[2], force[3];
	
	rhs = M * qdot0 + dt * force;
	VectorXf rhs_full(8);
	rhs_full << rhs(0), rhs(1), rhs(2), rhs(3), 0, 0, 0, 0;
	
	VectorXf qdot = M.inverse() * rhs_full;

	L3node * new_node = new L3node(*this);
	*(new_node->world_x_dot) = qdot;
	
	Vector3f myqdot = Vector3f(qdot(0), qdot(1), qdot(2));
	*(new_node->world_x) = *world_x + dt*myqdot;
	
    return new_node;
}

E0node::E0node(double m, Vector3f& w_x) {
    
    world_x = new Vector3f(w_x);
    isEnode = true;
    
}

E0node::~E0node() {
    
    delete world_x;

    if ( cons_node != NULL)
        delete cons_node;
    
}

E0node* E0node::update() {
    
    return this;
}

E1node::E1node(double m, double x_tilda, double x_tilda_dot, Vector3f& w_x) {
 
    world_x = new Vector3f(w_x);
    isEnode = true;

}

E1node::~E1node() {
    
    delete world_x;

    if ( cons_node != NULL)
        delete cons_node;
    
}

E1node * E1node::update() {
    
    return this;
}

E1node * E1node::update(node* q0, node* q1) {// applies gravity and interpolates E1node
   return this;
    
}

E2node::E2node(double m, Vector2f& x_tilda, Vector2f& x_tilda_dot, Vector3f& w_x) {
    
    world_x = new Vector3f(w_x);
    isEnode = true;

}

E2node::~E2node() {
    
    delete world_x;

    if ( cons_node != NULL)
        delete cons_node;
    
}

E2node * E2node::update() {
    
    return this;
}

E3node::E3node(double m, Vector3f& x_tilda, Vector3f& x_tilda_dot, Vector3f& w_x) {
    world_x = new Vector3f(w_x);
    isEnode = true;
    material_coordinate = m;
    
}

E3node::~E3node() {
    
    delete world_x;
    
    if ( cons_node != NULL)
        delete cons_node;
    
}