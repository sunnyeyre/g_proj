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
#include <Eigen/SVD>
#include <Eigen/Cholesky>

#define PI 3.14159265
#define g 9.80665
#define ks 8.0 //spring constant and length of string
#define mass 0.2
#define length 0.1
#define kd 2.0 // damping constant of spring

using namespace Eigen;
using namespace std;
extern float dt;

int NUM = 10; // number of Lagrangian nodes in between 2 rigid nodes
float l0 = 4.0/NUM; // rest length of spring
float maxL = 20 * l0;

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

    material_coordinate = m;
    world_x = new Vector3f(w_x);
    left = NULL;
    right = NULL;
    x = NULL;
    
}

void node::display() {
    
    if (cons == CIRCLE) {
        
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
node * node::linear_search(float s, node* rootNode) {
// given an s and its old rightmost node, find interval [s0, s1] and return Lnode with s0
   
    int i = floor(s * NUM);
    int j = ceil(s * NUM);
    if(abs((float) i - s*NUM) < 0.0001) {
        j = i+1;
    }
    
    node * iterator = rootNode;
    while(abs(iterator->material_coordinate - (float)j/(float)NUM) > 0.0001)
        iterator = iterator->right;
    
    return iterator;

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
    
}

L2node * L2node::update() {
    
    return this;
}

L3node::L3node(double m, Vector3f& w_x, Vector3f& w_x_dot) {
    material_coordinate = m;
    left = NULL;
    right = NULL;
//    x = new union type_of_parameter;
//    x->v3 = new Vector3f(x_tilda[0], x_tilda[1], x_tilda[2]);
//    x_dot = new union type_of_parameter;
//    x_dot->v3 = new Vector3f(x_tilda_dot[0], x_tilda_dot[1], x_tilda_dot[2]);
    world_x = new Vector3f(w_x);
    world_x_dot = new Vector3f(0.0, 0.0, 0.0);
    isEnode = false;
    force = new Vector3f(0.0, 0.0, 0.0);
    rho = 1.0;

}

L3node::~L3node() {
    
    delete x->v3;
    delete x;
    delete world_x;
//    if ( cons_node != NULL)
//        delete cons_node;

}

void L3node::left_force_accumulate(Vector3f& x, Vector3f& xdot, float alpha, const Vector3f& f) { // spring
    //testing: need to put in the right scaling factors, but this equation is for vertical dimension only
    //
    Vector3f dxLeft = x - *world_x;
    dxLeft.normalize();
    *force = dxLeft.dot(f) * f * (1.0-alpha);
    float restLengthRight = l0 * alpha;
    float restLengthLeft = l0;
    
    Vector3f deltaXDotRight = *(world_x_dot) - xdot;
    Vector3f deltaXDotLeft = *(world_x_dot) - *(left->world_x_dot);
    Vector3f deltaXRight = *(world_x) - x;
    Vector3f deltaXLeft = *world_x -*(left->world_x);
    float absDeltaXRight = sqrt(deltaXRight.dot(deltaXRight));
    float absDeltaXLeft = sqrt(deltaXLeft.dot(deltaXLeft));
    Vector3f constant = -(ks * (absDeltaXLeft - restLengthLeft) + kd * deltaXDotLeft.dot(deltaXLeft)/absDeltaXLeft)*deltaXLeft/absDeltaXLeft - (ks * (absDeltaXRight - restLengthRight) + kd * deltaXDotRight.dot(deltaXRight)/absDeltaXRight)*deltaXRight/absDeltaXRight;
   *force += constant;
    
        //     cout << "force" << *force << endl;
}

void L3node::right_force_accumulate(Vector3f& x, Vector3f& xdot, float alpha, const Vector3f& f) {

    Vector3f dxRight = x - *world_x;
    dxRight.normalize();
    *force = dxRight.dot(f) * f * (1.0-alpha);
    float restLengthLeft = l0 * alpha;
    float restLengthRight = l0;
    
    Vector3f deltaXDotLeft = *(world_x_dot) - xdot;
    Vector3f deltaXDotRight = *(world_x_dot) - *(right->world_x_dot);
    Vector3f deltaXLeft = *world_x - x;
    Vector3f deltaXRight = *world_x - *(right->world_x);
    float absDeltaXLeft = sqrt(deltaXLeft.dot(deltaXLeft));
    float absDeltaXRight = sqrt(deltaXRight.dot(deltaXRight));
    Vector3f constant = -(ks * (absDeltaXLeft - restLengthLeft) + kd * deltaXDotLeft.dot(deltaXLeft)/absDeltaXLeft)*deltaXLeft/absDeltaXLeft - (ks * (absDeltaXRight - restLengthRight) + kd * deltaXDotRight.dot(deltaXRight)/absDeltaXRight)*deltaXRight/absDeltaXRight;
   *force += constant;
    
        //   cout << "force" << *force << endl;

}

void E3node::force_accumulate(Vector4f& q0, Vector4f& q1) { //gravity only, resets force every time
	
    Vector3f gravity = Vector3f(0, g, 0);
    float constant = rho / 2.0;
//    float deltas = s1_node->material_coordinate - s0_node->material_coordinate;
//    Vector3f x0_plus_x1 = Vector3f(*(s1_node->world_x) + *(s0_node->world_x));
    float deltas = q1[3] - q0[3];
    Vector4f sumq = q1 + q0;
    Vector3f x0_plus_x1 = Vector3f(sumq[0], sumq[1], sumq[2]);
    
    float g_dot_sum = x0_plus_x1.dot(gravity);    
    
    force << gravity.x()*deltas, gravity.y()*deltas, gravity.z()*deltas, gravity.x()*deltas, gravity.y()*deltas, gravity.z()*deltas, g_dot_sum * (-1.0), g_dot_sum;
//    force << g* deltas, g * deltas, 
//            -g_dot_sum, g_dot_sum;
    force *= -constant;
}

void L3node::force_accumulate() { //gravity only, resets force every time
	
    Vector3f gravity = Vector3f(0, -g, 0);
    *force << gravity.x(), gravity.y(), gravity.z();
}

L3node* L3node::update() {

	node * s1_node = linear_search(material_coordinate, rootNode);
    node * s0_node = s1_node->left;

    Vector3f deltax = *(s0_node->world_x) - *world_x;
    force_accumulate();
    if(abs(deltax.dot(deltax)) > maxL)
       *force = Vector3f(0.0,0.0,0.0);

    Vector3f xDotNew = *force * dt/ rho + *world_x_dot;
    Vector3f xNew = xDotNew*dt + *world_x;

    L3node * newNode = new L3node(material_coordinate, xNew, xDotNew);

    L3node* leftPart = dynamic_cast<L3node *> (s0_node)->updateLeft(xNew, xDotNew, material_coordinate, *force);
    L3node* rightPart = dynamic_cast<L3node *> (s1_node)->updateRight(xNew, xDotNew, 1.0-material_coordinate, *force);
    leftPart->right = rightPart;
    rightPart->left = leftPart;
    node* root = leftPart;
    while(root->left != NULL)
        root =  root->left;
    newNode->rootNode = root;
    
    return newNode;
}

E3node* E3node::update() {

	node * s1_node = linear_search(material_coordinate, rootNode);
    node * s0_node = s1_node->left;
    
    Vector4f q0(s0_node->world_x->x(), s0_node->world_x->y(), s0_node->world_x->z(), s0_node->material_coordinate);
    Vector4f q1(s1_node->world_x->x(), s1_node->world_x->y(), s1_node->world_x->z(), s1_node->material_coordinate);
    
    force_accumulate(q0, q1);

        // input M matrix
	MatrixXf M(8,4);
//	float deltas = right->material_coordinate - left->material_coordinate;
//	Vector3f deltax = *(right->world_x) - *(left->world_x);
    Vector4f deltaq = q1 - q0;
    float deltas = q1[3] - q0[3];
    Vector3f deltax = Vector3f(deltaq[0], deltaq[1], deltaq[2]);
	float mss = deltax.dot(deltax) / deltas;
	
	M << 2 * deltas, deltas, -2 * deltax[0], -2*deltax[1], -2*deltax[2], -deltax[0], -deltax[1], -deltax[2],  
    deltas, 2*deltas, -deltax[0], -deltax[1], -deltax[2], -2*deltax[0], -2*deltax[1], -2*deltax[2],
    -2*deltax[0], -2*deltax[1], -2*deltax[2], -deltax[0], -deltax[1], -deltax[2], 2*mss, mss,
    -deltax[0], -deltax[1], -deltax[2], -2*deltax[0], -2*deltax[1], -2*deltax[2], mss, 2*mss;
    
    M /= rho/6.0;

        // input right hand matrix, M*q0 + dt*force

    float alpha = (material_coordinate - q0[3])/(q1[3]-q0[3]);
    Vector3f pos = (1.0-alpha)* *(s0_node->world_x) + alpha * *(s1_node->world_x);
    Vector3f velocity = (1.0-alpha) * *(s0_node->world_x_dot) + alpha * *(s1_node->world_x_dot);
    Vector4f q = Vector4f(pos[0], pos[1], pos[2], material_coordinate);
    Vector4f qdot = Vector4f(velocity[0], velocity[1], velocity[2], sdot);
    
    VectorXf rhs = M * qdot + dt * force;
    
        // inverting the M matrix using JacobiSVD
    JacobiSVD<MatrixXf> svd(M, ComputeThinU | ComputeThinV);
        // save the result in qDotNew
    Vector4f qDotNew = svd.solve(rhs);
    
    cout << "qDotNew : " << qDotNew << endl;
    
    Vector4f qNew = q + dt * qDotNew;
    cout << "qNew : " << qNew << endl;
    Vector3f xNew = Vector3f(qNew[0], qNew[1], qNew[2]);
    Vector3f xDotNew = Vector3f(qDotNew[0], qDotNew[1], qDotNew[2]);

        // make a new E3node for the new state with the udpated information
    
    E3node * newNode = new E3node(qNew[3], qDotNew[3], xNew, xDotNew, rootNode);
//	newNode->velocity = M.inverse() * rhs; //qdot update
//    newNode->pos = pos + dt * velocity; //position update

//    Vector3f x0dot = Vector3f(velocity[0], velocity[1], velocity[2]);
//    Vector3f x1dot = Vector3f(velocity[3], velocity[4], velocity[5]);
//    Vector3f x0 = Vector3f(pos[0], pos[1], pos[2]);
//    Vector3f x1 = Vector3f(pos[3], pos[4], pos[5]);
//        //post - stasbilization, TODO
//    *(newNode->world_x) = (1.0-alpha) * x0 + alpha*x1;
    
    L3node* leftPart = dynamic_cast<L3node *> (s0_node)->updateLeft(xNew, xDotNew, alpha, force);
    L3node* rightPart = dynamic_cast<L3node *> (s1_node)->updateRight(xNew, xDotNew, 1.0-alpha,force);
    leftPart->right = rightPart;
    rightPart->left = leftPart;
    node* root = leftPart;
    while(root->left != NULL)
        root =  root->left;
    newNode->rootNode = root;
    
    return newNode;
}

L3node* L3node::updateLeft(Vector3f& x, Vector3f& xdot, float alpha, const Vector3f& f) {
    
    if(cons != RIGID) {
        left_force_accumulate(x,xdot, alpha, f);
        
        Vector3f l = *world_x - x;
        Vector3f acceleration;
        if(abs(l.dot(l)) >= maxL)
           acceleration = Vector3f(0.0,0.0,0.0);
        else
            acceleration = *force/ (rho * sqrt(l.dot(l)));
        
        Vector3f newXDot = acceleration * dt + *world_x_dot;
        Vector3f newX = newXDot * dt + *world_x;
        L3node * newNode = new L3node(material_coordinate, newX, newXDot);
        //newNode->left = dynamic_cast<L3node *>(left)->updateLeft(newX, newXDot, 1.0, *force);
        newNode->left = dynamic_cast<L3node *>(left)->updateLeft(newX, newXDot, 1.0, f);
        newNode->left->right = newNode;
        newNode->cons = NORMAL;
        return newNode;
    }
    else {
        L3node * newNode = new L3node(material_coordinate, *world_x, *world_x_dot);
        newNode->left = NULL;
        newNode->cons = RIGID;
        return newNode;
    }
}

L3node* L3node::updateRight(Vector3f& x, Vector3f& xdot, float alpha,const Vector3f& f) {
    
    if(cons != RIGID) {
        right_force_accumulate(x,xdot, alpha, f);

        Vector3f l = *world_x - x;
        Vector3f acceleration;
        if(abs(l.dot(l)) >= maxL)
           acceleration = Vector3f(0.0,0.0,0.0);
        else
            acceleration = *force/ (rho * sqrt(l.dot(l)));
  
        Vector3f newXDot = acceleration * dt + *world_x_dot;
        Vector3f newX = newXDot * dt + *world_x;
        L3node * newNode = new L3node(material_coordinate, newX, newXDot);
        //newNode->right = dynamic_cast<L3node *>(right)->updateRight(newX, newXDot, 1.0, *force);
        newNode->right = dynamic_cast<L3node *>(right)->updateRight(newX, newXDot, 1.0, f);
        newNode->right->left = newNode;
        newNode->cons = NORMAL;
        return newNode;
    }
    else {
        L3node * newNode = new L3node(material_coordinate, *(world_x), *(world_x_dot));
        newNode->right = NULL;
        newNode->cons = RIGID;
        return newNode;
    }
}

E0node::E0node(double m, Vector3f& w_x) {
    
    world_x = new Vector3f(w_x);
    isEnode = true;
    
}

E0node::~E0node() {
    
    delete world_x;
//
//    if ( cons_node != NULL)
//        delete cons_node;
    
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

//    if ( cons_node != NULL)
//        delete cons_node;
    
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

//    if ( cons_node != NULL)
//        delete cons_node;
    
}

E2node * E2node::update() {
    
    return this;
}

E3node::E3node() {
    
    world_x = new Vector3f();
    world_x_dot = new Vector3f(0.0,0.0,0.0);
    isEnode = true;
//    pos = Vector4f(0.0,0,0,0);
//    velocity = Vector4f(0,0,0,0);
    force = VectorXf::Zero(8);
    
    
}

E3node::E3node(double m, double mdot, Vector3f& w_x, Vector3f& w_x_dot, node * root) {

    world_x = new Vector3f(w_x);
    world_x_dot = new Vector3f(w_x_dot);
    isEnode = true;
    material_coordinate = m;
    sdot = mdot;
    
//    node * right_of_this = this->linear_search(material_coordinate, root);
//    node * left_of_this = right_of_this->left;
    
//    pos = Vector4f(0,0,0,0);
//    pos << left_of_this->world_x->x(), left_of_this->world_x->y(), left_of_this->world_x->z(),
//            right_of_this->world_x->x(), right_of_this->world_x->y(), right_of_this->world_x->z(),
//    left_of_this->material_coordinate, right_of_this->material_coordinate;
    
    
//    velocity = Vector4f(0,0,0,0);
//    velocity << left_of_this->world_x_dot->x(), left_of_this->world_x_dot->y(), left_of_this->world_x_dot->z(),
//            right_of_this->world_x_dot->x(), right_of_this->world_x_dot->y(), right_of_this->world_x_dot->z(),
//    0, 0;
    
    force = VectorXf::Zero(8);
    
}

E3node::~E3node() {
    
    delete world_x;
    
//    if ( cons_node != NULL)
//        delete cons_node;
    
}
