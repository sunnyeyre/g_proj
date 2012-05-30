#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "algebra3.h"
#include "FreeImage.h"
#include "Energy_Func.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/LU>

#ifdef OSX
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#include <OpenGL/gl.h>
#include "glui.h"
#else
#include <GL/glut.h>
#include <GL/glu.h>
#include <OpenGL/gl.h>
#include "glui.h"
#endif

#include <time.h>
#include <math.h>

#define PI 3.14159265
#define g 9.80665 // acceleration due to gravity
#define l 2 //length of pendulum
#define ks 8.0 //spring constant and length of string
#define kd 2.0 // damping constant of spring

using namespace Eigen;
using namespace std;

//******************************************************************
//  Some Classes
//******************************************************************

class Viewport {
public:
	int w,h;
};

class State {
public:

    node * root;
    node * e1;
    display_node * d_root;
    State * next;
};

Viewport viewport;
State* initial_state;
State* present_state;

/* some config variables for frames/sec, size of node, init rot angle, etc.*/
float radius = 0.05;
float win_x, win_y, win_z;
float dt = 1.0/80.0;
extern float l0;
int main_window;
extern int NUM; //number of Lagrangian nodes in between the 2 rigid nodes


int iter = 0;


/***********************************************************************/

extern L1node* pendulum_update(L1node * old_node);
extern L3node* string_update(L3node * old_node);
extern E3node* mass_node(E3node * old_node);
extern node* linear_search(float s, node* rootNode);

void advanceT(int value) { //jump to the next state every 5 ms
    
    present_state = present_state->next;
    glutTimerFunc(10, advanceT, 0);
}

State * takeastep(State * previous_state) { // calculates the next state in the simulation
                                            // this is where the node::update funcs are called    

   State * state = new State;

   float deltas;
   node * ptr = previous_state->root->right;
   MatrixXf M = MatrixXf::Zero(3*(NUM-1)+1, 3*(NUM-1)+1); //genralized mass matrix for the entire system
   for(int i=0; i<3*(NUM-1); i+=3) {

      // diagonal entryies = 2*deltas, off diagonal = deltas
         deltas = ptr->right->material_coordinate - ptr->material_coordinate;
         M(i,i) += 2*deltas;
         M(i,i+1) += 2*deltas;
         M(i,i+2) += 2*deltas;
         M(i+1,i) += 2*deltas;
         M(i+2,i) += 2*deltas;
         M(i+1,i+1) += 2*deltas;
         M(i+1,i+2) += 2*deltas;
         M(i+1,i+2) += 2*deltas;
         M(i+2,i+2) += 2*deltas;

         if(i == 3*(NUM-2))
            break;

         M(i+3, i) += deltas;
         M(i+3, i+1) += deltas;
         M(i+3, i+2) += deltas;
         M(i+4, i) += deltas;
         M(i+5, i) += deltas;
         M(i+4, i+1) += deltas;
         M(i+4, i+1) += deltas;
         M(i+4, i+1) += deltas;
         M(i+5, i+2) += deltas;

         M(i, i+3) = deltas;
         M(i+1, i+3) = deltas;
         M(i+2, i+3) = deltas;
         M(i, i+4) = deltas;
         M(i, i+5) = deltas;
         M(i+1, i+4) = deltas;
         M(i+1, i+4) = deltas;
         M(i+1, i+4) = deltas;
         M(i+2, i+5) = deltas;

         M(i+3, i+3) += 2*deltas;
         M(i+3, i+4) += 2*deltas;
         M(i+3, i+5) += 2*deltas;
         M(i+4, i+3) += 2*deltas;
         M(i+5, i+3) += 2*deltas;
         M(i+4, i+4) += 2*deltas;
         M(i+4, i+5) += 2*deltas;
         M(i+5, i+4) += 2*deltas;
         M(i+5, i+5) += 2*deltas;

         ptr = ptr->right;
         i+=3;
   }
   
   node * right_of_e0 = linear_search(previous_state->e1->material_coordinate, previous_state->root);
   
   float alpha = (previous_state->e1->material_coordinate - right_of_e0->left->material_coordinate)/(right_of_e0->material_coordinate - right_of_e0->left->material_coordinate);
   Vector3f deltax = *(right_of_e0->world_x)-*(right_of_e0->left->world_x);
   float mss = deltax.dot(deltax)/(right_of_e0->material_coordinate - right_of_e0->left->material_coordinate);

   M(3*(NUM-1), 3*(NUM-1)) = 2*mss;
   M(3*(NUM-1)-3, 3*(NUM-1)) += -deltax.x();
   M(3*(NUM-1)-2, 3*(NUM-1)) += -deltax.y();
   M(3*(NUM-1)-1, 3*(NUM-1)) += -deltax.z();

   M(3*(NUM-1), 3*(NUM-1)-3) += -deltax.x();
   M(3*(NUM-1), 3*(NUM-1)-2) += -deltax.y();
   M(3*(NUM-1), 3*(NUM-1)-1) += -deltax.z();

   // force term : taking into account of gravity and spring potential energy
   VectorXf f = VectorXf(3*(NUM-1)+1);
   VectorXf qOldDot = VectorXf(3*(NUM-1)+1);
   VectorXf qOld = VectorXf(3*(NUM-1)+1);
   ptr = previous_state->root->right;
   Vector3f deltaXDotRight, deltaXDotLeft, deltaXRight, deltaXLeft, constant;
   float absDeltaXRight, absDeltaXLeft, restLengthRight, restLengthLeft;

   for(int i=0; i<3*(NUM-1); i+=3 ) {

      //putting in old values for Euler integration
      qOldDot[i] = ptr->world_x_dot->x();
      qOldDot[i+1] = ptr->world_x_dot->y();
      qOldDot[i+2] = ptr->world_x_dot->z();

      qOld[i] = ptr->world_x->x();
      qOld[i+1] = ptr->world_x->y();
      qOld[i+2] = ptr->world_x->z();

      //gravity
      f[i+1] = -g * ptr->rho * (ptr->material_coordinate - ptr->left->material_coordinate);

      //spring
      if(ptr != right_of_e0->left && ptr != right_of_e0) { //normal case

         deltaXDotRight = *(ptr->world_x_dot) - *(ptr->right->world_x_dot);
         deltaXDotLeft = *(ptr->world_x_dot) - *(ptr->left->world_x_dot);
         deltaXRight = *(ptr->world_x) - *(ptr->right->world_x);
         deltaXLeft = *(ptr->world_x) - *(ptr->left->world_x);

         absDeltaXRight = sqrt(deltaXRight.dot(deltaXRight));
         absDeltaXLeft = sqrt(deltaXLeft.dot(deltaXLeft));

         constant = -(ks * (absDeltaXLeft - l0) + kd * deltaXDotLeft.dot(deltaXLeft)/absDeltaXLeft)*deltaXLeft/absDeltaXLeft - (ks * (absDeltaXRight - l0) + kd * deltaXDotRight.dot(deltaXRight)/absDeltaXRight)*deltaXRight/absDeltaXRight;
         
      }
      else if(ptr == right_of_e0->left) { //ball's left node is ptr

         deltaXDotRight = *(ptr->world_x_dot) - *(previous_state->e1->world_x_dot);
         deltaXDotLeft = *(ptr->world_x_dot) - *(ptr->left->world_x_dot);
         deltaXRight = *(ptr->world_x) - *(previous_state->e1->world_x);
         deltaXLeft = *(ptr->world_x) - *(ptr->left->world_x);
         restLengthRight = l0 * alpha;
         restLengthRight = l0;
         
         if(deltaXRight.dot(deltaXRight) < 0.001) {
            deltaXDotRight = *(ptr->world_x_dot) - *(ptr->right->world_x_dot);
            deltaXRight = *(ptr->world_x) - *(ptr->right->world_x);
            restLengthLeft = l0;
         }
         
         absDeltaXRight = sqrt(deltaXRight.dot(deltaXRight));
         absDeltaXLeft = sqrt(deltaXLeft.dot(deltaXLeft));          
         constant = -(ks * (absDeltaXLeft - restLengthLeft) + kd * deltaXDotLeft.dot(deltaXLeft)/absDeltaXLeft)*deltaXLeft/absDeltaXLeft - (ks * (absDeltaXRight - restLengthRight) + kd * deltaXDotRight.dot(deltaXRight)/absDeltaXRight)*deltaXRight/absDeltaXRight;
      }
      else if(ptr == right_of_e0) { //ball's right node is ptr
         deltaXDotLeft = *(ptr->world_x_dot) - *(previous_state->e1->world_x_dot);
         deltaXDotRight = *(ptr->world_x_dot) - *(ptr->right->world_x_dot);
         deltaXLeft = *(ptr->world_x) - *(previous_state->e1->world_x);
         deltaXRight = *(ptr->world_x) - *(ptr->right->world_x);
         restLengthLeft = l0 * alpha;
         restLengthRight = l0;
         
         if(deltaXLeft.dot(deltaXLeft) < 0.001) {
            deltaXDotLeft = *(ptr->world_x_dot) - *(ptr->left->world_x_dot);
            deltaXLeft = *(ptr->world_x) - *(ptr->left->world_x);
            restLengthLeft = l0;
         }
         
         absDeltaXLeft = sqrt(deltaXLeft.dot(deltaXLeft));
         absDeltaXRight = sqrt(deltaXRight.dot(deltaXRight));
         constant = -(ks * (absDeltaXLeft - restLengthLeft) + kd * deltaXDotLeft.dot(deltaXLeft)/absDeltaXLeft)*deltaXLeft/absDeltaXLeft - (ks * (absDeltaXRight - restLengthRight) + kd * deltaXDotRight.dot(deltaXRight)/absDeltaXRight)*deltaXRight/absDeltaXRight;
      }
      
      f[i] += constant.x();
      f[i+1] += constant.y();
      f[i+2] += constant.z();

      f[i] = 0.0; f[i+2] = 0.0;
      ptr = ptr->right;
   }


   Vector3f x0_and_x1= (1-alpha)* *(right_of_e0->left->world_x) + alpha* *(right_of_e0->world_x);
    f[3*(NUM-1)] = previous_state->e1->rho * x0_and_x1.dot(Vector3f(0,-g,0)); 

   VectorXf rhs = M * qOldDot +  dt * f;
   VectorXf qNewDot;

   qNewDot = M.lu().solve(rhs);


   VectorXf qNew = qNewDot * dt + qOld;

   L3node * newRoot = new L3node(0.0, *(previous_state->root->world_x), *(previous_state->root->world_x_dot));
   state->root = newRoot;
   node * newPtr = state->root;
   ptr = previous_state->root->right;

   Vector3f newX, newXDot;

   //update the new state
   for(int i=0; i<3*(NUM-1); i+=3) {

      newX = Vector3f(qNew[i], qNew[i+1], qNew[i+2]);
      newXDot = Vector3f(qNewDot[i], qNewDot[i+1], qNewDot[i+2]);
      
      L3node * newNode = new L3node(ptr->material_coordinate, newX, newXDot);
      
      newPtr->right = newNode;
      newPtr->right->left = newPtr;
      
      ptr = ptr->right;
      newPtr = newPtr->right;
   }

   L3node* rightNode = new L3node(1.0, *(ptr->world_x), *(ptr->world_x_dot));
   newPtr->right = rightNode;
   newPtr->right->left = newPtr;

   cout << "qNew : " << qNew << "qDotNew: " << qNewDot << endl;
   right_of_e0 = linear_search(qNew[3*(NUM-1)], newRoot);
   alpha = (qNew[3*(NUM-1)] - right_of_e0->left->material_coordinate) / (right_of_e0->material_coordinate - right_of_e0->left->material_coordinate);
   newX = (1.0-alpha) * *(right_of_e0->left->world_x) + alpha * *(right_of_e0->world_x);
   newXDot = (1.0-alpha) * *(right_of_e0->left->world_x_dot) + alpha * *(right_of_e0->world_x_dot);

   E0node * newE = new E0node(qNew[3*(NUM-1)], qNewDot[3*(NUM-1)], newX, newXDot, newRoot);
   newE->rho = 2.0;
   state->e1 = newE;
   state->next = NULL;
   return state;

 /*
// state->e1 = dynamic_cast<E3node *> (previous_state->e1)->update();
 //state->e1 = dynamic_cast<L3node *> (previous_state->e1)->update();
 state->e1 = dynamic_cast<E0node *> (previous_state->e1)->update();
 state->root = dynamic_cast<E0node *> (state->e1)->rootNode;
 state->next = NULL;
 return state;
*/

}

void buildFrames(){ //builds frames into a cyclic finite state machine

    State * prevState = NULL;
//    if(prevState != NULL && abs(present_state->e1->material_coordinate - prevState->e1->material_coordinate) < 0.001 ) {
      if(iter > 50) {
 //      if(present_state->e1->material_coordinate > 0.99){
        present_state->next = initial_state;
       
        glutTimerFunc(1, advanceT, 0);
        return;
    }
    else {
      iter++;
        prevState = present_state;
        present_state->next = takeastep(present_state);
        
        assert(present_state->next != NULL); //in case of an error, make program fail
       
        present_state = present_state->next;
        
        buildFrames();
    }
}

void saveToFile(int value) { //outputs a bmp file
    
    glReadBuffer(GL_FRONT);
        
    unsigned char *pRGB = new unsigned char [ 3 * viewport.w * viewport.h + 3];
    glReadBuffer(GL_FRONT_AND_BACK);
    glReadPixels(0, 0, viewport.w, viewport.h, GL_RGB, GL_UNSIGNED_BYTE, pRGB);
        
    FreeImage_Initialise();
        
    FIBITMAP * bitmap = FreeImage_ConvertFromRawBits((BYTE *) pRGB, viewport.w, viewport.h, 3*viewport.w, 24, 0, 0, 0, true);
    
    bitmap = FreeImage_Rotate(bitmap, 180);
    
    stringstream ss;
    ss << "frame" << value << ".bmp";
    string mystring = ss.str();
    const char * pchar = mystring.c_str();
    FreeImage_Save(FIF_BMP, bitmap, pchar, 0);
        
    FreeImage_DeInitialise();	

    if (present_state->next != initial_state) {
        glutTimerFunc(1000/30, saveToFile, value+1);
    }
    else {
        printf("done!\n");
        printf("%d frames total\n", value);
    }
}

void myReshape(int w, int h) { // standard OpenGL idle function
	glViewport(0,0, (GLsizei) w, (GLsizei) h);// sets the rectangle that will be the window
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();				// loading the identity matrix for the screen
	gluPerspective(65.0, (GLfloat) w/(GLfloat) h, 1.0, 20.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0, 0.0, -5.0);
}

State * state_initialize() { //for now, hard coding in initial parameters
    
    State * state = new State;
    
    double s = 1.0/NUM;
    Vector3f pos0(-2.0, 0.0, 0.0);
    Vector3f pos1(2.0, 0.0, 0.0);
    Vector3f pos(0,0,0);
    Vector3f vel(0,0,0);
    L3node * first_node = new L3node(0.0, pos0, vel); //for now, NULL == RIGID, i.e. do nothing
    L3node * last_node = new L3node(1.0, pos1, vel);
    L3node * nth_node;
    node * prev_node = first_node;
    first_node->cons = RIGID;

    // automatically generates n nodes betwen the stationary nodes L0 and L1
   for(int i=1; i<NUM; i++) {
        
        pos = Vector3f(s*pos1 + (1.0-s)*pos0);
        nth_node = new L3node(s, pos, vel);
        nth_node->cons = NORMAL;
        s += 1.0/NUM;
        nth_node->left = prev_node;
      //  nth_node->cons_node = NULL;
        prev_node->right = nth_node;
        
        prev_node = nth_node;
    }
    nth_node->right = last_node;
 
//    first_node->right = last_node;
    last_node->left = nth_node;
   // last_node->left = first_node;
    last_node->cons = RIGID;

    /********** adding an E3_node ****************/
    s = 0.5/NUM;
    pos = Vector3f((s)*pos1 + (1.0-s)*pos0);
//    L3node * interesting_node = new L3node(s, pos, vel);
//    E3node* interesting_node = new E3node(s, 0.0, pos, vel, first_node);
    E0node * interesting_node = new E0node(s, 0.0, pos, vel, first_node);
    interesting_node->rho = 2.0;
    interesting_node->rootNode = first_node;
    /********************************************/
    
    state->root = first_node;
    state->next = NULL;
    state->e1 = interesting_node;
    return state;
}

void display_lines() { //connects all the lines using all nodes
    
    glBegin(GL_LINES);
    
    node* ptr = present_state->root;
    node* e1 = present_state->e1;
    node* left_of_e1 = (linear_search(e1->material_coordinate, ptr))->left;
    glVertex3f(ptr->world_x->operator[](0),ptr->world_x->operator[](1),ptr->world_x->operator[](2));
    
    while(ptr->right->right != NULL) {

       if(ptr == left_of_e1) {
         glVertex3f(ptr->world_x->operator[](0),ptr->world_x->operator[](1),ptr->world_x->operator[](2));
         glVertex3f(ptr->world_x->operator[](0),ptr->world_x->operator[](1),ptr->world_x->operator[](2));
         glVertex3f(e1->world_x->operator[](0),e1->world_x->operator[](1),e1->world_x->operator[](2));
         glVertex3f(e1->world_x->operator[](0),e1->world_x->operator[](1),e1->world_x->operator[](2));
       } 
         ptr = ptr->right;
        
        glVertex3f(ptr->world_x->operator[](0),ptr->world_x->operator[](1),ptr->world_x->operator[](2));
        glVertex3f(ptr->world_x->operator[](0),ptr->world_x->operator[](1),ptr->world_x->operator[](2));
    }
    ptr = ptr->right;
    glVertex3f(ptr->world_x->operator[](0),ptr->world_x->operator[](1),ptr->world_x->operator[](2));
	glEnd();
    
}

void display_nodes() { //display all massless nodes as spheres
    
    node * ptr = present_state->root;
    node* e1 = present_state->e1;
    node* left_of_e1 = (linear_search(e1->material_coordinate, ptr))->left;

    glTranslatef(ptr->world_x->operator[](0), ptr->world_x->operator[](1), ptr->world_x->operator[](2));
    glutSolidSphere(radius, 20, 20);
    while(ptr->right->right != NULL) {

       if(ptr == left_of_e1) {
         
         glTranslatef(e1->world_x->operator[](0)-ptr->world_x->operator[](0), e1->world_x->operator[](1)-ptr->world_x->operator[](1), e1->world_x->operator[](2)-ptr->world_x->operator[](2));

         glutSolidSphere(radius * 2, 20, 20); 
         ptr = ptr->right;

         glTranslatef(ptr->world_x->operator[](0)-e1->world_x->operator[](0), ptr->world_x->operator[](1)-e1->world_x->operator[](1), -ptr->world_x->operator[](2)-e1->world_x->operator[](2));
         glutSolidSphere(radius, 20, 20);
       } 
      else {
        
            ptr = ptr->right;
            glTranslatef(ptr->world_x->operator[](0)-ptr->left->world_x->operator[](0), ptr->world_x->operator[](1)-ptr->left->world_x->operator[](1), ptr->world_x->operator[](2)-ptr->left->world_x->operator[](2));
            glutSolidSphere(radius, 20, 20);
         } 
   }
    ptr = ptr->right;
    glTranslatef(ptr->world_x->operator[](0)-ptr->left->world_x->operator[](0), ptr->world_x->operator[](1)-ptr->left->world_x->operator[](1), ptr->world_x->operator[](2)-ptr->left->world_x->operator[](2));
    glutSolidSphere(radius, 20, 20);
}
   
void initScene(){ //standard OpenGL function to start scene, calls state_initialize
    
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f); // Clear to black, fully transparent
	glClearDepth(0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glTranslatef(0.0, 0.0, -5.0);
	glPushMatrix();
	glColor3f(0.0,0.0,1.0);
	
    /****************initializes the first frame and calls buildFrames to calculate all Frames **/
    initial_state = state_initialize();
    present_state = initial_state;

    Vector4f q0(-2.0,0,0,0);
    Vector4f q1(-1.6,0,0,0.1);

    buildFrames();
    /*****************************************************************************************/
    
	glPushMatrix();    
    /***************** these functions make this code input-independent **********************/
    display_lines();
    display_nodes();
    /*****************************************************************************************/
	glPopMatrix();

	glShadeModel (GL_FLAT);
	glutSwapBuffers();
    
    glutTimerFunc(1000/30, saveToFile, 0);
}

void initLights() { // std OpenGL func to put in lights
	glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	GLfloat global_ambient[] = {.1f, .1f, .1f};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);
	
	GLfloat ambient[] = {0.1f, 0.1f, 0.1f};
	GLfloat diffuse[] = {1.0f, 1.0f, 1.0f};
	GLfloat specular[] = {1.0f, 1.0f, 1.0f};
	GLfloat pos[] = {0, 0, 1, 0};
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT0, GL_POSITION, pos);
	glEnable(GL_LIGHT0);
	
	GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_diffuse[] = {1, 0, 1, 1.0};
	GLfloat mat_ambient[] = {.1, .1, .1, 1.0};
	GLfloat mat_shininess[] = {50.0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
    
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
}

void mouseTracker(int button, int state, int x, int y){
	if(button == 0 && state == GLUT_DOWN){
		win_x = (x/viewport.w - 0.5)* 10.0;
		win_y = (viewport.h/2.0 - y) / viewport.h * 10.0;
		win_z = 0.0;
		printf("%d %d, %f %f\n", x, y, win_x, win_y);
		glutPostRedisplay();
	}
}

void myDisplay() { //std OpenGL display func, calls display() method of all nodes	
    glClear(GL_COLOR_BUFFER_BIT);				// clear the color buffer (sets everything to black)
	glViewport(0,0, (GLsizei) viewport.w, (GLsizei) viewport.h);// sets the rectangle that will be the window
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();				// loading the identity matrix for the screen
	gluPerspective(65.0, (GLfloat) viewport.w/(GLfloat) viewport.h, 1.0, 20.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    
	glColor3f(0,0,1);
	glTranslatef(0.0, 0.0, -5.0);

	glPushMatrix();   
    /***************** these functions make this code input-independent **********************/
    //present_state = present_state->next;
    display_lines();
    display_nodes();
    /*****************************************************************************************/	
	glPopMatrix();
    
	glShadeModel(GL_FLAT);
	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set 
}

void processNormalKeys(unsigned char key, int x, int y) {
    
	if (key == 32) //the space key
		glutPostRedisplay();
    
	if (key == 27) {
        // when receiving the escape key the function escapes
		exit(0);
	}
}

void myFrameMove() {
    
	glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}

void myGlutIdle( void) {
    glutSetWindow(main_window);
    glutPostRedisplay();
}

int main(int argc, char *argv[]) {	
    
  	glutInit(&argc, argv);
  	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);

  	viewport.w = 500;
  	viewport.h = 500;
   
    /***** this opens a dialog box to input the right nodes ****/
    
//    GLUI *glui = GLUI_Master.create_glui( "GLUI", 0 );
//    glui->add_statictext( "Simple GLUI Example" );
//    glui->add_separator();
//    //glui->add_checkbox( "Wireframe", &wireframe, 1, control_cb );
//    //GLUI_Spinner *segment_spinner =
//    //glui->add_spinner( "Segments:",GLUI_SPINNER_INT, &segments );
//    //segment_spinner->set_int_limits( 3, 60, GLUI_LIMIT_WRAP );
//    //GLUI_EditText *edittext =
//    //glui->add_edittext( "Text:", GLUI_EDITTEXT_TEXT, text );
//    glui->add_column(true);       /**  Begin new column - 'true' indicates   **/
//    /* *  a vertical bar should be drawn        **/
//    GLUI_Panel *obj_panel = glui->add_panel( "Object Type" );
//    //GLUI_RadioGroup *group1 =
//   // glui->add_radiogroup_to_panel(obj_panel,&obj,3,control_cb);
//    //glui->add_radiobutton_to_group( group1, "Sphere" );
//    //glui->add_radiobutton_to_group( group1, "Torus" );
//    glui->add_button( "Quit", 0,(GLUI_Update_CB)exit );
//    /** Tell GLUI window which other window to recognize as the main gfx window **/
//    glui->set_main_gfx_window( main_window );
//    /** Register the Idle callback with GLUI (instead of with GLUT)     **/
//    GLUI_Master.set_glutIdleFunc( myGlutIdle );    
    
    /***********************************************************/

    
    glutInitWindowSize(viewport.w, viewport.h);
  	glutInitWindowPosition(0, 0);
  	main_window = glutCreateWindow("Toy Model");
    
  	initScene();							// quick function to set up scene
	glutDisplayFunc(myDisplay);
    
	glutKeyboardFunc(processNormalKeys);
    //	glutSpecialFunc(processSpecialKeys);
  	glutReshapeFunc(myReshape);				// function to run when the window gets resized
  	initLights();
	glEnable(GL_DEPTH_TEST);
	glutMouseFunc(mouseTracker);
	glutIdleFunc(myFrameMove);				// function to run when not handling any other task
    
    glutMainLoop();							// infinite loop that will keep drawing and resizing and whatever else
    
  	return 0;
}

