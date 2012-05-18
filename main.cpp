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
int main_window;
extern int NUM; //number of Lagrangian nodes in between the 2 rigid nodes


int iter = 0;


/***********************************************************************/

extern L1node* pendulum_update(L1node * old_node);
extern L3node* string_update(L3node * old_node);
extern E3node* mass_node(E3node * old_node);

void advanceT(int value) { //jump to the next state every 5 ms
    
    present_state = present_state->next;
    glutTimerFunc(10, advanceT, 0);
}

State * takeastep(State * previous_state) { // calculates the next state in the simulation
                                            // this is where the node::update funcs are called    
 State * state = new State;
// state->e1 = dynamic_cast<E3node *> (previous_state->e1)->update();
 //state->e1 = dynamic_cast<L3node *> (previous_state->e1)->update();
 state->e1 = dynamic_cast<E0node *> (previous_state->e1)->update();
 state->root = dynamic_cast<E0node *> (state->e1)->rootNode;
 state->next = NULL;
 return state;

//    State * state = new State;
//    L3node * ptr = present_state->root;
//    E1node * interesting_node = dynamic_cast<E1node*> (present_state->e1);
//    L3node * right_of_e1, prev_node;
//    double s0, s1;
//    
//    while(ptr->right != NULL) { //linearly traverse and finds [s0,s1] for E1node, adds mg to 
//                                // these 2 nodes
//        if(interesting_node->material_coordinate >= ptr->material_coordinate && 
//           interesting_node->material_coordinate <= ptr->right->material_coordinate) {
//            
//            s0 = ptr->material_coordinate;
//            s1 = ptr->right->material_coordinate;
//            
//            right_of_e1 = ptr->right;
//                //        ptr = ptr->update(-mg);
//                //    right_of_e1 = right_of_e1->update(-mg);
//            ptr->right = right_of_e1;
//            right_of_e1->left = ptr;
//            
//            while(ptr->left != NULL) {
//                prev_node = ptr;
//                ptr= ptr->left;
//                    //    ptr = ptr->update(0);
//                prev_node->left = ptr;
//                ptr->right = prev_node;
//            }
//            
//            while(right_of_e1->right != NULL) {
//                prev_node = right_of_e1;
//                right_of_e1 = right_of_e1->right;
//                    //               right_of_e1 = right_of_e1->udpate(0);
//                prev_node->right = right_of_e1;
//                right_of_e1->left = prev_node;
//            }
//
//                //    interesting_node= interesting_node->update(ptr, ptr->right);
//            
//            right_of_e1->right = interesting_node;
//            interesting_node->left = right_of_e1;
//            
//            break; //found interval, break out of while loop
//        }
//        else
//            ptr = ptr->right;
//    }
//    
//    if(ptr->left != NULL) {
//        printf("error! interval was never found\n");
//        return NULL;
//    }
//    
//    state->root = ptr;
//   
//    state->next = NULL;
//    return state;

}

void buildFrames(){ //builds frames into a cyclic finite state machine

    State * prevState = NULL;
//    if(prevState != NULL && abs(present_state->e1->material_coordinate - prevState->e1->material_coordinate) < 0.001 ) {
//      if(iter > 10000) {
       if(present_state->e1->material_coordinate > 0.99){
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
    node* left_of_e1 = (ptr->linear_search(e1->material_coordinate, ptr))->left;
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
    node* left_of_e1 = (ptr->linear_search(e1->material_coordinate, ptr))->left;

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

