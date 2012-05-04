#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "algebra3.h"
#include "FreeImage.h"
#include "Energy_Func.h"

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
    display_node * d_root;
    State * next;
};

Viewport viewport;
State* initial_state;
State* present_state;

/* some config variables for frames/sec, size of node, init rot angle, etc.*/
float radius = 0.1;
float win_x, win_y, win_z;
float dt = 1.0/30.0;
int main_window;
/***********************************************************************/

extern L1node* pendulum_update(L1node * old_node);

void advanceT(int value) { //jump to the next state every 5 ms
    
    present_state = present_state->next;
    glutTimerFunc(5, advanceT, 0);
}

State * takeastep(State * previous_state) { // calculates the next state in the simulation
    
    State * state = new State;
    state->root = dynamic_cast<L0node*>(previous_state->root)->update();
        //  state->d_root = dynamic_cast<bead*>(previous_state->d_root)->update(*(state->root->right->world_x));
   
    state->next = NULL;
    return state;

}

void buildFrames(){ //builds frames into a cyclic finite state machine

    L1node * myNode = dynamic_cast<L1node*> (present_state->root->right);
    L1node * prevNode = dynamic_cast<L1node*> (initial_state->root->right);
    
    if(abs((float)myNode->x->d - prevNode->x->d) < 0.01 && initial_state->next != NULL && abs((float) myNode->x_dot->d - prevNode->x_dot->d) < 0.01) {
        
        present_state->next = initial_state;
       
        glutTimerFunc(5 , advanceT, 0);
        return;
    }
    else {
    
        present_state->next = takeastep(present_state);
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
    
    vec3 pos0 = vec3(0.0, 0.0, 0.0);
    
    L0node * first_node = new L0node(0.0, pos0); //for now, NULL == RIGID, i.e. do nothing
    
    vec3 pos1 = vec3(l*sin(-PI/3.0), l*cos(-PI/3.0), 0.0);
    
    L1node * second_node = new L1node(1.0, 30.0 * PI/180.0, 0.0, pos1);
    second_node->left = first_node;
    second_node->cons_node = new circle_constraint(l, 0.0, second_node->world_x->operator[](2), second_node);
    
    first_node->right = second_node;
    first_node->func = NULL;
    second_node->func = &pendulum_update; 
    
    state->root = first_node;
//    state->d_root = new bead(radius, NULL, *(first_node->world_x));
//    state->d_root->next = new bead(radius, second_node->cons_node, *(second_node->world_x));
//    
    state->next = NULL;
    
    return state;
}
    
void initScene(){ //standard OpenGL function to start scene, calls state_initialize
    
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f); // Clear to black, fully transparent
	glClearDepth(0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glTranslatef(0.0, 0.0, -5.0);
	glPushMatrix();
	
    initial_state = state_initialize();
    present_state = initial_state;
    buildFrames();
    
	glColor3f(0.0,0.0,1.0);
	
	glutSolidSphere(radius,20,20);
    
	glPushMatrix();
	glBegin(GL_LINES);
	glVertex3f(0,0,0);
	glVertex3f(0,-2,0);
	glEnd();
	
	glTranslatef(0,-2,0);
	glutSolidSphere(radius,20,20);
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
    glRotatef(present_state->root->right->x->d * 180.0/PI, 0, 0, 1);
	glPushMatrix();
    
	glutSolidSphere(radius,20,20);
    
	glBegin(GL_LINES);
	glVertex3f(0,0,0);
	glVertex3f(0,-2,0);
	glEnd();
	
	glTranslatef(0,-2,0);
	glutSolidSphere(radius,20,20);
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

