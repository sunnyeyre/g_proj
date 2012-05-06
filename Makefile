#!gmake 

#-----------------------------------------
#Basic Stuff -----------------------------
BITTINESS = -m32
CC          = g++ -Wall ${BITTINESS} -ggdb
cc          = gcc -Wall ${BITTINESS}

#-----------------------------------------
#Optimization ----------------------------
OPTMAC = -fast

SOURCES = main.cpp
OBJECTS = Energy_Func.o main.o
TARGETS = test

#-----------------------------------------
# Mac specific stuff
FRAMEWORK = -framework GLUT
FRAMEWORK += -framework OpenGL
MACLIBS = -lglut -lGLU -lGL -lXmu -lX11 -lm -lstdc++ -lfreeimage -lGLUI
MACINCS = -L"/System/Library/Frameworks/OpenGL.framework/Libraries" -I"/System/Library/Frameworks/GLUT.framework/Headers" -L ./ -I ./ -L /opt/local/lib/ -I /opt/local/include/ -L /usr/local/lib/ -I /usr/local/include/ -I/Developer-3.0/SDKs/MacOSX10.6.sdk/usr/X11/include -L/Developer-3.0/SDKs/MacOSX10.6.sdk/usr/X11/lib/ -L eigen/ -I eigen/ -I eigen/Eigen/

#-----------------------------------------
CCOPTSMAC = $(OPTMAC) $(MACINCS) -DOSX
LDOPTSMAC = $(OPTMAC) $(MACINCS) $(MACLIBS) -DOSX

#-----------------------------------------
#-----------------------------------------

default: $(TARGETS)

clean:
	/bin/rm -f *.o $(TARGETS)

#-----------------------------------------
#-----------------------------------------

test: $(OBJECTS)
	$(CC) $(LDOPTS) -L./ $(FRAMEWORK) $(LDOPTSMAC) $(OBJECTS) -o $(TARGETS)

main.o: main.cpp
	$(CC) main.cpp -c $(CCOPTSMAC) -o main.o
Energy_Func.o: Energy_Func.cpp
	$(CC) Energy_Func.cpp -c $(CCOPTSMAC) -I./ -I eigen/ -o Energy_Func.o

