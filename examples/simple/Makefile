INC_DIR= ../../include
LIB_DIR= ../../lib
BIN_DIR= .

LDFLAG=-pthread -lgstreamer-0.10 -lgobject-2.0 -lgmodule-2.0 -lgthread-2.0 -lxml2 -lglib-2.0 -L/usr/X11R6/lib -L$(LIB_DIR)
LIBS= -lARgsub -lARvideo -lAR -lpthread -lglut -lGLU -lGL -lXi -lX11 -lm -pthread -lgstreamer-0.10 -lgobject-2.0 -lgmodule-2.0 -lgthread-2.0 -lxml2 -lglib-2.0

#with optimization
#CFLAG= -O -pthread -I/usr/include/gstreamer-0.10 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/libxml2 -I/usr/X11R6/include -g -I$	(INC_DIR)

#without optimization
CFLAG= -O0 -pthread -I/usr/include/gstreamer-0.10 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/libxml2 -I/usr/X11R6/include -g -I $(INC_DIR)

INC_DIR_OTHER= -lGLEW -lboost_thread -lboost_system
INC_SOME= /usr/local/include
TOTAL=  $(INC_DIR_OTHER) -pthread `pkg-config --cflags opencv``pkg-config --libs opencv` -L/usr/local/lib
# options to pass to the compiler
CXXFlags3= -Wall -g -I. -I$(INC_SOME)


all: $(BIN_DIR)/simpleTest

$(BIN_DIR)/simpleTest: clean random1.o simpleTest.o $(OBJS)
	g++ -o $(BIN_DIR)/simpleTest simpleTest.o random1.o maya.o $(TOTAL) $(OBJS) $(LDFLAG) $(LIBS)

simpleTest.o: simpleTest.c $(HEADDERS)
	gcc -c $(CFLAG) simpleTest.c

random1.o: random1.cpp maya.cpp
	g++ -c  $(CXXFlags3) random1.cpp maya.cpp  
cross.o:
	g++ -c cross.cpp

clean:
	rm -f *.o
	rm -f $(BIN_DIR)/simpleTest

allclean:
#	rm -f *.o
#	rm -f $(BIN_DIR)/simpleTest
#	rm -f Makefile
