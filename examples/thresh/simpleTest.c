#ifdef _WIN32
#include <windows.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glut.h>
#else
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif
#include <AR/gsub.h>
#include <AR/video.h>
#include <AR/param.h>
#include <AR/ar.h>
#include <string.h>
#include <time.h>
#include <math.h>
static int w,h;
float dbx[4], dby[4];
double dbx1[4], dby1[4];
static double testvar1 = 0.0;
int a;
char **b;
static struct timespec curr, last;
static double xm[2];
static double ym[2];


//
// Camera configuration.
//
#define Rishabh
#ifdef _WIN32
char			*vconf = "Data\\WDM_camera_flipV.xml";
#else
char			*vconf = "";
#endif

int             xsize, ysize;
int             thresh = 100;
int             count = 0;

char           *cparam_name    = "Data/camera_para.dat";
ARParam         cparam;

char           *patt_name      = "Data/patt.hiro";
int             patt_id;
double          patt_width     = 80.0;
double          patt_center[2] = {0.0, 0.0};
double          patt_trans[3][4];

static void   init(void);
static void   cleanup(void);
static void   keyEvent( unsigned char key, int x, int y);
static void   mainLoop(void);
static void   draw( void );
static void   draw1( void );
static void   draw2( void );
static void   draw3( void );

int main(int argc, char **argv)
{
	clock_gettime(CLOCK_MONOTONIC, &last);
	setbuf(stdout,NULL);
	a = argc;
	b = argv;
	glutInit(&argc, argv);
	init();

	arVideoCapStart();
	argMainLoop( NULL, keyEvent, mainLoop );
	return (0);
}

static void   keyEvent( unsigned char key, int x, int y)
{ 	
	static char stringtemp[] = "";
	static char command;
	static int val;
	//stringtemp = "";
	/* quit if the ESC key is pressed */
	if( key == 0x1b ) {
		printf("*** %f (frame/sec)\n", (double)count/arUtilTimer());
		cleanup();
		exit(0);
	}
	else if(key == '|'){
		cleanup();
		glutInit(&a,b);
		init();
	       	arVideoCapStart();
       		argMainLoop( NULL, keyEvent, mainLoop );
		
	}	
	else if(key == (char)13){
	/*	if(testvar1 == 1)
			testvar1 = 0;
		else
			testvar1 = 1;
		printf("Pressed |: test = %f \n",(double)testvar1); */
		printf("\n%s \n",stringtemp);
		sscanf(stringtemp,"%c %d",&command,&val);
		printf("%d",val);
		switch(command)
		{
			case 'w':
			w = val;
			break;
			case 'h':
			h = val;
			break;
			case 'd':
			printf("w=%d, h=%d",w,h);
		}	
		strcpy(stringtemp, "");
		printf("\n");
	}
	else{
	printf("%c",key);
	strncat(stringtemp, &key, 1);
	}

}

/* main loop */
static void mainLoop(void)
{
	ARUint8         *dataPtr;
	ARUint8		*threshptr;
	ARUint8		*clipped_image;
	ARMarkerInfo    *marker_info;
	ARMarkerInfo	*marker_info_1;
	int             marker_num;
	int		marker_num_1;
	int             j, k;
	int		a,b;
	long i;
	int 		w_clipped,h_clipped;
//	for(i = 0;i<100000;i++);
	/* grab a vide frame */
	if( (dataPtr = (ARUint8 *)arVideoGetImage()) == NULL ) {
		arUtilSleep(2);
		return;
	}

	if( count == 0 ) arUtilTimerReset();
	count++;
	
	argDrawMode2D();
	threshptr = malloc(xsize*ysize*3);
	arCreateThresh(dataPtr, threshptr,&thresh,xsize,ysize,25);
	printf("%d\n",thresh);
	/* detect the markers in the video frame */
	if( arDetectMarker(dataPtr, thresh, &marker_info, &marker_num) < 0 ) {
		cleanup();
		exit(0);
	}
	argDispImage( threshptr, 0, 0);
	arVideoCapNext();

	/* check for object visibility */
	k = -1;
	for( j = 0; j < marker_num; j++ ) {
		if( patt_id == marker_info[j].id ) {
			if( k == -1 ) k = j;
			else if( marker_info[k].cf < marker_info[j].cf ) k = j;
		}
	}
	if( k == -1 ) {
		argSwapBuffers();
		return;
	}

	
	dbx[0] = marker_info[k].vertex[0][0];
	dby[0] = marker_info[k].vertex[0][1];

	dbx[1] = marker_info[k].vertex[1][0];
	dby[1] = marker_info[k].vertex[1][1];

	dbx[2] = marker_info[k].vertex[2][0];
	dby[2] = marker_info[k].vertex[2][1];

	dbx[3] = marker_info[k].vertex[3][0];
	dby[3] = marker_info[k].vertex[3][1];






	minmax(dbx,dby,xm,ym);
	printf("%f %f, %f %f",xm[0],xm[1],ym[0],ym[1]);
	clipped_image = extract(dataPtr,xm,ym,xsize,&w_clipped,&h_clipped);	
	printf("\t w=%d h=%d \t",w_clipped,h_clipped);
	argDispImageDrawPixels1( clipped_image , w_clipped, h_clipped);

	arDetect_inverted(clipped_image, thresh,dbx1,dby1 );

	dbx1[0]+=xm[0];
	dbx1[1]+=xm[0];
	dbx1[2]+=xm[0];
	dbx1[3]+=xm[0];
 
 	dby1[0]+=ym[0];
	dby1[1]+=ym[0];
	dby1[2]+=ym[0];
	dby1[3]+=ym[0];

/*
	marker_info[k].line[0][0] = 0;	
	marker_info[k].line[0][1] = 0;	
	marker_info[k].line[0][2] = 0;	
	marker_info[k].line[1][0] = 0;	
	marker_info[k].line[1][1] = 0;	
	marker_info[k].line[1][2] = 0;	
	marker_info[k].line[0][0] = 0;	
*/
	
	scale(dbx1,dby1);
	
	for(a = 0; a < 4; a++)
	{
		double dist = (marker_info[k].vertex[a][0]-dbx1[0])*(marker_info[k].vertex[a][0]-dbx1[0])+
					(marker_info[k].vertex[a][1]-dby1[0])*(marker_info[k].vertex[a][1]-dby1[0]);
		double dist2;
		int flag = 0;
		for(b = 1; b < 4; b++)
		{
			dist2=(marker_info[k].vertex[a][0]-dbx1[b])*(marker_info[k].vertex[a][0]-dbx1[b])+
				(marker_info[k].vertex[a][1]-dby1[b])*(marker_info[k].vertex[a][1]-dby1[b]);
			if(dist>dist2)
			{dist=dist2; flag = b;}
		}

		marker_info[k].vertex[a][0]+=dbx1[flag];
		marker_info[k].vertex[a][1]+=dby1[flag];       

		marker_info[k].vertex[a][0]/=2;
		marker_info[k].vertex[a][1]/=2;


	}

	for(a = 0; a < 4; a++)
	{
		double p,r;
		p = - (marker_info[k].vertex[a][1] - marker_info[k].vertex[(a+1)%4][1])/
				 (marker_info[k].vertex[a][0] - marker_info[k].vertex[(a+1)%4][0]);
		
		marker_info[k].line[a][0] = p/sqrt(1+p*p);

		marker_info[k].line[a][1] = 1/sqrt(1+p*p);

		marker_info[k].line[a][2] = -marker_info[k].line[a][0]*marker_info[k].vertex[a][0]- 
				marker_info[k].line[a][1]*marker_info[k].vertex[a][1];
	}


/*
	marker_info[k].vertex[0][0]+=dbx1[0];
	marker_info[k].vertex[0][1]+=dby1[0];

	marker_info[k].vertex[1][0]+=dbx1[1];
	marker_info[k].vertex[1][1]+=dby1[1];

	marker_info[k].vertex[2][0]+=dbx1[2];
	marker_info[k].vertex[2][1]+=dby1[2];
	
	marker_info[k].vertex[3][0]+=dbx1[3];
	marker_info[k].vertex[3][1]+=dby1[3];
		


	marker_info[k].vertex[0][0]/=2;
	marker_info[k].vertex[0][1]/=2;

	marker_info[k].vertex[1][0]/=2;
	marker_info[k].vertex[1][1]/=2;

	marker_info[k].vertex[2][0]/=2;
	marker_info[k].vertex[2][1]/=2;
	
	marker_info[k].vertex[3][0]/=2;
	marker_info[k].vertex[3][1]/=2;
*/		

	dbx[0] = marker_info[k].vertex[0][0];
	dby[0] = marker_info[k].vertex[0][1];

	dbx[1] = marker_info[k].vertex[1][0];
	dby[1] = marker_info[k].vertex[1][1];

	dbx[2] = marker_info[k].vertex[2][0];
	dby[2] = marker_info[k].vertex[2][1];

	dbx[3] = marker_info[k].vertex[3][0];
	dby[3] = marker_info[k].vertex[3][1];



	/* get the transformation between the marker and the real camera */
	arGetTransMat(&marker_info[k], patt_center, patt_width, patt_trans);
	

	//printf("daf %f \n",testvar1);
	draw1();
//	draw2();
//	draw3();	//reconstructed
	draw();
	argSwapBuffers();
//	free(threshptr);


}

static void init( void )
{
	ARParam  wparam;

	/* open the video path */
	if( arVideoOpen( vconf ) < 0 ) exit(0);
	/* find the size of the window */
	if( arVideoInqSize(&xsize, &ysize) < 0 ) exit(0);
	printf("Image size (x,y) = (%d,%d)\n", xsize, ysize);

	/* set the initial camera parameters */
	if( arParamLoad(cparam_name, 1, &wparam) < 0 ) {
		printf("Camera parameter load error !!\n");
		exit(0);
	}
	arParamChangeSize( &wparam, xsize, ysize, &cparam );
	arInitCparam( &cparam );
	printf("*** Camera Parameter ***\n");
	arParamDisp( &cparam );

	if( (patt_id=arLoadPatt(patt_name)) < 0 ) {
		printf("pattern load error !!\n");
		exit(0);
	}

	/* open the graphics window */
	argInit( &cparam, 1.0, 0, 0, 0, 0 );
}

/* cleanup function called when program exits */
static void cleanup(void)
{
	arVideoCapStop();
	arVideoClose();
	argCleanup();
}

static void draw1( void )
{
	/*    double    gl_para[16];
	      GLfloat   mat_ambient[]     = {0.0, 0.0, 1.0, 1.0};
	      GLfloat   mat_flash[]       = {0.0, 0.0, 1.0, 1.0};
	      GLfloat   mat_flash_shiny[] = {50.0};
	      GLfloat   light_position[]  = {100.0,-200.0,200.0,0.0};
	      GLfloat   ambi[]            = {0.1, 0.1, 0.1, 0.1};
	      GLfloat   lightZeroColor[]  = {0.9, 0.9, 0.9, 0.1};

	      argDrawMode3D();
	      argDraw3dCamera( 0, 0 );
	      glClearDepth( 1.0 );
	      glClear(GL_DEPTH_BUFFER_BIT);
	      glEnable(GL_DEPTH_TEST);
	      glDepthFunc(GL_LEQUAL);


	/* load the camera transformation matrix *
	argConvGlpara(patt_trans, gl_para);
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixd( gl_para );

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambi);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_flash);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_flash_shiny);	
	glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
	glMatrixMode(GL_MODELVIEW);

	//	glTranslatef( 0.0, 0.0, 25.0 );
	//	glutSolidCube(50.0);
	//	glLoadIdentity();
	//	glTranslatef(0.0,0.0,100.0);
	//	glutSolidCube(5.0);
	glLoadIdentity();
	glColor3f(1.0f,0.0f,0.0f);
	glPointSize(5.0f);
	glBegin(GL_POINTS);
	glVertex3f(debug1,debug2,50.0);
	glEnd();
	glDisable( GL_LIGHTING );
	glDisable( GL_DEPTH_TEST );
	 */    
	//const XSize = 640, YSize = 480;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho (0, xsize, ysize, 0, 0, 1);
	glMatrixMode (GL_MODELVIEW);
	glDisable(GL_DEPTH_TEST);
	glColor3f(1.0f,0.0f,0.0f);
	glPointSize(5.0f);

//	printf("points: %f, %f \t",dbx[0],dby[0]);
//	printf("%f, %f \t",dbx[1],dby[1]);
//	printf("%f, %f \t",dbx[2],dby[2]);
//	printf("%f, %f \t",dbx[3],dby[3]);

	if(testvar1 == 0.0){

		glBegin(GL_LINE_LOOP);
		glVertex2f(dbx[0],dby[0]);
		glVertex2f(dbx[1],dby[1]);
		glVertex2f(dbx[2],dby[2]);
		glVertex2f(dbx[3],dby[3]);
	}

	else if(testvar1 == 1.0){
		glBegin(GL_LINES);
		glVertex2f(0,dby[0]);
		glVertex2f(dbx[0],0);

		glVertex2f(0,dby[1]);
		glVertex2f(dbx[1],0);

		glVertex2f(0,dby[2]);
		glVertex2f(dbx[2],0);

		glVertex2f(0,dby[3]);
		glVertex2f(dbx[3],0);

		glVertex2f(0,dby[4]);
		glVertex2f(dbx[4],0);
	}

	glEnd();


}
static void draw( void )
{
    double    gl_para[16];
    GLfloat   mat_ambient[]     = {0.0, 0.0, 1.0, 1.0};
    GLfloat   mat_flash[]       = {0.0, 0.0, 1.0, 1.0};
    GLfloat   mat_flash_shiny[] = {50.0};
    GLfloat   light_position[]  = {100.0,-200.0,200.0,0.0};
    GLfloat   ambi[]            = {0.1, 0.1, 0.1, 0.1};
    GLfloat   lightZeroColor[]  = {0.9, 0.9, 0.9, 0.1};

    argDrawMode3D();
    argDraw3dCamera( 0, 0 );
    glClearDepth( 1.0 );
    glClear(GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    /* load the camera transformation matrix */
    argConvGlpara(patt_trans, gl_para);
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd( gl_para );

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambi);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_flash);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_flash_shiny);
    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMatrixMode(GL_MODELVIEW);
    glTranslatef( 0.0, 0.0, 25.0 );
    glutSolidCube(50.0);
    glDisable( GL_LIGHTING );

    glDisable( GL_DEPTH_TEST );
}
static void draw3( void )
{
    double    gl_para[16];
    GLfloat   mat_ambient[]     = {0.0, 0.0, 1.0, 1.0};
    GLfloat   mat_flash[]       = {0.0, 0.0, 1.0, 1.0};
    GLfloat   mat_flash_shiny[] = {50.0};
    GLfloat   light_position[]  = {100.0,-200.0,200.0,0.0};
    GLfloat   ambi[]            = {0.1, 0.1, 0.1, 0.1};
    GLfloat   lightZeroColor[]  = {0.9, 0.9, 0.9, 0.1};

    argDrawMode3D();
    argDraw3dCamera( 0, 0 );
    glClearDepth( 1.0 );
    glClear(GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);

    /* load the camera transformation matrix */
    argConvGlpara(patt_trans, gl_para);
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd( gl_para );

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambi);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightZeroColor);
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_flash);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_flash_shiny);
    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
    glMatrixMode(GL_MODELVIEW);
//    glTranslatef( 0.0, 0.0, 25.0 );
// z = -50.0f;
    glBegin(GL_LINE_LOOP);

float z = 0;
    	glVertex3f(40.0f, 40.0f, z);
    	glVertex3f(40.0f, -40.0f, z);
    	glVertex3f(-40.0f, -40.0f, z);
    	glVertex3f(-40.0f, 40.0f, z);

	glEnd();
    glDisable( GL_LIGHTING );

    glDisable( GL_DEPTH_TEST );
}
static void draw2( void )
{
    
	//const XSize = 640, YSize = 480;
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	glOrtho (0, xsize, ysize, 0, 0, 1);
	glMatrixMode (GL_MODELVIEW);
	glDisable(GL_DEPTH_TEST);

//	printf("points: %f, %f \t",dbx[0],dby[0]);
//	printf("%f, %f \t",dbx[1],dby[1]);
//	printf("%f, %f \t",dbx[2],dby[2]);
//	printf("%f, %f \t",dbx[3],dby[3]);

	glColor3ub(0,255,0);
	glBegin(GL_LINE_LOOP);
	glVertex2f(dbx1[0],dby1[0]);
	glVertex2f(dbx1[1],dby1[1]);
	glVertex2f(dbx1[2],dby1[2]);
	glVertex2f(dbx1[3],dby1[3]);
	glEnd();

}
