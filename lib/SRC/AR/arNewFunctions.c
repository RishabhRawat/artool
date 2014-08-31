#include<AR/ar.h>
#include <AR/config.h>
#include <AR/param.h>
#include <AR/gsub.h>
#ifndef __APPLE__
#  include <GL/glut.h>
#  ifdef GL_VERSION_1_2
#    include <GL/glext.h>
#  endif
#else
#  include <GLUT/glut.h>
#  include <OpenGL/glext.h>
#endif

#ifndef GL_ABGR
#  define GL_ABGR GL_ABGR_EXT
#endif
#ifndef GL_BGRA
#  define GL_BGRA GL_BGRA_EXT
#endif
#ifndef GL_BGR
#  define GL_BGR GL_BGR_EXT
#endif
#ifndef GL_RGBA
#  define GL_RGBA GL_RGBA_EXT
#endif
#ifndef GL_RGB
#  define GL_RGB GL_RGB_EXT
#endif

#ifdef _WIN32
#  include <windows.h>
#  define put_zero(p,s) ZeroMemory(p, s)
#else
#  include <string.h>
#  define put_zero(p,s) memset((void *)p, 0, s)
#endif




static int check_square( int area, ARMarkerInfo2 *marker_info2, double factor );

static int get_vertex( int x_coord[], int y_coord[], int st, int ed,
                       double thresh, int vertex[], int *vnum );

static ARInt16 *labeling2_Clip( ARUint8 *image, int thresh,int *label_num, int **area, double **pos, int **clip, int **label_ref, int LorR );

static ARMarkerInfo2    marker_info2_test[AR_SQUARE_MAX];


static ARMarkerInfo2          *marker_info2;
static ARMarkerInfo           *wmarker_info;
static int                    wmarker_num = 0;

#define WORK_SIZE   1024*32

#define HARDCODED_BUFFER_WIDTH  1024
#define HARDCODED_BUFFER_HEIGHT 1024

static ARInt16      l_imageL[HARDCODED_BUFFER_WIDTH*HARDCODED_BUFFER_HEIGHT];
static ARInt16      l_imageR[HARDCODED_BUFFER_WIDTH*HARDCODED_BUFFER_HEIGHT];
/*****************************************************************************/
	
static int			arImXsize_Clip,arImYsize_Clip;

static int          workL[WORK_SIZE];
static int          workR[WORK_SIZE];
static int          work2L[WORK_SIZE*7];
static int          work2R[WORK_SIZE*7];	

static int          wlabel_numL;
static int          wlabel_numR;
static int          wareaL[WORK_SIZE];
static int          wareaR[WORK_SIZE];
static int          wclipL[WORK_SIZE*4];
static int          wclipR[WORK_SIZE*4];
static double       wposL[WORK_SIZE*2];
static double       wposR[WORK_SIZE*2];



void arCreateThresh(ARUint8 *image,ARUint8 *threshptr,int *thresh,int xmax,int ymax,int percentile)
{
        ARUint8 *ptr,*ptr2;
        ptr = threshptr;
        ptr2 = image;
        int hist[256];
        int i,sum;
        for(i = 0;i < 256; i++)
                hist[i]=0;
        sum = 0;
        for(i = 0;i<xmax*ymax;i++ ){
                sum = *(ptr2) + *(ptr2+1) + *(ptr2+2);
                hist[(int)(sum/3)]++;
                ptr2+=3;
        }
        sum = 0;
//      printf("sum = %d",sum);
        for(i = 0;i<256;i++){

                sum+=hist[i];

                if(sum>(xmax*ymax*percentile/100)){
                *thresh = i;
                break;
        }
        }
//      printf("sum = %d thresh = %d\n",sum,*thresh);
        ptr2 = image;
        for(i = 0;i<xmax*ymax;i++ ){
                sum = *(ptr2) + *(ptr2+1) + *(ptr2+2);
                hist[sum/3]++;
                if( sum  > *thresh*3 ){
                        *ptr = 255;
                        *(ptr+1) = 255;
                        *(ptr+2) = 255;
                }
                ptr+=3;
                ptr2+=3;
        }
}                          

void arPrintMat(double trans[3][4]){

printf("[%f %f %f %f ; %f %f %f %f ; %f %f %f %f ; ] \n ", trans[0][0],trans[0][1],trans[0][2],trans[0][3],trans[1][0],trans[1][1],trans[1][2],trans[1][3],trans[2][0],trans[2][1],trans[2][2],trans[2][3]);

}

void minmax (float dbx[4], float dby[4], double *xm,double *ym){

	int i;
	double minx,miny,maxx,maxy;
	minx = maxx = dbx[0];
	miny = maxy = dby[0];
	for(i = 1; i<4;i++){
	
	if(dbx[i]<minx) minx = dbx[i];
	
	if(dby[i]<miny) miny = dby[i];
	
	if(dbx[i]>maxx) maxx = dbx[i];

	if(dby[i]>maxy) miny = dby[i];
	}	
	xm[0] = minx;
	xm[1] = maxx;
	ym[0] = miny;
	ym[1] = maxy;
	return;

}

void scale(double *x,double *y){
	
	int i;	
	double xav,yav;
	xav=yav=0;

	for(i = 0;i < 4;i++)
	{
		xav+=x[i];
		yav+=y[i];
	}
	xav/=4;
	yav/=4;

	for(i = 0;i < 4;i++)
	{
		x[i] = 1.6*(x[i]-xav)+xav;
		y[i] = 1.6*(y[i]-yav)+yav;
	}

}








ARUint8* extract(ARUint8* threshptr, double *xm, double *ym, int xsize, int* w, int* h ){
	
	ARUint8* clipped;
	
	int minx,miny,maxx,maxy;
	int i,j;	

	minx = (int)xm[0];
	maxx = (int)xm[1];
	miny = (int)ym[0];
	maxy = (int)ym[1];
	*w = maxx-minx+1;
	*h = maxy-miny+1;
	arImXsize_Clip = *w;
	arImYsize_Clip = *h;
	clipped = malloc((*w)*(*h)*3);

	//Transfromation as (minx,miny) -> (0,0)
	for(j = miny;j <= maxy;j++)
	for(i = minx;i <= maxx;i++){
	
	
	*(clipped+3*((j-miny)*(*w)+(i-minx))) =*(threshptr+3*(j*xsize+i));
	*(clipped+3*((j-miny)*(*w)+(i-minx))+1) =*(threshptr+3*(j*xsize+i)+1);
	*(clipped+3*((j-miny)*(*w)+(i-minx))+2) =*(threshptr+3*(j*xsize+i)+2);

	}
	return clipped;

}

void argDispImageDrawPixels1( ARUint8 *image ,int gImXsize, int gImYsize)
{
    float    sx, sy;
    GLfloat  zoom;

        zoom = 1;
        sx = 0;
        sy = 320 - 0.5;
    glDisable(GL_TEXTURE_2D);
    glPixelZoom( zoom, -zoom);
    glRasterPos3f( sx, sy, -1.0 );

#if (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_ARGB)
#  ifdef AR_BIG_ENDIAN
    glDrawPixels( gImXsize, gImYsize, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV, image );
#  else
    glDrawPixels( gImXsize, gImYsize, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8, image );
#  endif        
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_ABGR)
    glDrawPixels( gImXsize, gImYsize, GL_ABGR, GL_UNSIGNED_BYTE, image );
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_BGRA)
    glDrawPixels( gImXsize, gImYsize, GL_BGRA, GL_UNSIGNED_BYTE, image );
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_BGR)
    glDrawPixels( gImXsize, gImYsize, GL_BGR, GL_UNSIGNED_BYTE, image );
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_RGBA)
    glDrawPixels( gImXsize, gImYsize, GL_RGBA, GL_UNSIGNED_BYTE, image );
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_RGB)
    glDrawPixels( gImXsize, gImYsize, GL_RGB, GL_UNSIGNED_BYTE, image );
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_MONO)
    glDrawPixels( gImXsize, gImYsize, GL_LUMINANCE, GL_UNSIGNED_BYTE, image );
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_2vuy)
#  ifdef AR_BIG_ENDIAN
    glDrawPixels( gImXsize, gImYsize, GL_YCBCR_422_APPLE, GL_UNSIGNED_SHORT_8_8_REV_APPLE, image );
#  else
    glDrawPixels( gImXsize, gImYsize, GL_YCBCR_422_APPLE, GL_UNSIGNED_SHORT_8_8_APPLE, image );
#  endif
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_yuvs)
#  ifdef AR_BIG_ENDIAN
        glDrawPixels( gImXsize, gImYsize, GL_YCBCR_422_APPLE, GL_UNSIGNED_SHORT_8_8_APPLE, image );
#  else
        glDrawPixels( gImXsize, gImYsize, GL_YCBCR_422_APPLE, GL_UNSIGNED_SHORT_8_8_REV_APPLE, image );
#  endif
#else
#  error Unknown default pixel format defined in config.h
#endif
}
int arDetect_inverted( ARUint8 *dataPtr, int thresh, double *dbx1,double *dby1 )
{
    ARInt16                *limage;
    ARMarkerInfo	   marker_info_1[1];
    int                    label_num;
    int                    *area, *clip, *label_ref;
    double                 *pos;
    int                    i;

    limage = labeling2_Clip( dataPtr, thresh,
                         &label_num, &area, &pos, &clip, &label_ref ,1);
    if( limage == 0 )    return -1;

    marker_info2 = arDetectMarker2_Clip( limage, label_num, label_ref,
                                    area, pos, clip, AR_AREA_MAX, AR_AREA_MIN,
                                    1.0, &wmarker_num);
    if( marker_info2 == 0 ) return -1;

	if( wmarker_num != 1)
	printf("error\n");


	arGetLine(marker_info2->x_coord, marker_info2->y_coord, marker_info2->coord_num, marker_info2->vertex,
                      marker_info_1[0].line, marker_info_1[0].vertex);

/*
	dbx1[0]=marker_info2->x_coord[marker_info2->vertex[0]];
	dbx1[1]=marker_info2->x_coord[marker_info2->vertex[1]];
	dbx1[2]=marker_info2->x_coord[marker_info2->vertex[2]];
	dbx1[3]=marker_info2->x_coord[marker_info2->vertex[3]];
	
	dby1[0]=marker_info2->y_coord[marker_info2->vertex[0]];
	dby1[1]=marker_info2->y_coord[marker_info2->vertex[1]];
	dby1[2]=marker_info2->y_coord[marker_info2->vertex[2]];
	dby1[3]=marker_info2->y_coord[marker_info2->vertex[3]];
  */
	
	dbx1[0]=marker_info_1[0].vertex[0][0];
	dbx1[1]=marker_info_1[0].vertex[1][0];
	dbx1[2]=marker_info_1[0].vertex[2][0];
	dbx1[3]=marker_info_1[0].vertex[3][0];
	
	dby1[0]=marker_info_1[0].vertex[0][1];
	dby1[1]=marker_info_1[0].vertex[1][1];
	dby1[2]=marker_info_1[0].vertex[2][1];
	dby1[3]=marker_info_1[0].vertex[3][1];
 


  
	return 0;
}

static ARInt16 *labeling2_Clip( ARUint8 *image, int thresh,int *label_num, int **area, double **pos, int **clip, int **label_ref, int LorR )
{
    ARUint8   *pnt;                     /*  image pointer       */
    ARInt16   *pnt1, *pnt2;             /*  image pointer       */
    int       *wk;                      /*  pointer for work    */
    int       wk_max;                   /*  work                */
    int       m,n;                      /*  work                */
    int       i,j,k;                    /*  for loop            */
    int       lxsize, lysize;
    int       poff;
    ARInt16   *l_image;
    int       *work, *work2;
    int       *wlabel_num;
    int       *warea;
    int       *wclip;
    double    *wpos;
	int		  pnt2_index;   // [tp]
	int		  thresht3 = thresh * 3;

	if (LorR) {
        l_image = &l_imageL[0];
        work    = &workL[0];
        work2   = &work2L[0];
        wlabel_num = &wlabel_numL;
        warea   = &wareaL[0];
        wclip   = &wclipL[0];
        wpos    = &wposL[0];
    } else {
        l_image = &l_imageR[0];
        work    = &workR[0];
        work2   = &work2R[0];
        wlabel_num = &wlabel_numR;
        warea   = &wareaR[0];
        wclip   = &wclipR[0];
        wpos    = &wposR[0];
    }

    if (arImageProcMode == AR_IMAGE_PROC_IN_HALF) {
        lxsize = arImXsize_Clip / 2;
        lysize = arImYsize_Clip / 2;
    } else {
        lxsize = arImXsize_Clip;
        lysize = arImYsize_Clip;
    }

    pnt1 = &l_image[0]; // Leftmost pixel of top row of image.
    pnt2 = &l_image[(lysize - 1)*lxsize]; // Leftmost pixel of bottom row of image.

// 4x loop unrolling
	for (i = 0; i < lxsize - (lxsize%4); i += 4) {
        *(pnt1++) = *(pnt2++) = 0;
        *(pnt1++) = *(pnt2++) = 0;
        *(pnt1++) = *(pnt2++) = 0;
        *(pnt1++) = *(pnt2++) = 0;
    }
    pnt1 = &l_image[0]; // Leftmost pixel of top row of image.
    pnt2 = &l_image[lxsize - 1]; // Rightmost pixel of top row of image.

// 4x loop unrolling
    for (i = 0; i < lysize - (lysize%4); i += 4) {
		*pnt1 = *pnt2 = 0;
        pnt1 += lxsize;
        pnt2 += lxsize;

		*pnt1 = *pnt2 = 0;
        pnt1 += lxsize;
        pnt2 += lxsize;

		*pnt1 = *pnt2 = 0;
        pnt1 += lxsize;
        pnt2 += lxsize;

		*pnt1 = *pnt2 = 0;
        pnt1 += lxsize;
        pnt2 += lxsize;
    }

    wk_max = 0;
    pnt2 = &(l_image[lxsize+1]);
    if (arImageProcMode == AR_IMAGE_PROC_IN_HALF) {
        pnt = &(image[(arImXsize_Clip*2+2)*AR_PIX_SIZE_DEFAULT]);
        poff = AR_PIX_SIZE_DEFAULT*2;
    } else {
        pnt = &(image[(arImXsize_Clip+1)*AR_PIX_SIZE_DEFAULT]);
        poff = AR_PIX_SIZE_DEFAULT;
    }
    for (j = 1; j < lysize - 1; j++, pnt += poff*2, pnt2 += 2) {
        for(i = 1; i < lxsize-1; i++, pnt+=poff, pnt2++) {
#if (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_ARGB)
            if( *(pnt+1) + *(pnt+2) + *(pnt+3) > thresht3 )
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_ABGR)
            if( *(pnt+1) + *(pnt+2) + *(pnt+3) > thresht3 )
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_BGRA)
            if( *(pnt+0) + *(pnt+1) + *(pnt+2) > thresht3 )
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_BGR)
            if( *(pnt+0) + *(pnt+1) + *(pnt+2) > thresht3 )
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_RGBA)
            if( *(pnt+0) + *(pnt+1) + *(pnt+2) > thresht3 )
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_RGB)
            if( *(pnt+0) + *(pnt+1) + *(pnt+2) > thresht3 )
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_MONO)
			if( *(pnt) > thresh )
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_2vuy)
			if( *(pnt+1) > thresh )
#elif (AR_DEFAULT_PIXEL_FORMAT == AR_PIXEL_FORMAT_yuvs)
			if( *(pnt+0) > thresh )
#else
#  error Unknown default pixel format defined in config.h
#endif
			{
                pnt1 = &(pnt2[-lxsize]);
                if( *pnt1 > 0 ) {
                    *pnt2 = *pnt1;

		// OPTIMIZED CODE [tp]
		// ((*pnt2)-1)*7 should be treated as constant, since
		//  work2[n] (n=0..xsize*ysize) cannot overwrite (*pnt2)
		pnt2_index = ((*pnt2)-1) * 7;
                work2[pnt2_index+0]++;
                work2[pnt2_index+1]+= i;
                work2[pnt2_index+2]+= j;
                work2[pnt2_index+6] = j;
		// --------------------------------

                }
                else if( *(pnt1+1) > 0 ) {
                    if( *(pnt1-1) > 0 ) {
                        m = work[*(pnt1+1)-1];
                        n = work[*(pnt1-1)-1];
                        if( m > n ) {
                            *pnt2 = n;
                            wk = &(work[0]);
                            for(k = 0; k < wk_max; k++) {
                                if( *wk == m ) *wk = n;
                                wk++;
                            }
                        }
                        else if( m < n ) {
                            *pnt2 = m;
                            wk = &(work[0]);
                            for(k = 0; k < wk_max; k++) {
                                if( *wk == n ) *wk = m;
                                wk++;
                            }
                        }
                        else *pnt2 = m;

#ifndef USE_OPTIMIZATIONS
						// ORIGINAL CODE
						work2[((*pnt2)-1)*7+0] ++;
                        work2[((*pnt2)-1)*7+1] += i;
                        work2[((*pnt2)-1)*7+2] += j;
                        work2[((*pnt2)-1)*7+6] = j;
#else
						// PERFORMANCE OPTIMIZATION:
						pnt2_index = ((*pnt2)-1) * 7;
						work2[pnt2_index+0]++;
						work2[pnt2_index+1]+= i;
						work2[pnt2_index+2]+= j;
						work2[pnt2_index+6] = j;
#endif

                    }
                    else if( *(pnt2-1) > 0 ) {
                        m = work[*(pnt1+1)-1];
                        n = work[*(pnt2-1)-1];
                        if( m > n ) {
                            *pnt2 = n;
                            wk = &(work[0]);
                            for(k = 0; k < wk_max; k++) {
                                if( *wk == m ) *wk = n;
                                wk++;
                            }
                        }
                        else if( m < n ) {
                            *pnt2 = m;
                            wk = &(work[0]);
                            for(k = 0; k < wk_max; k++) {
                                if( *wk == n ) *wk = m;
                                wk++;
                            }
                        }
                        else *pnt2 = m;

#ifndef USE_OPTIMIZATIONS
						// ORIGINAL CODE
                        work2[((*pnt2)-1)*7+0] ++;
                        work2[((*pnt2)-1)*7+1] += i;
                        work2[((*pnt2)-1)*7+2] += j;
#else
						// PERFORMANCE OPTIMIZATION:
						pnt2_index = ((*pnt2)-1) * 7;
						work2[pnt2_index+0]++;
						work2[pnt2_index+1]+= i;
						work2[pnt2_index+2]+= j;
#endif

                    }
                    else {
                        *pnt2 = *(pnt1+1);

#ifndef USE_OPTIMIZATIONS
						// ORIGINAL CODE
                        work2[((*pnt2)-1)*7+0] ++;
                        work2[((*pnt2)-1)*7+1] += i;
                        work2[((*pnt2)-1)*7+2] += j;
                        if( work2[((*pnt2)-1)*7+3] > i ) work2[((*pnt2)-1)*7+3] = i;
                        work2[((*pnt2)-1)*7+6] = j;
#else
						// PERFORMANCE OPTIMIZATION:
						pnt2_index = ((*pnt2)-1) * 7;
						work2[pnt2_index+0]++;
						work2[pnt2_index+1]+= i;
						work2[pnt2_index+2]+= j;
                        if( work2[pnt2_index+3] > i ) work2[pnt2_index+3] = i;
						work2[pnt2_index+6] = j;
#endif
                    }
                }
                else if( *(pnt1-1) > 0 ) {
                    *pnt2 = *(pnt1-1);

#ifndef USE_OPTIMIZATIONS
						// ORIGINAL CODE
                    work2[((*pnt2)-1)*7+0] ++;
                    work2[((*pnt2)-1)*7+1] += i;
                    work2[((*pnt2)-1)*7+2] += j;
                    if( work2[((*pnt2)-1)*7+4] < i ) work2[((*pnt2)-1)*7+4] = i;
                    work2[((*pnt2)-1)*7+6] = j;
#else
					// PERFORMANCE OPTIMIZATION:
					pnt2_index = ((*pnt2)-1) * 7;
					work2[pnt2_index+0]++;
					work2[pnt2_index+1]+= i;
					work2[pnt2_index+2]+= j;
                    if( work2[pnt2_index+4] < i ) work2[pnt2_index+4] = i;
					work2[pnt2_index+6] = j;
#endif
                }
                else if( *(pnt2-1) > 0) {
                    *pnt2 = *(pnt2-1);

#ifndef USE_OPTIMIZATIONS
						// ORIGINAL CODE
                    work2[((*pnt2)-1)*7+0] ++;
                    work2[((*pnt2)-1)*7+1] += i;
                    work2[((*pnt2)-1)*7+2] += j;
                    if( work2[((*pnt2)-1)*7+4] < i ) work2[((*pnt2)-1)*7+4] = i;
#else
					// PERFORMANCE OPTIMIZATION:
					pnt2_index = ((*pnt2)-1) * 7;
					work2[pnt2_index+0]++;
					work2[pnt2_index+1]+= i;
					work2[pnt2_index+2]+= j;
                    if( work2[pnt2_index+4] < i ) work2[pnt2_index+4] = i;
#endif
				}
                else {
                    wk_max++;
                    if( wk_max > WORK_SIZE ) {
                        return(0);
                    }
                    work[wk_max-1] = *pnt2 = wk_max;
                    work2[(wk_max-1)*7+0] = 1;
                    work2[(wk_max-1)*7+1] = i;
                    work2[(wk_max-1)*7+2] = j;
                    work2[(wk_max-1)*7+3] = i;
                    work2[(wk_max-1)*7+4] = i;
                    work2[(wk_max-1)*7+5] = j;
                    work2[(wk_max-1)*7+6] = j;
                }
            }
            else {
                *pnt2 = 0;
            }
        }
        if (arImageProcMode == AR_IMAGE_PROC_IN_HALF) pnt += arImXsize_Clip*AR_PIX_SIZE_DEFAULT;
    }

    j = 1;
    wk = &(work[0]);
    for(i = 1; i <= wk_max; i++, wk++) {
        *wk = (*wk==i)? j++: work[(*wk)-1];
    }
    *label_num = *wlabel_num = j - 1;
    if( *label_num == 0 ) {
        return( l_image );
    }

    put_zero( (ARUint8 *)warea, *label_num *     sizeof(int) );
    put_zero( (ARUint8 *)wpos,  *label_num * 2 * sizeof(double) );
    for(i = 0; i < *label_num; i++) {
        wclip[i*4+0] = lxsize;
        wclip[i*4+1] = 0;
        wclip[i*4+2] = lysize;
        wclip[i*4+3] = 0;
    }
    for(i = 0; i < wk_max; i++) {
        j = work[i] - 1;
        warea[j]    += work2[i*7+0];
        wpos[j*2+0] += work2[i*7+1];
        wpos[j*2+1] += work2[i*7+2];
        if( wclip[j*4+0] > work2[i*7+3] ) wclip[j*4+0] = work2[i*7+3];
        if( wclip[j*4+1] < work2[i*7+4] ) wclip[j*4+1] = work2[i*7+4];
        if( wclip[j*4+2] > work2[i*7+5] ) wclip[j*4+2] = work2[i*7+5];
        if( wclip[j*4+3] < work2[i*7+6] ) wclip[j*4+3] = work2[i*7+6];
    }

    for( i = 0; i < *label_num; i++ ) {
        wpos[i*2+0] /= warea[i];
        wpos[i*2+1] /= warea[i];
    }

    *label_ref = work;
    *area      = warea;
    *pos       = wpos;
    *clip      = wclip;
    return (l_image);
}

ARMarkerInfo2 *arDetectMarker2_Clip( ARInt16 *limage, int label_num, int *label_ref,int *warea, double *wpos, int *wclip,int area_max, int area_min, double factor, int *marker_num )
{
    ARMarkerInfo2     *pm;
    int               xsize, ysize;
    int               marker_num2;
    int               i, j, ret;
    double            d;

    if( arImageProcMode == AR_IMAGE_PROC_IN_HALF ) {
        area_min /= 4;
        area_max /= 4;
        xsize = arImXsize_Clip / 2;
        ysize = arImYsize_Clip / 2;
    }
    else {
        xsize = arImXsize_Clip;
        ysize = arImYsize_Clip;
    }
    marker_num2 = 0;
    for(i=0; i<label_num; i++ ) {
        if( warea[i] < area_min || warea[i] > area_max ) continue;
        if( wclip[i*4+0] == 1 || wclip[i*4+1] == xsize-2 ) continue;
        if( wclip[i*4+2] == 1 || wclip[i*4+3] == ysize-2 ) continue;

        ret = arGetContour_1( limage, label_ref, i+1,
                            &(wclip[i*4]), &(marker_info2_test[marker_num2]));
        if( ret < 0 ) continue;

        ret = check_square( warea[i], &(marker_info2_test[marker_num2]), factor );
        if( ret < 0 ) continue;

        marker_info2_test[marker_num2].area   = warea[i];
        marker_info2_test[marker_num2].pos[0] = wpos[i*2+0];
        marker_info2_test[marker_num2].pos[1] = wpos[i*2+1];
        marker_num2++;
        if( marker_num2 == AR_SQUARE_MAX ) break;
    }

    for( i=0; i < marker_num2; i++ ) {
        for( j=i+1; j < marker_num2; j++ ) {
            d = (marker_info2_test[i].pos[0] - marker_info2_test[j].pos[0])
              * (marker_info2_test[i].pos[0] - marker_info2_test[j].pos[0])
              + (marker_info2_test[i].pos[1] - marker_info2_test[j].pos[1])
              * (marker_info2_test[i].pos[1] - marker_info2_test[j].pos[1]);
            if( marker_info2_test[i].area > marker_info2_test[j].area ) {
                if( d < marker_info2_test[i].area / 4 ) {
                    marker_info2_test[j].area = 0;
                }
            }
            else {
                if( d < marker_info2_test[j].area / 4 ) {
                    marker_info2_test[i].area = 0;
                }
            }
        }
    }
    for( i=0; i < marker_num2; i++ ) {
        if( marker_info2_test[i].area == 0.0 ) {
            for( j=i+1; j < marker_num2; j++ ) {
                marker_info2_test[j-1] = marker_info2_test[j];
            }
            marker_num2--;
        }
    }

    if( arImageProcMode == AR_IMAGE_PROC_IN_HALF ) {
        pm = &(marker_info2_test[0]);
        for( i = 0; i < marker_num2; i++ ) {
            pm->area *= 4;
            pm->pos[0] *= 2.0;
            pm->pos[1] *= 2.0;
            for( j = 0; j< pm->coord_num; j++ ) {
                pm->x_coord[j] *= 2;
                pm->y_coord[j] *= 2;
            }
            pm++;
        }
    }

    *marker_num = marker_num2;
    return( &(marker_info2_test[0]) );
}

int arGetContour_1( ARInt16 *limage, int *label_ref,
                  int label, int clip[4], ARMarkerInfo2 *marker_info2 )
{
    static int      xdir[8] = { 0, 1, 1, 1, 0,-1,-1,-1};
    static int      ydir[8] = {-1,-1, 0, 1, 1, 1, 0,-1};
    static int      wx[AR_CHAIN_MAX];
    static int      wy[AR_CHAIN_MAX];
    ARInt16         *p1;
    int             xsize, ysize;
    int             sx, sy, dir;
    int             dmax, d, v1;
    int             i, j;

    if( arImageProcMode == AR_IMAGE_PROC_IN_HALF ) {
        xsize = arImXsize_Clip / 2;
        ysize = arImYsize_Clip / 2;
    }
    else {
        xsize = arImXsize_Clip;
        ysize = arImYsize_Clip;
    }
    j = clip[2];
    p1 = &(limage[j*xsize+clip[0]]);
    for( i = clip[0]; i <= clip[1]; i++, p1++ ) {
        if( *p1 > 0 && label_ref[(*p1)-1] == label ) {
            sx = i; sy = j; break;
        }
    }
    if( i > clip[1] ) {
        printf("??? 1\n"); return(-1);
    }

    marker_info2->coord_num = 1;
    marker_info2->x_coord[0] = sx;
    marker_info2->y_coord[0] = sy;
    dir = 5;
    for(;;) {
        p1 = &(limage[marker_info2->y_coord[marker_info2->coord_num-1] * xsize
                    + marker_info2->x_coord[marker_info2->coord_num-1]]);
        dir = (dir+5)%8;
        for(i=0;i<8;i++) {
            if( p1[ydir[dir]*xsize+xdir[dir]] > 0 ) break;
            dir = (dir+1)%8;
        }
        if( i == 8 ) {
            printf("??? 2\n"); return(-1);
        }
        marker_info2->x_coord[marker_info2->coord_num]
            = marker_info2->x_coord[marker_info2->coord_num-1] + xdir[dir];
        marker_info2->y_coord[marker_info2->coord_num]
            = marker_info2->y_coord[marker_info2->coord_num-1] + ydir[dir];
        if( marker_info2->x_coord[marker_info2->coord_num] == sx
         && marker_info2->y_coord[marker_info2->coord_num] == sy ) break;
        marker_info2->coord_num++;
        if( marker_info2->coord_num == AR_CHAIN_MAX-1 ) {
            printf("??? 3\n"); return(-1);
        }
    }

    dmax = 0;
    for(i=1;i<marker_info2->coord_num;i++) {
        d = (marker_info2->x_coord[i]-sx)*(marker_info2->x_coord[i]-sx)
          + (marker_info2->y_coord[i]-sy)*(marker_info2->y_coord[i]-sy);
        if( d > dmax ) {
            dmax = d;
            v1 = i;
        }
    }

    for(i=0;i<v1;i++) {
        wx[i] = marker_info2->x_coord[i];
        wy[i] = marker_info2->y_coord[i];
    }
    for(i=v1;i<marker_info2->coord_num;i++) {
        marker_info2->x_coord[i-v1] = marker_info2->x_coord[i];
        marker_info2->y_coord[i-v1] = marker_info2->y_coord[i];
    }
    for(i=0;i<v1;i++) {
        marker_info2->x_coord[i-v1+marker_info2->coord_num] = wx[i];
        marker_info2->y_coord[i-v1+marker_info2->coord_num] = wy[i];
    }
    marker_info2->x_coord[marker_info2->coord_num] = marker_info2->x_coord[0];
    marker_info2->y_coord[marker_info2->coord_num] = marker_info2->y_coord[0];
    marker_info2->coord_num++;

    return 0;
}



static int check_square( int area, ARMarkerInfo2 *marker_info2, double factor )
{
    int             sx, sy;
    int             dmax, d, v1;
    int             vertex[10], vnum;
    int             wv1[10], wvnum1, wv2[10], wvnum2, v2;
    double          thresh;
    int             i;


    dmax = 0;
    v1 = 0;
    sx = marker_info2->x_coord[0];
    sy = marker_info2->y_coord[0];
    for(i=1;i<marker_info2->coord_num-1;i++) {
        d = (marker_info2->x_coord[i]-sx)*(marker_info2->x_coord[i]-sx)
          + (marker_info2->y_coord[i]-sy)*(marker_info2->y_coord[i]-sy);
        if( d > dmax ) {
            dmax = d;
            v1 = i;
        }
    }

    thresh = (area/0.75) * 0.01 * factor;
    vnum = 1;
    vertex[0] = 0;
    wvnum1 = 0;
    wvnum2 = 0;
    if( get_vertex(marker_info2->x_coord, marker_info2->y_coord, 0,  v1,
                   thresh, wv1, &wvnum1) < 0 ) {
        return(-1);
    }
    if( get_vertex(marker_info2->x_coord, marker_info2->y_coord,
                   v1,  marker_info2->coord_num-1, thresh, wv2, &wvnum2) < 0 ) {
        return(-1);
    }

    if( wvnum1 == 1 && wvnum2 == 1 ) {
        vertex[1] = wv1[0];
        vertex[2] = v1;
        vertex[3] = wv2[0];
    }
    else if( wvnum1 > 1 && wvnum2 == 0 ) {
        v2 = v1 / 2;
        wvnum1 = wvnum2 = 0;
        if( get_vertex(marker_info2->x_coord, marker_info2->y_coord,
                       0,  v2, thresh, wv1, &wvnum1) < 0 ) {
            return(-1);
        }
        if( get_vertex(marker_info2->x_coord, marker_info2->y_coord,
                       v2,  v1, thresh, wv2, &wvnum2) < 0 ) {
            return(-1);
        }
        if( wvnum1 == 1 && wvnum2 == 1 ) {
            vertex[1] = wv1[0];
            vertex[2] = wv2[0];
            vertex[3] = v1;
        }
        else {
            return(-1);
        }
    }
    else if( wvnum1 == 0 && wvnum2 > 1 ) {
        v2 = (v1 + marker_info2->coord_num-1) / 2;
        wvnum1 = wvnum2 = 0;
        if( get_vertex(marker_info2->x_coord, marker_info2->y_coord,
                   v1, v2, thresh, wv1, &wvnum1) < 0 ) {
            return(-1);
        }
        if( get_vertex(marker_info2->x_coord, marker_info2->y_coord,
                   v2, marker_info2->coord_num-1, thresh, wv2, &wvnum2) < 0 ) {
            return(-1);
        }
        if( wvnum1 == 1 && wvnum2 == 1 ) {
            vertex[1] = v1;
            vertex[2] = wv1[0];
            vertex[3] = wv2[0];
        }
        else {
            return(-1);
        }
    }
    else {
        return(-1);
    }

    marker_info2->vertex[0] = vertex[0];
    marker_info2->vertex[1] = vertex[1];
    marker_info2->vertex[2] = vertex[2];
    marker_info2->vertex[3] = vertex[3];
    marker_info2->vertex[4] = marker_info2->coord_num-1;

    return(0);
}

static int get_vertex( int x_coord[], int y_coord[], int st,  int ed,
                       double thresh, int vertex[], int *vnum)
{
    double   d, dmax;
    double   a, b, c;
    int      i, v1;

    a = y_coord[ed] - y_coord[st];
    b = x_coord[st] - x_coord[ed];
    c = x_coord[ed]*y_coord[st] - y_coord[ed]*x_coord[st];
    dmax = 0;
    for(i=st+1;i<ed;i++) {
        d = a*x_coord[i] + b*y_coord[i] + c;
        if( d*d > dmax ) {
            dmax = d*d;
            v1 = i;
        }
    }
    if( dmax/(a*a+b*b) > thresh ) {
        if( get_vertex(x_coord, y_coord, st,  v1, thresh, vertex, vnum) < 0 )
            return(-1);

        if( (*vnum) > 5 ) return(-1);
        vertex[(*vnum)] = v1;
        (*vnum)++;

        if( get_vertex(x_coord, y_coord, v1,  ed, thresh, vertex, vnum) < 0 )
            return(-1);
    }

    return(0);
}
