#include <iostream>
#include <math.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <maya.hpp>
using namespace std;
using namespace cv;

float dist(Point p1, Point p2)
{
 float r=(p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
 return sqrt(r);
}

float ang(Point p1, Point p2, Point p3){
  Point first, second;
  first.x=p1.x-p3.x; first.y=p1.y-p3.y;
  second.x=p2.x-p3.x; second.y=p2.y-p3.y;
  float mag1=sqrt(pow(first.x,2)+pow(first.y,2)), mag2=sqrt(pow(second.x,2)+pow(second.y,2));
  return (first.x*second.x+first.y*second.y)/(mag1*mag2);
}

void *gestureREC(void * argument){
	ComFlag=0;
	ComOutput=0;
	VideoCapture cap(0);
	Mat frame1, frame2;
	Point Start1, Start2, End1, End2; Start1.x=Start1.y=Start2.x=Start2.y=End1.x=End1.y=End2.x=End2.y=0;
	int moveoutput=0, n=10, counter=0;
	int i[n];
	Point avMc2; avMc2.x=avMc2.y=0;
	Point avMc1; avMc1.x=0; avMc1.y=0;
	namedWindow("Window1",CV_WINDOW_OPENGL);
// namedWindow("Window2",CV_WINDOW_AUTOSIZE);
 bool flag1_1= false, flag1_2= false, flag2_1= false, flag2_2= false, buffer=true;
 while(1)
 {
  if(counter==n) {counter=0; avMc1.x=avMc1.y=0; buffer=false;}
  cap.read(frame1);
  Mat dElement= getStructuringElement( MORPH_ELLIPSE, Size(5,5));
  Mat eElement= getStructuringElement( MORPH_ELLIPSE, Size(3,3));
  medianBlur(frame1,frame1,5);
  cvtColor(frame1,frame2,COLOR_BGR2YCrCb);
  inRange(frame2,Scalar(0,135,85),Scalar(255,180,135),frame2);
  dilate(frame2, frame2, eElement);
  erode(frame2, frame2, dElement);
//  medianBlur(frame2, frame2, 7);
 
 vector<vector<Point> > contours;
 vector<Vec4i> hierarchy, conDef;
 int largestContour=0, largestIndex;
 findContours( frame2, contours, hierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, Point(0,0)); //finding contour
 if (contours.size()>0)
 {
 for (int j=0; j< contours.size(); j++)// getting the largest contour
 {
  if( contours[j].size()> largestContour)
  {
   largestContour= contours[j].size();
   largestIndex=j;
  }
 }
 
  vector<Point> poly;
  vector<int>  poly1;
  drawContours( frame1, contours, largestIndex, Scalar(0,0,0), 1, 8);// drawing contour
  convexHull(contours[largestIndex], poly, 1, 1);// getting convex hull
  convexHull(contours[largestIndex], poly1, 1, 0);
//for drawing hull
/* for (int i=0; i< (poly.size()-1); i++)
 {
  line(frame, poly[i], poly[i+1], Scalar(0,0,0),3,8);
 }
 line(frame, poly[0], poly[poly.size()-1], Scalar(0,0,0),3,8);*/

 Moments mu[n];
 mu[counter]= moments(contours[largestIndex],false);
// double hu[7];
// HuMoments(mu,hu);
 Point2f mc[n];
 mc[counter]= Point2f(mu[counter].m10/mu[counter].m00,mu[counter].m01/mu[counter].m00);// getting centroid of area
 avMc1.x+=mc[counter].x; avMc1.y+=mc[counter].y;
 if(poly.size() > 3 )
  convexityDefects(contours[largestIndex], poly1, conDef);// convexity defects
 vector<Vec4i>:: iterator d= conDef.begin();
 Point list[50];
 i[counter]=0;
 while(d!=conDef.end())//getting deefect points
 {
  Vec4i& v=(*d);
  int startidx=v[0]; Point ptStart(contours[largestIndex][startidx]);
  int endidx=v[1]; Point ptEnd(contours[largestIndex][endidx]);
  int faridx=v[2]; Point ptFar(contours[largestIndex][faridx]);
  float depth=v[3]/256;

/*  line( frame, ptStart, ptFar, Scalar(200,182,98),3);
  line( frame, ptFar, ptEnd, Scalar(200,182,98),3);*/
  
  if((depth>40)&&(dist(ptStart,ptEnd)>20)&&(ang(ptStart,ptEnd,ptFar)>0))
  {
   list[3*i[counter]]=ptStart; 
   list[3*i[counter]+1]=ptFar;
   list[3*i[counter]+2]=ptEnd;
   i[counter]++;
  }
  
  d++;
 } 
 float avI=0;
 if(counter==(n-1))
 {
  for( int j=0; j<n; j++) {avI+=i[j];} 
  avI/=n; avMc1.x/=n; avMc1.y/=n;
// Point Start1, Start2, End1, End2; Start1.x=Start1.y=Start2.x=Start2.y=End1.x=End1.y=End2.x=End2.y=0; 
  
 if (avI>3) {flag1_1=true;} else flag1_1= false;
 if ((int (avI+0.5))==1) flag2_1=true; else flag2_1= false;
 
 if(flag1_1== true)// for moving right or left
 {
  Start1.x=avMc2.x; Start1.y= avMc2.y;
  End1.x=avMc1.x; End1.y=avMc1.y;
		line(frame1, Start1, End1, Scalar(225,0,0), 2, 8,0); 
		if((End1.x-Start1.x)>100){
			cout<<"right"<<endl;// RIGHT
			ComFlag=1;
			ComOutput=1;	
		}
		else if((End1.x-Start1.x)<(-100)){
			cout<<"left"<<endl;// LEFT
			ComFlag=1;
			ComOutput=2;
		}
 }

 if(flag2_1== true)
 {
  line(frame1, list[0], list[2], Scalar(25,65,200),2,8, 0);
  Start1.x=list[0].x; Start1.y=list[0].y; End1.x=list[2].x; End1.y=list[2].y;
  if(flag2_2==false) { Start2.x=list[0].x; Start2.y=list[0].y; End2.x=list[2].x; End2.y=list[2].y;}
	if(fabs(dist(Start1, End1)/dist(Start2, End2))>1.3){
		cout<<"zoom-out"<<endl;// ZOOM OUT
	}
  if(fabs(dist(Start1, End1)/dist(Start2, End2))<0.7){ cout<<"zoom-in"<<endl;}// ZOOM IN
  Start2.x=list[0].x; Start2.y=list[0].y; End2.x=list[2].x; End2.y=list[2].y;
  //cout<<dist(Start1, End1)<<" "<<dist(Start2, End2)<<endl;
  flag2_2=true;
 }
 avMc2.x=avMc1.x; avMc2.y=avMc1.y;
  } 
 circle(frame1, avMc2, 4, Scalar(0,225,0), -1, 8,0);
 } 
 imshow("Window1",frame1);
// imshow("Window2", frame2);
 if(waitKey(30)==27)
 {
  break;}
 counter++;
 }
	//cout<<moveoutput<<endl;
	return 0;
}

 
