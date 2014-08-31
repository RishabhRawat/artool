#include <time.h>
#include <stdio.h>

/**                                                                             
 * sleep for `sec' seconds, without relying on the wall clock of time(2)        
 * or gettimeofday(2).                                                          
 *                                                                              
 * under ideal conditions is accurate to one microsecond. To get nanosecond     
 * accuracy, replace sleep()/usleep() with something with higher resolution     
 * like nanosleep() or ppoll().                                          
 */
int main()
{
        int sec = 1;
	struct timespec ts_start;
        struct timespec ts_end;

	clock_gettime(CLOCK_MONOTONIC, &ts_start);
	printf("%lu.%lu \n",(long)ts_start.tv_sec,(long)ts_start.tv_nsec);
        ts_end = ts_start;
}
