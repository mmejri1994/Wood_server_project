
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "Log.h"
#include "Branch.h"
using namespace WoodSeer;

#define SAVE_CSV

int main(int argc, char *argv[]) {
    int num = 1;
    if (argc > 1) {
        num = atoi(argv[1]);
    }

    unsigned int N = 128;
    double scale = 1.0/N;

    //srand48(time(NULL));
    while (num > 0) {
        Log L(1+random()%8);
        // Generate surface
        printf("%d\nSurface %f\n",num,L.getMaxSurface());
#ifdef SAVE_CSV
        char volume[128];
        sprintf(volume,"volume_%03d.csv",num);
        FILE * fpv=fopen(volume,"w");
        fprintf(fpv,"x,y,z,density\n");
#endif
        for (unsigned int k=0;k<N;k++) {
            double z = k*scale;
            //cv::Mat_<float> density(N,N,0.0);
            for (unsigned int i=0;i<N;i++) {
                double x = (double(i) - N/2.)*scale;
                for (unsigned int j=0;j<N;j++) {
                    double y = (double(j)-N/2.)*scale;
                    double r = hypot(x,y);
                    double theta = atan2(y,x);
                    double d = L.density(r,theta,z);
                    //double d_base = L.density_simple(r,theta,z);
#ifdef SAVE_CSV
                    fprintf(fpv,"%f,%f,%f,%f\n",x,y,z,d);
#endif
                }
            }

        }
#ifdef SAVE_CSV
        fclose(fpv);
#endif

        // Generate rz_plane
        //printf("RZ Plane %f\n",L.getMaxDensity());
#ifdef SAVE_CSV
        char surf[128];
        sprintf(surf,"surface_%03d.csv",num);
        FILE * fp=fopen(surf,"w");
        fprintf(fp,"x,y,z,r\n");
#endif
        for (unsigned int k=0;k<N;k++) {
            double z = k*scale;
            for (unsigned int i=0;i<N;i++) {
                double x = (double(i) - N/2.)*scale;
                for (unsigned int j=0;j<N;j++) {
                    double y = (double(j)-N/2.)*scale;
                    double theta = atan2(y,x);
                    double r = L.surface(theta,z);
                    double R = hypot(x,y);
                    double output=2*R-r;
                    fprintf(fp,"%f,%f,%f,%f\n",x,y,z,output);
                    }
             }
       }

#ifdef SAVE_CSV
        fclose(fp);
#endif

        num -= 1;
    }

    return 0;
};


