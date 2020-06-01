
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "Log.h"
#include "Branch.h"
using namespace WoodSeer;


int main() {
    unsigned int N = 128;
    double scale = 1.0/N;

    srand48(time(NULL));
    Log L(2+random()%3);

    // Generate surface 
    printf("Surface %f\n",L.getMaxSurface());
    FILE * fp=fopen("volume.csv","w");
    fprintf(fp,"density\n");
    for (unsigned int i=0;i<N;i++) {
        double x = (double(i) - N/2.)*scale;
        for (unsigned int j=0;j<N;j++) {
            double y = (double(j)-N/2.)*scale;
            double r = hypot(x,y);
            double theta = atan2(y,x);
            for (unsigned int k=0;k<N;k++) {
                double z = k*scale;
                fprintf(fp,"%f\n",L.density(r,theta,z));
            }
        }
    }
    fclose(fp);

    // Generate rz_plane 
    printf("RZ Plane %f\n",L.getMaxDensity());
    cv::Mat_<float> surface(N,5*N,0.0);
    fp=fopen("surface.csv","w");
    fprintf(fp,"r\n");
    for (unsigned int i=0;i<5*N;i++) {
        double theta = -M_PI + i * 2 * M_PI / (5*N);
        for (unsigned int k=0;k<N;k++) {
            double z = k * scale;
            double r = L.surface(theta,z);
            surface(k,i) = r;
            double x = r * cos(theta);
            double y = r * sin(theta);
            fprintf(fp,"%f\n",r);
        }
    }
    fclose(fp);

    return 0;
};


