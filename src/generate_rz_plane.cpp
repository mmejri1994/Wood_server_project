
#include <stdio.h>
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "Log.h"
#include "Branch.h"
using namespace WoodSeer;

// #define SAVE_CSV

int main(int argc, char *argv[]) {

    int num = 4000;
    if (argc > 1) {
        num = atoi(argv[1]);
    }

    unsigned int N = 128;
    double scale = 1.0/N;
    double step = 1.0/num;

    srand48(time(NULL));
    Log L(2+random()%3);
    while (num > 0) {
            double theta = -M_PI + num*2*step*M_PI;

            cv::Mat_<float> surface(N,N,0.0);
            cv::Mat_<float> rz_plane(N,N/2,0.0);

            // Generate surface 
            printf("%d\nSurface %f\n",num,L.getMaxSurface());
#ifdef SAVE_CSV
            FILE * fp=fopen("surf.csv","w");
#endif
            for (unsigned int i=0;i<N;i++) {
                for (unsigned int j=0;j<N;j++) {
                    double th = theta - M_PI/12 + i * M_PI / 6 / N;
                    double z = j * scale;
                    surface(j,i) = L.surface(th,z);
#ifdef SAVE_CSV
                    fprintf(fp,"%f,%f,%f,%f\n",th,z,surface(j,i),
                            surface(j,i)/L.getMaxSurface());
#endif
                }
            }
#ifdef SAVE_CSV
            fclose(fp);
#endif

            // Generate rz_plane 
            printf("RZ Plane %f\n",L.getMaxDensity());
#ifdef SAVE_CSV
            fp=fopen("rzplane.csv","w");
#endif
            for (unsigned int i=0;i<N/2;i++) {
                for (unsigned int j=0;j<N;j++) {
                    double r = i*scale;
                    double z = j * scale;
                    rz_plane(j,i) = L.density(r,theta,z);
#ifdef SAVE_CSV
                    fprintf(fp,"%f,%f,%f,%f\n",r,z,rz_plane(j,i),
                            rz_plane(j,i)/L.getMaxDensity());
#endif
                }
            }
#ifdef SAVE_CSV
            fclose(fp);
#endif

            cv::Mat_<uint8_t> surf_out(surface.size()); 
            cv::Mat_<uint8_t> rz_out(rz_plane.size()); 
            surface = 255*(surface - L.getRadius())/(L.getMaxSurface()-L.getRadius());
            rz_plane = 255*(rz_plane / L.getMaxDensity());
            surface.convertTo(surf_out,CV_8U);
            rz_plane.convertTo(rz_out,CV_8U);
            char tmp[128];
            sprintf(tmp,"%06d_",num);
            cv::imwrite(std::string(tmp)+"surface.png",surf_out);
            cv::imwrite(std::string(tmp)+"rzplane.png",rz_out);
            num -= 1;
        
    }

    return 0;
};


