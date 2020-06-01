
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

//#define SAVE_CSV

int main(int argc, char *argv[]) {

    //FILE * fpb = fopen("branches.csv","w");
    //fprintf(fpb,"nb_branches\n","w");
for (unsigned int s=6;s<7;s++){

    int num = 200;
    if (argc > 1) {
        num = atoi(argv[1]);
    }



    while (num > 0) {
            srand48(time(NULL));
            unsigned int nb_branches=s;
            Log L(nb_branches);
           // fprintf(fpb,"%0d\n",nb_branches);
            unsigned int nb_slices=64;
            unsigned int N = 64;
            double scale = 1.0/N;
        // Generate surface
       // printf("%d\nSurface %f\n",num,L.getMaxSurface());
#ifdef SAVE_CSV
        //FILE * fp=fopen("volume.csv","w");
        //fprintf(fp,"x,y,z,density\n");
#endif
        char dirname[128];
        sprintf(dirname,"nb_branches_%02d_%06d",s,num);
        mkdir(dirname,0755);
        assert(chdir(dirname)==0);
        for (unsigned int k=0;k<nb_slices;k++){
             
            double theta = -M_PI + k * 2 * M_PI / nb_slices;
            cv::Mat_<float> density(N,N/2,0.0);
            for (unsigned int i=0;i<N/2;i++) {
                for (unsigned int j=0;j<N;j++) {
                    double r = i*scale;
                    double z = j * scale;
                    double d = L.density(r,theta,z);
                    density(j,i) = d;
#ifdef SAVE_CSV
       //             fprintf(fp,"%f,%f,%f,%f\n",x,y,z,d);
#endif
                }
            }
            cv::Mat_<uint8_t> d_out(density.size());
            density = 255*(density / L.getMaxDensity());
            density.convertTo(d_out,CV_8U);
            char tmp1[128];
            sprintf(tmp1,"layer_%06d.png",k);
            cv::imwrite(tmp1,d_out);
        
#ifdef SAVE_CSV
        //fclose(fp);
#endif

        // Generate rz_plane
        //printf("RZ Plane %f\n",L.getMaxDensity());
        cv::Mat_<float> surface(N,N,0.0);
#ifdef SAVE_CSV

#endif
            //printf(tmp1);
            for (unsigned int i=0;i<N;i++) {
                for (unsigned int j=0;j<N;j++) {
                    double th = theta - M_PI/6 + i * M_PI / 3 / N;
                    double z = j * scale;
                    surface(j,i) = L.surface(th,z);
#ifdef SAVE_CSV

#endif
            }
        }
#ifdef SAVE_CSV
        //fclose(fp);
#endif
        
        cv::Mat_<uint8_t> surf_out(surface.size());
        surface = 255*surface;
        surface.convertTo(surf_out,CV_8U);
        char tmp2[128];
        sprintf(tmp2,"surface_%06d.png",k);
        cv::imwrite(tmp2,surf_out);
        
       }
//-----------------------------------------3D-Woodseer----------------------------

        // Generate surface
        printf("%d\nSurface %f\n",num,L.getMaxSurface());

        char vox[128];
        sprintf(vox,"volume_%03d.csv",num);
        FILE * fpv=fopen(vox,"w"); 
        fprintf(fpv,"density\n"); 
        for (unsigned int k=0;k<N;k++) {
            double z = k*scale;
            cv::Mat_<float> density(N,N,0.0);
            for (unsigned int i=0;i<N;i++) {
                double x = (double(i) - N/2.)*scale;
                for (unsigned int j=0;j<N;j++) {
                    double y = (double(j)-N/2.)*scale;
                    double r = hypot(x,y);
                    double theta = atan2(y,x);
                    double d = L.density(r,theta,z);
                    density(j,i) = d;

                    fprintf(fpv,"%f\n",d);

                }
            }

        }


        fclose(fpv);


        // Generate rz_plane
        //printf("RZ Plane %f\n",L.getMaxDensity());
        //cv::Mat_<float> surface(N,N*N,0.0);

        char surf[128];
        sprintf(surf,"surface_%03d.csv",num);
        FILE * fp=fopen(surf,"w");
        fprintf(fp,"r\n");

        for (unsigned int k=0;k<N;k++) {
            double z = k*scale;
            for (unsigned int i=0;i<N;i++) {
                double x = (double(i) - N/2.)*scale;
                for (unsigned int j=0;j<N;j++) {
                    double y = (double(j)-N/2.)*scale;
                    double theta = atan2(y,x);
                    double r = L.surface(theta,z);
                    double R = hypot(x,y);
                    double output=R-r;
                    fprintf(fp,"%f\n",output);
                    }
             }
       }

#ifdef SAVE_CSV
        fclose(fp);
#endif




        num -= 1;
        assert(chdir("..")==0);

    }

   }
//    fclose(fpb);
    return 0;

};



