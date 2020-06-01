
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
    int num = 1000;
    if (argc > 1) {
        num = atoi(argv[1]);
    }

    unsigned int N = 128;
    double scale_z = 0.5/N;
    double scale = 1.4/N;
    double nb_branches;
    srand48(time(NULL));
    while (num > 0) {
        nb_branches = 1+random()%4;
        Log L(nb_branches);
        printf("nb branches is %f\n",nb_branches);

        // Generate surface
        printf("%d\nSurface %f\n",num,L.getMaxSurface());
        /*char dirname[128];
        sprintf(dirname,"%06d",num);
        mkdir(dirname,0755);
        assert(chdir(dirname)==0);*/
#ifdef SAVE_CSV
        char volume[128];
        sprintf(volume,"volume_%03d.csv",num);
        FILE * fpv=fopen(volume,"w");
        fprintf(fpv,"x,y,z,density\n");
        char surf[128];
        sprintf(surf,"surface_%03d.csv",num);
        FILE * fp=fopen(surf,"w");
        fprintf(fp,"x,y,z,r\n");
        double output;
#endif
        for (unsigned int k=0;k<N;k++) {
            double z = k*scale_z;
            cv::Mat_<float> density(N,N,0.0);
            cv::Mat_<float> surface1(N,N,0.0);
            for (unsigned int i=0;i<N;i++) {
                double x = (double(i) - N/2.)*scale;
                for (unsigned int j=0;j<N;j++) {
                    double y = (double(j)-N/2.)*scale;
                    double r = hypot(x,y);
                    double theta = atan2(y,x);
                    double d = L.density(r,theta,z);
                    density(j,i) = d;
                    double r_s = L.surface(theta,z);
                    //printf("%f,\n",r_s);
                    output=2*r-r_s;
                    surface1(j,i)=output;
                    if (output<-0.023){surface1(j,i)=0;}
                    else if (output<0 && output-2e-2){surface1(j,i)=1;}
                    //else if (output<2e-3 && output>-2e-3){surface1(i,j)=0;}
                    //else if (r-0.5>-1e-3 && r-0.5>1e-3){surface1(i,j)=2;}
                    else{surface1(j,i)=0;}
                    fprintf(fpv,"%f,%f,%f,%f\n",x,y,z,density(j,i));
                    fprintf(fp,"%f,%f,%f,%f\n",x,y,z,surface1(j,i));
                    //printf("%f\n",output);
                   // if(r<0.3){
                   // surface1(i,j)=1;
                   // }
#ifdef SAVE_CSV
                    //fprintf(fpv,"%f,%f,%f,%f\n",x,y,z,d);
                    //fprintf(fp,"%f,%f,%f,%f\n",x,y,z,output);
#endif
                }
            }
        /*    cv::Mat_<uint8_t> d_out(density.size());
            density = 255*(density);
            density.convertTo(d_out,CV_8U);
            char tmp[128];
            sprintf(tmp,"layer_%03d.png",k);
            cv::imwrite(tmp,d_out);

        cv::Mat_<uint8_t> surf_out(surface1.size());
        //cv::minMaxLoc(surface1, &mini, &maxi);
        //cv::Mat_<uint8_t> output(surface.size());
        surface1 = 255*surface1 ;
        surface1.convertTo(surf_out,CV_8U);
        char surf_name[128];
        sprintf(surf_name,"surface_%03d.png",k);
        cv::imwrite(surf_name,surf_out);*/





        }
//#ifdef SAVE_CSV
        fclose(fpv);
        fclose(fp);
        // Generate rz_plane
        //printf("RZ Plane %f\n",L.getMaxDensity());
//        char surf[128];
//        sprintf(surf,"surface_%03d.csv",num);
//        FILE * fp=fopen(surf,"w");
//        fprintf(fp,"x,y,z,r\n");
//#endif
//        for (unsigned int k=0;k<N;k++) {
//            double z = k*scale;
//            cv::Mat_<float> surface1(N,N,0.0);
//            for (unsigned int i=0;i<N;i++) {
//                double x = (double(i) - N/2.)*scale;
//                for (unsigned int j=0;j<N;j++) {
//                    double y = (double(j)-N/2.)*scale;
//                    double theta = atan2(y,x);
//                    double r = L.surface(theta,z);
//                    double R = hypot(x,y);
//                    double output=3*R-r;
//                    surface1(i,j)=output;
//
//                    //fprintf(fp,"%f,%f,%f,%f\n",x,y,z,output);
//                    }
//             }
//
//
//        cv::Mat_<uint8_t> surf_out(surface1.size());
//        //cv::Mat_<uint8_t> output(surface.size());
//        surface1 = 255*(surface1 );
//        surface1.convertTo(surf_out,CV_8U);
//        char surf_name[128];
//        sprintf(surf_name,"surface_%03d.png",k);
//        cv::imwrite(surf_name,surf_out);
//       }
//
//#ifdef SAVE_CSV
//        fclose(fp);
//#endif
       /* cv::Mat_<uint8_t> surf_out(surface.size());
        cv::Mat_<uint8_t> output(surface.size());
        surface = 255*(surface - L.getRadius())/(L.getMaxSurface()-L.getRadius());
        surface.convertTo(surf_out,CV_8U);
        char surf_name[128];
        sprintf(surf_name,"surface_%03d.png",num);
        cv::imwrite(surf_name,output);*/
       // assert(chdir("..")==0);
        num -= 1;
    }

    return 0;
};


