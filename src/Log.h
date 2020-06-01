#ifndef WOODSEER_LOG_H
#define WOODSEER_LOG_H

#include <stdlib.h>
#include <math.h>
#include <vector>
#include "Branch.h"

namespace WoodSeer {
    class Log {
        protected:
            Branch b;
            double radius;
            std::vector<Branch> branches;
            double density_out;
            double density_base;
            double density_noise;
            double surface_noise;
            double surface_amplitude;
            double density_amplitude;
            double eps ;
            double nb_cercles;
        public:
            Log(unsigned int n, double r=0.5) {
                b= Branch(true,0.0);
                surface_noise = 0.00;
                density_noise = 0.00;
                density_base = 0.05;
                density_out = 0.05;
                radius = r+0.1;
                eps=0.18;
                nb_cercles = 4+random()%4;
                for (unsigned int i=0;i<n;i++) {
                    branches.push_back(Branch(true,
                                -M_PI + i*2*M_PI/n + Branch::randt()*M_PI/n));
                }
            }

            double getRadius() const {
                return radius;
            }

            double getMaxSurface() const {
                double rmax=0;
                for (size_t i=0;i<branches.size();i++) {
                    double r = branches[i].getMaxSurface(radius);
                    if (r > rmax) rmax = r;
                }
                return surface_noise + rmax;
            }

            double getMaxDensity() const {
                double dmax=0;
                for (size_t i=0;i<branches.size();i++) {
                    double d = branches[i].getMaxDensity();
                    if (d > dmax) dmax = d;
                }
                return std::max(density_base,density_out) + density_noise + dmax;
            }

            double density(double r, double theta, double z) {
                double dmax = 0;

                for (size_t i=0;i<branches.size();i++) {
                    double d = branches[i].density(r,theta,z);

                    if (d > dmax) dmax = d;
                }

                double R = radius*b.radius(r,theta,z);
                //double Rc = radius*b.contour_cercle(r,theta,z);
                    if (r > R*(0.5-eps)) {
                    return std::max(density_out+ Branch::randt()*density_noise,0.0);
                    }
             /*       double ext_R=6*Rc*(0.5-eps)/9;
                    double int_R=0.95*ext_R;
                    for(unsigned int i=0;i<nb_cercles;i++){
                        if(r<ext_R and r>int_R){
                            //printf("done R %03f\n",r);
                            return std::max(0.4 + density_base ,0.0);
                             }
                     ext_R = ext_R*0.75;
                     int_R = int_R*0.75;

                    } */
                dmax = 0;
                for (size_t i=0;i<branches.size();i++) {
                    double d = branches[i].density(r,theta,z);

                    if (d > dmax) dmax = d;
                }
                //printf("theta_max = %f,theta = %f\n",theta_max,theta);
                return std::max(dmax + density_base + Branch::randt()*density_noise,0.0);
            }


            double density_simple(double r, double theta, double z) {

                double R = radius*b.radius(r,theta,z);
                /*double Rc = radius*b.contour_cercle(r,theta,z);
                    if (r > R*(0.5-eps)) {
                    return std::max(density_out,0.0);
                    }*/

                double dmax = 0;
                for (size_t i=0;i<branches.size();i++) {
                    double d = branches[i].density(r,theta,z);

                    if (d > dmax) dmax = d;
                }
                return std::max(dmax + density_base + Branch::randt()*density_noise,0.0);
            }

            double surface(double theta, double z) {

                double rmax=0;
                for (size_t i=0;i<branches.size();i++) {
                    double r = branches[i].surface(radius,theta,z);
                    if (r > rmax) rmax = r;
                }

                    //double R = radius*b.radius(radius,theta,z);
                //double Rc = radius*b.contour_cercle(r,theta,z);
                return rmax + Branch::randt()*surface_noise;

            }
    };
};


#endif // WOODSEER_LOG_H
