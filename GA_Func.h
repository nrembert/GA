#include <iostream>
#include <stdio.h>
#include <string>
#include <cstring>
#include <vector> 
#include <array>
#include <chrono>
#include <ctime>
#include <ratio>
#include <pthread.h>
#include <signal.h>
#include <mutex>
#include <stack>
#include <condition_variable>
#include <unistd.h>
#include <semaphore.h>
#include <thread>
#include <random>
#include <numeric>
#include <math.h>
#include <fstream>
//#include <boost/numeric/odeint.hpp>

class satellite{
    public:
        double inc_deg;
        double ecc;
        double a_km;
        double O_deg;
        double w_deg;
        double TA_deg;
        double pos[3];
        double vel[3];
};

class constellation{
    public:
        int Planes;
        int TotalSats;
        int SatsPerPlane;
        std::vector<satellite> satellites;
};

class Individual 
{ 
public: 
    std::vector<double> chromosome; 
    double fitness;
}; 
std::vector<Individual> population; 
std::vector<Individual> new_generation;

class Access{
    public:
        int Revisit;
        std::vector<int> AccessTime;
};

class ScenarioParameters{
    public:
        std::string epoch;
        int Timestep;
        double SensorHalfAngle;
        double Length;
        std::vector<double> GroundTargets_X;
        std::vector<double> GroundTargets_Y;
};

double fitness_function(std::vector<double> chromosome){
    double x = 0;
    int len = chromosome.size();
    for (int i = 0; i < len; i++) {
        x = x +chromosome[i];
    }
    return x;
}

double random_num(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

const double mu = 3.98618e14; // [m3/s2]


satellite eph2RV (satellite SeedSat){
    double a_m = SeedSat.a_km*1e3;
    double ecc = SeedSat.ecc;
    double RAAN_rad = SeedSat.O_deg*M_PI/180;
    double w_rad = SeedSat.w_deg*M_PI/180;
    double i_rad = SeedSat.inc_deg*M_PI/180;
    double TA_rad = SeedSat.TA_deg*M_PI/180;
    double E_rad = atan(sqrt(1-pow(SeedSat.ecc,2))*sin(TA_rad)/(SeedSat.ecc+cos(TA_rad)));
    double v = TA_rad;
    
    // Radius
    double r = a_m*(1-pow(ecc,2))/(1+ecc*cos(v));
    //std::cout << "R: " << r << std::endl;
    //Angular Momentum
    double h = pow((mu*a_m*(1-pow(ecc,2))),0.5);
    double p = a_m*(1-pow(ecc,2));
    //std::cout << "h: " << h << " p: " << p << std::endl;
    //Position
    SeedSat.pos[0] = r*(cos(RAAN_rad)*cos(w_rad+v) - sin(RAAN_rad)*sin(w_rad+v)*cos(i_rad));
    SeedSat.pos[1] = r*(sin(RAAN_rad)*cos(w_rad+v) + cos(RAAN_rad)*sin(w_rad+v)*cos(i_rad));
    SeedSat.pos[2] = r*(sin(i_rad)*sin(w_rad+v));
    
    SeedSat.vel[0] = SeedSat.pos[0]*h*ecc/(r*p)*sin(v)-(h/r)*(cos(RAAN_rad)*sin(w_rad+v)+sin(RAAN_rad)*cos(w_rad+v)*cos(i_rad));
    SeedSat.vel[1] = SeedSat.pos[1]*h*ecc/(r*p)*sin(v)-(h/r)*(sin(RAAN_rad)*sin(w_rad+v)-cos(RAAN_rad)*cos(w_rad+v)*cos(i_rad));
    SeedSat.vel[2] = SeedSat.pos[2]*h*ecc/(r*p)*sin(v) + h/r*sin(i_rad)*cos(w_rad+v);
    for(int i = 0; i < 3; i++){
        SeedSat.pos[i] *= 1e-3;
        SeedSat.vel[i] *= 1e-3;
    }
    //std::cout << "Pos: " << SeedSat.pos[0] << " " << SeedSat.pos[1] << " " << SeedSat.pos[2] << std::endl;
    //std::cout << "Vel: " << SeedSat.vel[0] << " " << SeedSat.vel[1] << " " << SeedSat.vel[2] << std::endl;
    return SeedSat;

}

satellite CreateSat ( double inc_deg, double ecc, double a_km, double O_deg, double w_deg, double TA_deg){
    satellite Sat;
    Sat.inc_deg = inc_deg;
    Sat.ecc = ecc;
    Sat.a_km = a_km;
    Sat.O_deg = O_deg;
    Sat.w_deg = w_deg;
    Sat.TA_deg = TA_deg;
    Sat = eph2RV (Sat);
    return Sat;
}


constellation CreateConstellation (satellite SeedSat, int SatsPerPlane, int Planes, bool InterPlaneSapcing){
    constellation SatConst;
        double TAdiff = 360/SatsPerPlane;
        double RAANdiff = 360/Planes;
        int TotalSats = SatsPerPlane*Planes;
        int SatCount = 1;
        SatConst.satellites.push_back(SeedSat);
        for (int i = 0; i<Planes; i++){
            for (int j = 0; j<SatsPerPlane; j++){
                if (i==0 && j==0){
                    continue;
                }
                satellite NewSat;
                NewSat = SeedSat;
                //Change the RAAN of the satellite
                NewSat.O_deg += (i)*RAANdiff;
                while (NewSat.O_deg > 360 || NewSat.O_deg < 0){
                    if (NewSat.O_deg > 360){
                        NewSat.O_deg -= 360;
                    }
                    else{
                        NewSat.O_deg += 360;
                    }
                }
                // Change the TA of the satellite
                NewSat.TA_deg +=(j)*TAdiff;
                while (NewSat.TA_deg > 360 || NewSat.TA_deg < 0){
                    if (NewSat.TA_deg > 360){
                        NewSat.TA_deg -= 360;
                    }
                    else{
                        NewSat.TA_deg += 360;
                    }
                } 
                
                // Add this to the constellation
                
                if (!InterPlaneSapcing){
                    NewSat = eph2RV (NewSat);
                }
                //std::cout << "Added a satellite to constellation \n";
                SatConst.satellites.push_back(NewSat);
                SatCount += 1;
            }
            if (InterPlaneSapcing){
                //std::cout << "In Here\n";
                for (int k = (SatCount - SatsPerPlane); k < SatCount; k++){
                    std::cout << k << std::endl;
                        double planespacediff = (i)*360/TotalSats;
                        SatConst.satellites[k].TA_deg += planespacediff;

                        while (SatConst.satellites[k].TA_deg > 360 || SatConst.satellites[k].TA_deg < 0){
                            if (SatConst.satellites[k].TA_deg > 360){
                                SatConst.satellites[k].TA_deg -= 360;
                            }
                            else{
                                SatConst.satellites[k].TA_deg += 360;
                            }
                        }
                        SatConst.satellites[k] = eph2RV (SatConst.satellites[k]);
                }

            }

        }

return SatConst;
}


Access Satellite_Access (constellation SatConst, ScenarioParameters Scenario){
Access SatAccess;
    const int dt = 1;
    //runge_kutta_dopri5<state_type> stepper;

return SatAccess;
}


