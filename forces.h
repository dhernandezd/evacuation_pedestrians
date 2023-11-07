#ifndef FORCES_H
#define FORCES_H

#include "vector2d.h"
#include "Particle.h"
#include "wall.h"

class Forces
{
private:
    double stiffness, dissipation; // contact parameters
    double A, B, tau;          // social paramters
    double dT; //driving force parameters
    double lambda; //visual effect parameter
public:
    Forces(/* args */);
    ~Forces();
    void drivingForce(Particle *p);

    //![Particle-Particle interactions]
    void contactForceHellbing(Particle *p, Particle *j, Vec2D *contactPoint, Vec2D *contactForce);
    void contactForceHellbing(Particle *p, Particle *j, Vec2D *contactPoint, Vec2D *contactForce, Vec2D *shearForce);
    void contactForceHellbing(Particle *p, Particle *j, Vec2D *contactPoint, Vec2D *contactForce, Vec2D *shearForce, Vec2D distance);
    
    //circular specification
    void socialForceCS(Particle *p, Particle *j, Vec2D *socialForce, Vec2D distance);
    
    // elliptical specification 1(ES1)
    void socialForceES1(Particle *p, Particle *j, Vec2D *socialForce);

    // elliptical specification 2(ES2)
    void socialForceES2(Particle *p, Particle *j, Vec2D *socialForce);
    
    // new elliptical specification (NES) //TODO: verify force
    void socialForceNES(Particle *p, Particle *j, Vec2D *socialForce);

    //![Particle-Wall interactions]
    void contactForceHellbing(Particle *p, Wall *j, Vec2D *contactPoint, Vec2D *contactForce);
    
    //circular specification
    void socialForceCS(Particle *p, Wall *j, Vec2D *socialForce);
    
    // elliptical specification 1(ES1)
    void socialForceES1(Particle *p, Wall *j, Vec2D *socialForce);

    // elliptical specification 2(ES2)
    void socialForceES2(Particle *p, Wall *j, Vec2D *socialForce);
    
    // new elliptical specification (NES) //TODO: verify force
    void socialForceNES(Particle *p, Wall *j, Vec2D *socialForce);


    void visualEffect(Particle *p, Particle *j, Vec2D desiredVel, double *wij);
    void visualEffect(Particle *p, Particle *j, double *wij);

    void visualEffect(Particle *p, Wall *j, Vec2D desiredVel, double *wij);
    void visualEffect(Particle *p, Wall *j, double *wij);

    void setSocialParams(const double Aconst,const double Bconst){
        A = Aconst;
        B = Bconst;
    }
    
    void setContactParams(const double stiffnessIn,const double dissipationIn){
        stiffness = stiffnessIn;
        dissipation = dissipationIn;
    }
    
    void setSocialParams(const double Aconst,const double Bconst,const double tauConst){
        A = Aconst;
        B = Bconst;
        tau = tauConst;
    }

    void setDrivingDT(const double dTin){
        dT = dTin;
    }

    void setVisualParam(const double lambdaIn){
        lambda = lambdaIn;
    }
};


#endif
