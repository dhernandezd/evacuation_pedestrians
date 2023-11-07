#ifndef PARTICLE_H
#define PARTICLE_H

#include "vector2d.h"

class Particle
{
private:
    int index;
    bool fixed;
    bool domainSel;
    double m;     // mass
    double R;     // radius
    double omega; // angular velocity
    double phi;   // angular position
    double desiredVelocity;
    //double angAcc;
    double torque;      // torque
    double oldTorque;
    Vec2D position;     // position
    //Vec2D acceleration; // position
    Vec2D objective;
    Vec2D velocity;     // velocity
    Vec2D contactForce; // granular force
    Vec2D socialForce;  // social force
    Vec2D drivingForcePart;  // social force
    Vec2D F_old;        // old force
    double inertia;

public:
    Particle(){
        index = 0;
        fixed = false;
	domainSel = false;
        m =0;
        R =0;
        omega= 0;
        phi = 0;
        torque = 0;
        oldTorque = 0;
        inertia = 0;
    }
    
    void setPosition(Vec2D positionIn)
    {
        position = positionIn;
    }

    void setPosition(double x, double y)
    {
        position.set(x, y);
    }

    void setPositionComponent(int i,double positionIn)
    {
        position.setComponent(i, positionIn);
    }

    void setRadius(double Radius)
    {
        R = Radius;
    }

    void setVelocity(Vec2D velocityIn)
    {
        velocity = velocityIn;
    }

    void setContactF(Vec2D contactForceIn)
    {
        contactForce = contactForceIn;
    }


    void setSocialF(Vec2D socialForceIn)
    {
        socialForce = socialForceIn;
    }

    void setDrivingForcePart(Vec2D drivingForcePartIn)
    {
        drivingForcePart = drivingForcePartIn;
    }

    void setObjective(Vec2D obj)
    {
        objective = obj;
    }

    void setDesiredVel(double dV)
    {
        desiredVelocity = dV;
    }

    void setMass(double Mass)
    {
        m = Mass;
    }

    void setInertia(double Inertia){
        inertia = Inertia;
    }

    void setTorque(double torqueIn)
    {
        torque = torqueIn;
    }

    void setOldTorque(double torqueIn)
    {
        oldTorque = torqueIn;
    }

    void setAngularPos(double ang)
    {
        phi = ang;
    }

    void setAngularVel(double Omega)
    {
        omega = Omega;
    }

    void setOldF(Vec2D oldF)
    {
        F_old = oldF;
    }

    void setIndex(int I)
    {
        index = I;
    }

    double getMass()
    {
        return m;
    }

    double getRadius()
    {
        return R;
    }

    void setFixed(){
        fixed = true;
    }

    double getInertia()
    {
        return inertia;
    }


    Vec2D getPosition()
    {
        return position;
    }

    Vec2D getVelocity()
    {
        return velocity;
    }


    Vec2D getSocialF()
    {
        return socialForce;
    }

    Vec2D getContactF()
    {
        return contactForce;
    }

    Vec2D getDrivingF()
    {
        return drivingForcePart;
    }

    Vec2D getObjective( )
    {
        return objective;
    }

    double getDesiredVel()
    {
        return desiredVelocity;
    }

    Vec2D getOldF()
    {
        return F_old;
    }

    int getIndex(){
        return index;
    }

    double getAngularPos()
    {
        return phi;
    }

    double getAngularVel()
    {
        return omega;
    }

    double getTorque()
    {
        return torque;
    }

    double getOldTorque()
    {
        return oldTorque;
    }

    bool isFixed(){
        return fixed;
    }
    
    bool isDomainSel(){
	return domainSel;
    }
    
    void setDomainSel(bool domainSel_i){
        domainSel = domainSel_i;
    }


};

#endif
