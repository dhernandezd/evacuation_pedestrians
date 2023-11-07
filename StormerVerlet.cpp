#include "StormerVerlet.h"

void updateX(Particle *p, double delta_t)
{
    double a = delta_t * .5 / p->getMass();
    double b = delta_t * .5 / p->getInertia();
    //std::cout<<p->getMass()<<" integrator "<<p->getIndex()<<std::endl;
    Vec2D F;
    F = p->getContactF() + p->getSocialF() + p->getDrivingF(); 
    p->setPosition(p->getPosition() + delta_t * (p->getVelocity() + a * F));
    p->setOldF(F);
    p->setAngularPos(p->getAngularPos() + delta_t * (p->getAngularVel() + b * p->getTorque()));
    p->setOldTorque(p->getTorque());
}

void updateV(Particle *p, double delta_t)
{
    Vec2D F;
    F = p->getContactF() + p->getSocialF() + p->getDrivingF();
    double a = delta_t * .5 / p->getMass();
    double b = delta_t * .5 / p->getInertia();
    p->setVelocity(p->getVelocity() + a * (F + p->getOldF()));
    p->setAngularVel(p->getAngularVel() + b * (p->getTorque() + p->getOldTorque()));
}
