#include "forces.h"

void Forces::drivingForce(Particle *p)
{
    Vec2D e_io = p->getObjective() - p->getPosition();
    double length = e_io.getLength();
    if (length != 0)
    {
	double length = e_io.getLength();
        e_io = e_io / length;
        p->setDrivingForcePart(p->getMass() * (p->getDesiredVel() * e_io - p->getVelocity()) / dT);
    }/*else{
	e_io = Vec2D(7.5,0.) - p->getPosition();
	e_io = e_io / length;
	p->setDrivingForcePart(p->getMass() * (p->getDesiredVel() * e_io - p->getVelocity()) / dT);
	}*/
    
}

void Forces::contactForceHellbing(Particle *p, Particle *j, Vec2D *contactPoint, Vec2D *contactForce)
{
    Vec2D distance = p->getPosition() - j->getPosition();
    double Rij = p->getRadius() + j->getRadius();
    double distanceLength = distance.getLength();
    if (distanceLength < Rij)
    {
        Vec2D norm = distance / distanceLength;
        Vec2D tangentialVelRel = (j->getVelocity() - p->getVelocity()) - Vec2D::dot((j->getVelocity() - p->getVelocity()), norm) * norm;
        double overlap = (Rij - distanceLength);
        Vec2D normalForce = overlap * stiffness * norm;
        Vec2D tangentialForce = overlap * dissipation * tangentialVelRel;
        *contactForce = normalForce + tangentialForce;
        p->setContactF(p->getContactF() + *contactForce);
        *contactPoint = (p->getPosition() + p->getRadius() * norm);
        p->setTorque(p->getTorque() + Vec2D::cross(*contactPoint, p->getContactF()));
    }
}

void Forces::contactForceHellbing(Particle *p, Particle *j, Vec2D *contactPoint, Vec2D *contactForce, Vec2D *shearForce, Vec2D distance)
{
    //Vec2D distance = p->getPosition() - j->getPosition();
    double Rij = p->getRadius() + j->getRadius();
    double distanceLength = distance.getLength();
    if (distanceLength < Rij)
    {
        Vec2D norm = distance / distanceLength;
        Vec2D tangentialVelRel = (j->getVelocity() - p->getVelocity()) - Vec2D::dot((j->getVelocity() - p->getVelocity()), norm) * norm;
        double overlap = (Rij - distanceLength);
        Vec2D normalForce = overlap * stiffness * norm;
        Vec2D tangentialForce = overlap * dissipation * tangentialVelRel;
	*shearForce = tangentialForce;
        *contactForce = normalForce + tangentialForce;
        p->setContactF(p->getContactF() + *contactForce);
        *contactPoint = (p->getPosition() + p->getRadius() * norm);
        p->setTorque(p->getTorque() + Vec2D::cross(*contactPoint, p->getContactF()));
    }
}

void Forces::contactForceHellbing(Particle *p, Particle *j, Vec2D *contactPoint, Vec2D *contactForce, Vec2D *shearForce)
{
    Vec2D distance = p->getPosition() - j->getPosition();
    double Rij = p->getRadius() + j->getRadius();
    double distanceLength = distance.getLength();
    if (distanceLength < Rij)
    {
        Vec2D norm = distance / distanceLength;
        Vec2D tangentialVelRel = (j->getVelocity() - p->getVelocity()) - Vec2D::dot((j->getVelocity() - p->getVelocity()), norm) * norm;
        double overlap = (Rij - distanceLength);
        Vec2D normalForce = overlap * stiffness * norm;
        Vec2D tangentialForce = overlap * dissipation * tangentialVelRel;
	*shearForce = tangentialForce;
        *contactForce = normalForce + tangentialForce;
        p->setContactF(p->getContactF() + *contactForce);
        *contactPoint = (p->getPosition() + p->getRadius() * norm);
        p->setTorque(p->getTorque() + Vec2D::cross(*contactPoint, p->getContactF()));
    }
}



void Forces::socialForceCS(Particle *p, Particle *j, Vec2D *socialForce, Vec2D distance)
{
    //Vec2D distance = p->getPosition() - j->getPosition();
    double Rij = p->getRadius() + j->getRadius();
    double distanceLength = distance.getLength();
   // if ((distanceLength < Rij)) distanceLength = Rij; 
    if ((distanceLength - Rij) < 5 * B){
        Vec2D norm = Vec2D::getUnitVector(distance);
        *socialForce = A * exp((Rij - distanceLength) / B) * norm;
        p->setSocialF(p->getSocialF() + *socialForce);
    }
}

void Forces::socialForceES1(Particle *p, Particle *j, Vec2D *socialForce)
{
    Vec2D distance = p->getPosition() - j->getPosition();
    double distanceLength = distance.getLength();
    Vec2D norm = Vec2D::getUnitVector(distance);
    double Rij = p->getRadius() + j->getRadius();
    Vec2D y_ij = distance - j->getVelocity() * tau;
    double auxv_ij = j->getVelocity().getLength() * tau;
    double normAuxDiff = y_ij.getLength();
    double b_ij = sqrt((distanceLength + y_ij.getLength()) * (distanceLength + y_ij.getLength()) - auxv_ij * auxv_ij) / 2.;
    if (normAuxDiff != 0)
    {
        *socialForce = A * exp((-b_ij) / B) * (norm + j->getVelocity() * tau / normAuxDiff) * (distanceLength + normAuxDiff) / (4 * b_ij);
        p->setSocialF(p->getSocialF() + *socialForce);
    }
    else
    {
        *socialForce = A * exp((-b_ij) / B) * (norm) * (distanceLength) / (4 * b_ij);
        p->setSocialF(p->getSocialF() + *socialForce);
    }
}

// elliptical specification 2(ES2)
void Forces::socialForceES2(Particle *p, Particle *j, Vec2D *socialForce)
{
    Vec2D distance = p->getPosition() - j->getPosition();
    Vec2D v_ij = (p->getVelocity() - j->getVelocity());
    double distanceLength = distance.getLength();
    Vec2D norm = Vec2D::getUnitVector(distance);
    double Rij = p->getRadius() + j->getRadius();
    Vec2D y_ij = distance - v_ij * tau;
    double normAuxDiff = y_ij.getLength();
    double auxv_ij = v_ij.getLength() * tau;

    double b_ij = sqrt((distanceLength + y_ij.getLength()) * (distanceLength + y_ij.getLength()) - auxv_ij * auxv_ij) / 2.;
    // if (normAuxDiff > 1e-8)
    *socialForce = A * exp((-b_ij) / B) * (norm + y_ij / normAuxDiff) * (distanceLength + normAuxDiff) / (4 * b_ij);
    p->setSocialF(p->getSocialF() + *socialForce);
    // else
    // return A*exp((-b_ij)/B)*(distance/d_ij)*(d_ij)/(4*b_ij);
}

void Forces::socialForceNES(Particle *p, Particle *j, Vec2D *socialForce)
{
    Vec2D distance = p->getPosition() - j->getPosition();
    Vec2D v_ij = (p->getVelocity() - j->getVelocity());
    double distanceLength = distance.getLength();
    Vec2D norm = Vec2D::getUnitVector(distance);
    double Rij = p->getRadius() + j->getRadius();
    Vec2D y_ij = distance - v_ij * tau;
    double normAuxDiff = y_ij.getLength();
    double auxv_ij = v_ij.getLength() * tau;
    double auxHeadOn = 1. / sqrt(1. + p->getVelocity().getLength());
    double b_ij = sqrt((distanceLength + y_ij.getLength()) * (distanceLength + y_ij.getLength()) - auxv_ij * auxv_ij) / 2. * auxHeadOn;
    if (normAuxDiff != 0)
    {
        *socialForce = A * exp((-b_ij) / B) * (norm + y_ij / normAuxDiff) * (distanceLength + normAuxDiff) / (4 * b_ij) * auxHeadOn;
        p->setSocialF(p->getSocialF() + *socialForce);
    }
    else
    {
        *socialForce = A * exp((-b_ij) / B) * (norm) * (distanceLength) / (4 * b_ij) * auxHeadOn;
        p->setSocialF(p->getSocialF() + *socialForce);
    }
}

void Forces::visualEffect(Particle *p, Particle *j, Vec2D desiredVel, double *wij)
{
    Vec2D distance = p->getPosition() - j->getPosition();
    double cosineAng = Vec2D::dot(desiredVel, -distance) / (desiredVel.getLength() * distance.getLength());
    *wij = lambda + (1 - lambda) * (1 + cosineAng) / 2.;
}

void Forces::contactForceHellbing(Particle *p, Wall *j, Vec2D *contactPoint, Vec2D *contactForce)
{
    Vec2D distance = p->getPosition() - j->getNearestPoint(p->getPosition());
    double Rij = p->getRadius();
    double distanceLength = distance.getLength();
    if (distanceLength < Rij)
    {
        Vec2D norm = j->getNormal();
        Vec2D tangentialVelRel = (p->getVelocity()) - Vec2D::dot((p->getVelocity()), norm) * norm;
        double overlap = (Rij - distanceLength);
        Vec2D normalForce = overlap * stiffness * norm;
        Vec2D tangentialForce = -overlap * dissipation * tangentialVelRel;
        *contactForce = normalForce + tangentialForce;
        p->setContactF(p->getContactF() + *contactForce);
        *contactPoint = (p->getPosition() + p->getRadius() * norm);
        p->setTorque(p->getTorque() + Vec2D::cross(*contactPoint, p->getContactF()));
    }
}

void Forces::socialForceCS(Particle *p, Wall *j, Vec2D *socialForce)
{
    Vec2D distance = p->getPosition() - j->getNearestPoint(p->getPosition());
    double Rij = p->getRadius();
    double distanceLength = distance.getLength();
    //if ((distanceLength <= Rij)) distanceLength = Rij;
    if ((distanceLength - Rij) < 5 * B)
    {
        Vec2D norm =  Vec2D::getUnitVector(distance);
        *socialForce = A * exp((Rij - distanceLength) / B) * norm;
        p->setSocialF(p->getSocialF() + *socialForce);
    }
}

void Forces::socialForceES1(Particle *p, Wall *j, Vec2D *socialForce)
{
    socialForceCS(p, j, socialForce);
}

void Forces::socialForceES2(Particle *p, Wall *j, Vec2D *socialForce)
{
    Vec2D distance = p->getPosition() - j->getNearestPoint(p->getPosition());
    double Rij = p->getRadius();
    double distanceLength = distance.getLength();
    Vec2D v_ij = (p->getVelocity());
    Vec2D norm = Vec2D::getUnitVector(distance);
    Vec2D y_ij = distance - v_ij * tau;
    double normAuxDiff = y_ij.getLength();
    double auxv_ij = v_ij.getLength() * tau;

    double b_ij = sqrt((distanceLength + y_ij.getLength()) * (distanceLength + y_ij.getLength()) - auxv_ij * auxv_ij) / 2.;
    // if (normAuxDiff > 1e-8)
    *socialForce = A * exp((-b_ij) / B) * (norm + y_ij / normAuxDiff) * (distanceLength + normAuxDiff) / (4 * b_ij);
    p->setSocialF(p->getSocialF() + *socialForce);
    // else
    // return A*exp((-b_ij)/B)*(distance/d_ij)*(d_ij)/(4*b_ij);
}

void Forces::socialForceNES(Particle *p, Wall *j, Vec2D *socialForce)
{
    Vec2D distance = p->getPosition() - j->getNearestPoint(p->getPosition());
    Vec2D v_ij = (p->getVelocity());
    double distanceLength = distance.getLength();
    Vec2D norm = Vec2D::getUnitVector(distance);
    double Rij = p->getRadius();
    Vec2D y_ij = distance - v_ij * tau;
    double normAuxDiff = y_ij.getLength();
    double auxv_ij = v_ij.getLength() * tau;
    double auxHeadOn = 1. / sqrt(1. + p->getVelocity().getLength());
    double b_ij = sqrt((distanceLength + y_ij.getLength()) * (distanceLength + y_ij.getLength()) - auxv_ij * auxv_ij) / 2. * auxHeadOn;
    if (normAuxDiff != 0)
    {
        *socialForce = A * exp((-b_ij) / B) * (norm + y_ij / normAuxDiff) * (distanceLength + normAuxDiff) / (4 * b_ij) * auxHeadOn;
        p->setSocialF(p->getSocialF() + *socialForce);
    }
    else
    {
        *socialForce = A * exp((-b_ij) / B) * (norm) * (distanceLength) / (4 * b_ij) * auxHeadOn;
        p->setSocialF(p->getSocialF() + *socialForce);
    }
}

void Forces::visualEffect(Particle *p, Wall *j, Vec2D desiredVel, double *wij)
{
    Vec2D distance = p->getPosition() - j->getNearestPoint(p->getPosition());
    double cosineAng = Vec2D::dot(desiredVel, -distance) / (desiredVel.getLength() * distance.getLength());
    *wij = (lambda + (1 - lambda) * (1 + cosineAng) / 2.);
}

Forces::Forces()
{
}

Forces::~Forces()
{
}
