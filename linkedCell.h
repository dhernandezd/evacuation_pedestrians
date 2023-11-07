#ifndef LINKEDCELL_H
#define LINKEDCELL_H

//#include "vector2d.h"
#include "forces.h"
#include "wall.h"
#include <iostream>
#include "Particle.h"
#include "StormerVerlet.h"
#include <vector>
#include <fstream>
#include <cstdlib>
//#include <random>
//std::random_device rd_LC;
//std::mt19937 gen_LC(rd_LC());
//std::uniform_real_distribution<> random_gen(0.,1.);
//#include "helpers.h"

#define index(ic, nc) ((ic)[0] + (nc)[0] * (ic)[1]) // macro to iterate the cells

typedef struct ParticleList
{
    Particle p;
    struct ParticleList *next;
} ParticleList;

typedef struct socialForceHandler
{
    Vec2D socialForce;
    int indexI, indexJ;

} socialForceHandler;

typedef struct contactForceHandler
{
    Vec2D contactForce;
    Vec2D shearForce;
    Vec2D contactPoint;
    int indexI, indexJ;

} contactForceHandler;

typedef ParticleList *Cell;

class List
{
private:
    bool PBC[2];
    double l[2];
    double pbcLength[2];
    int nc[2];
    Vec2D minDom;
    Vec2D maxDom;
    Cell *grid_;
    std::vector<contactForceHandler> contactList;
    std::vector<socialForceHandler> socialList;
    std::vector<Wall> walls;
    std::vector<Particle> particles;
    Forces force;
    double time, delta_t, t_end;
    double r_cut, long_axis;
    std::string nameFileSF, nameFileCF, nameFilePos, nameFileOut;
    std::ofstream fileSF, fileCF, filePos, fileOut;
    int frequence_print, nStep;
    int NP, NW; // noumber of particles and walls, respectively
    int Nout;
 

public:
    List(){PBC[0] = false; PBC[1] = false;
    nameFileSF = "social_force.txt"; nameFileCF = "contact_force.txt"; nameFilePos = "data_pos.txt";
    nameFileOut = "flux.txt";
    nStep =1; NP = 0; NW = -1;
    }
    List(Vec2D minDomain, Vec2D maxDomain, Vec2D deltaL, bool pbcX, bool pbcY)
    {
        PBC[0] = pbcX;
        PBC[1] = pbcY;
        l[0] = deltaL.X;
        l[1] = deltaL.Y;
        minDom = minDomain;
        maxDom = maxDomain;
        pbcLength[0] = (maxDomain.X - minDomain.X);
        pbcLength[1] = (maxDomain.Y - minDomain.Y);
        nc[0] = (int)floor(pbcLength[0] / l[0]);
        nc[1] = (int)floor(pbcLength[1] / l[1]);
        int pnc=1;
        for (int d=0; d<2; d++)
            pnc *= nc[d];
        //grid_ = (Cell*)malloc(pnc*sizeof(*grid));
        time = 0.;
        nameFileSF = "social_force.txt"; nameFileCF = "contact_force.txt"; nameFilePos = "data_pos.txt";
        nameFileOut = "flux.txt";
        nStep =1;    
        NP=0;
        NW=-1;
        Nout =0;
    }

    void setup(Vec2D minDomain, Vec2D maxDomain, Vec2D deltaL, bool pbcX, bool pbcY)
    {
        PBC[0] = pbcX;
        PBC[1] = pbcY;
        l[0] = deltaL.X;
        l[1] = deltaL.Y;
        minDom = minDomain;
        maxDom = maxDomain;
        pbcLength[0] = (maxDomain.X - minDomain.X);
        pbcLength[1] = (maxDomain.Y - minDomain.Y);
        nc[0] = (int)ceil(pbcLength[0] / l[0]) ;
        nc[1] = (int)ceil(pbcLength[1] / l[1]) ;
        time = 0.;
        nameFileSF = "social_force.txt"; nameFileCF = "contact_force.txt"; nameFilePos = "data_pos.txt";
        nameFileOut = "flux.txt";
        nStep =1;
        int pnc=1;
        for (int d=0; d<2; d++)
            pnc *= nc[d];
        grid_ = (Cell*)malloc(pnc*sizeof(*grid_));
        long_axis =0.7 + 0.001;
        //r_cut = (long_axis>5*0.08)?long_axis:5*0.08;
        r_cut = 0.6*5 + 0.7 + 0.001;
        NP = 0;
        NW = -1;
        Nout = 0;
    }


    void insertList(ParticleList **root_list, ParticleList *i)
    {
        i->next = *root_list;
        *root_list = i;
    }

    void deleteList(ParticleList **q)
    {
        *q = (*q)->next; // (*q)->next points to element to be removed
    }

    void compForceLC(Cell *grid);

    void moveParticlesLC(Cell *grid);

    void initLC(Cell *grid);

    void freeLC()
    {
        free(grid_);
    }

    void compXLC(Cell *grid, double delta_t)
    {
        int ic[2];
        for (ic[0] = 0; ic[0] < nc[0]; ic[0]++)
            for (ic[1] = 0; ic[1] < nc[1]; ic[1]++)
                for (ParticleList *i = grid[index(ic, nc)]; NULL != i; i = i->next)
                {
                    updateX(&i->p, delta_t);
                    //std::cout<< i->p.getIndex() <<" "<<i->p.getContactF()<<std::endl;
                }
        moveParticlesLC(grid);
    }
    void compVLC(Cell *grid, double delta_t)
    {
        int ic[2];
        for (ic[0] = 0; ic[0] < nc[0]; ic[0]++)
            for (ic[1] = 0; ic[1] < nc[1]; ic[1]++)
                for (ParticleList *i = grid[index(ic, nc)]; NULL != i; i = i->next)
                    updateV(&i->p, delta_t);
    }
    void outputResultsLC(){
        filePos << time << " "<< NP <<std::endl;
        int ic[2];
        for (ic[0] = 0; ic[0] < nc[0]; ic[0]++)
            for (ic[1] = 0; ic[1] < nc[1]; ic[1]++)
                for (ParticleList *i = grid_[index(ic, nc)]; NULL != i; i = i->next){
                    //filePos <<i->p.getIndex()<<" "<<i->p.getPosition() << " "<< i->p.getVelocity()<< 
                    //" "<< i->p.getAngularPos() <<" "<<i->p.getAngularVel() <<" "<<i->p.getRadius()<<std::endl;
                    filePos <<i->p.getIndex()<<" "<<i->p.getPosition() << " "<< i->p.getVelocity()<<" "<<i->p.getRadius()<<std::endl; 
                }
    }

    void outputTerminal(){
        std::cout << time << " "<< NP <<std::endl;
        int ic[2];
        for (ic[0] = 0; ic[0] < nc[0]; ic[0]++)
            for (ic[1] = 0; ic[1] <= nc[1]; ic[1]++)
                for (ParticleList *i = grid_[index(ic, nc)]; NULL != i; i = i->next){
                    std::cout<<i->p.getIndex()<<" "<<i->p.getRadius()<<" "<<i->p.getPosition() << " "<< i->p.getVelocity()<< 
                    " "<< i->p.getAngularPos() <<" "<<i->p.getAngularVel() << std::endl; 
                }
    }

    void outputSocialForce(){
        fileSF << time << " "<< socialList.size()<<std::endl;
        for(auto &sf: socialList){
            fileSF << sf.indexI <<" "<< sf.indexJ<<" "<< sf.socialForce <<std::endl;
        }
    }

    void outputContactForce(){
        fileCF << time<< " "<< contactList.size()<<std::endl;
        for(auto &cf: contactList){
            //fileCF << cf.indexI <<" "<< cf.indexJ<<" "<< cf.contactForce <<" "<<cf.contactPoint<<std::endl;
            fileCF << cf.indexI <<" "<< cf.indexJ<<" "<< cf.contactForce<<" "<<cf.shearForce <<std::endl;
        }
    }

    void timeIntegrationLC();

    virtual void setupInitialConditions(){} 

    virtual void actionsAfterTimeStep(){}

    virtual void actionsBeforeTimeStep(){}

    void clearForceHandlers(){
        contactList.clear();
        socialList.clear();
    }

    void setTime(double timeIn){
        time =  timeIn;
    }    

    void setTimeMax(double timeIn){
        t_end =  timeIn;
    }                                                                                   

    void setTimeStep(double dt){
        delta_t = dt;
    }

    void setFreqPrint(int Nref){
        frequence_print = Nref;
    }

    double getTime(){
        return time;
    }

    int getTimeStep(){
        return delta_t;
    }

    int getNStep(){
        return nStep;
    }

    int getFreqPrint(){
        return frequence_print;
    }

    void addParticle(Particle *part){
        part->setIndex(NP);
        particles.push_back(*part);
        NP++;
    }

    void addWall(Wall *wall){
        wall->setIndex(NW);
        walls.push_back(*wall);
        NW--;
    }

    void setDrivingDT(double DT)
    {
        force.setDrivingDT(DT);
    }

    bool checkIsInContact(Particle *part);

    ~List()
    {
        freeLC(); 
        }

    void setSocialParams(const double Aconst,const double Bconst){
        force.setSocialParams(Aconst,Bconst);
        //r_cut = 5*Bconst + 0.7;
        r_cut = 5*Bconst + 0.5 + 0.001;
        //r_cut = (0.7>5*0.08)?(0.7+ 0.001):5*0.08;
        //r_cut = 0.7;
    }
    void setSocialParams(const double Aconst,const double Bconst, const double inTau){
        force.setSocialParams(Aconst,Bconst, inTau);
        //r_cut = 5*Bconst + 0.7;
        r_cut = 5*Bconst + 0.5 + 0.001;
        //r_cut = (0.7>5*0.08)?(0.7+ 0.001):5*0.08;
        //r_cut = 0.7;
    }
    void setVisualParam(const double inLambda){
    	force.setVisualParam(inLambda);
    }
    void setContactParams(const double stiffnessIn,const double dissipationIn){
        force.setContactParams(stiffnessIn, dissipationIn);
    }
};

#endif
    
