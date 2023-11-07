#include "Particle.h"
#include "wall.h"
#include "linkedCell.h"
#include "helpers.h"

class Bottleneck: public List
{
public:

    void setupInitialConditions() override
    {
        std::ofstream config;
        config.open("configuration.txt");
        int N_part = 300;// 300; //number of particles
        double r_max = 0.275, r_min = 0.225; //maximum and minimum radius of pedestrians expressed in meters
        double mass = 80; // mass of pedestrians in Kg, in this case, all the pedestrians have the same mass
        double xMax = 0, xMin = -15;
        double yMax = 7.5, yMin = -7.5;
	double xMax2 = 7.5;  
        //double yMax = 0.8, yMin = -0.8;
        //setup(Vec2D(xMin, yMin), Vec2D(xMax, yMax), Vec2D(0.08*5 + 0.5,0.08*5 + 0.5), false, false); // 15X15 hall, in next steps, it is necessary to locate walls whith an entrance
        setup(Vec2D(xMin, yMin), Vec2D(xMax2, yMax),Vec2D((0.08*5 + 0.7),(0.08*5 + 0.7)), true, false); // 15X15 hall, in next steps, it is necessary to locate walls whith an entrance 
        
        config<<Vec2D(xMin, yMin)<<" "<<Vec2D(xMax2, yMax)<<std::endl;
        Particle p; //particle
        p.setMass(80);
        p.setDesiredVel(desiredVelocity);
        p.setObjective(Vec2D(0., 0.));


        //generate initial  configuration
        int particles_inserted = 0;
        for(int i = 0; i < N_part; i++)
        {
            double radius = random_local(r_min, r_max); // uniform distributed radius
            //std::cout<<radius<<std::endl;
            //radius = r_min;
            p.setRadius(radius);
            p.setInertia(p.getMass()*radius*radius*radius*radius/2.);
            Vec2D pos = random2D(Vec2D(xMin + radius, yMin + radius), Vec2D(xMax - radius, yMax - radius)); //select random position
            //pos  = Vec2D(xMin + (2*i + 1)*radius, 0.);
            p.setPosition(pos);

            int insert_tries = 0;
            while((insert_tries < 2000)){
                if(checkIsInContact(&p)){
                    insert_tries++;
                    pos = random2D(Vec2D(xMin + radius, yMin + radius), Vec2D(xMax - radius, yMax - radius));
                    p.setPosition(pos);
                }else{
                     insert_tries = 2000;
                     particles_inserted++;
                     //std::cout<<particles_inserted<<" "<<p.getPosition()<<" "<< p.getRadius()<<std::endl;
                     addParticle(&p);
                }
            }    
        }

        std::cout<<"particles inserted: "<<particles_inserted<<std::endl;
        
        Wall upWall(xMin, yMax, xMax2, yMax);
        upWall.setNormal(Vec2D(0, -1));
        addWall(&upWall);

        config<<xMin<<" "<<yMax<<" "<<xMax<<" "<<yMax<<std::endl;

        Wall downWall(xMin, yMin, xMax2, yMin);
        downWall.setNormal(Vec2D(0, 1));
        addWall(&downWall);

        config<<xMin<<" "<<yMin<<" "<<xMax<<" "<<yMin<<std::endl;

        Wall leftWall(xMin, yMin, xMin, yMax);
        leftWall.setNormal(Vec2D(1, 0));
        addWall(&leftWall);

        config<<xMin<<" "<<yMin<<" "<<xMin<<" "<<yMax<<std::endl;

        Wall orificeTopWall(xMax, orificeD/2., xMax, yMax);
        orificeTopWall.setNormal(Vec2D(-1, 0));
        addWall(&orificeTopWall);

        config<<xMax<<" "<< orificeD/2.<<" "<<xMax<<" "<<yMax<<std::endl;

        Wall orificeBottomWall(xMax, yMin, xMax, -orificeD/2.);
        orificeBottomWall.setNormal(Vec2D(-1, 0));
        addWall(&orificeBottomWall);

        config<<xMax<<" "<<yMin<<" "<<xMax<<" "<<-orificeD/2.<<std::endl;
        
    
        config.close();
        
    }

    void setOrificeD(double d)
    {
        orificeD = d;
    }

    void setDesiredVelocity(double dV)
    {
        desiredVelocity = dV;
    }

private:
    double orificeD, desiredVelocity;

};




int main()
{   
    Bottleneck problem;

    double desiredVelocity, orifice;
    std::ifstream input;
    input.open("input_parameters.in");
    input>> desiredVelocity >> orifice;
    input.close();
    problem.setOrificeD(0.5*orifice);
    problem.setDesiredVelocity(desiredVelocity);
    problem.setTimeMax(1000);
    problem.setTimeStep(1.e-3);
    problem.setFreqPrint(100);
    problem.setSocialParams(2000.,0.08);
    //problem.setSocialParams(0.,0.08);
    problem.setDrivingDT(0.5);
    problem.setContactParams(1.2e5, 2.4e5);
    //problem.setContactParams(2.2e5, 0.0);
    problem.timeIntegrationLC();
    return 0;
}
