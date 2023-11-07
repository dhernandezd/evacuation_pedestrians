#include "io.h"

void readPositions(std::string nameFile, std::vector<Particle> *particles){
    double x,y;
    Particle p;
    std::ifstream positionsFile(nameFile);
    int I = 0;
    while(positionsFile >> x >> y){
        p.setPosition(x,y);
        p.setIndex(I);
        particles->push_back(p);
        I++;
    }
    positionsFile.close();

};

void readFixedParticles(std::string nameFile, std::vector<Particle> *particles){
    double x,y;
    Particle p;
    std::ifstream positionsFile(nameFile);
    int I = 0;
    while(positionsFile >> x >> y){
        p.setPosition(x,y);
        p.setIndex(I);
        p.setFixed();
        particles->push_back(p);
        I--;
    }
    positionsFile.close();

};

void readRadii(std::string nameFile, std::vector<Particle> particles){
    double R;
    std::ifstream radiiFile(nameFile);
    for(auto &p: particles){
        radiiFile >> R;
        p.setRadius(R);
    }
    radiiFile.close();
};

void readMasses(std::string nameFile, std::vector<Particle> particles){
    double mass;
    std::ifstream massFile(nameFile);
    for(auto &p: particles){
        massFile >> mass;
        p.setMass(mass);
    }
    massFile.close();
};

void readWalls(std::string nameFile, std::vector<Wall> *walls){
    double x1, x2, y1, y2, n1, n2;
    std::ifstream wallFile(nameFile);
    int I = -1;
    while(wallFile >> x1 >> x2 >> y1 >> y2 >> n1 >> n2)
    {
        Wall w(x1, x2, y1, y2);
        w.setNormal(Vec2D(n1, n2));
        w.setIndex(I);
        walls->push_back(w);
        I--;
    }
    wallFile.close();
};
