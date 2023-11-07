#include "linkedCell.h"

void List::compForceLC(Cell *grid)
{
    int ic[2], kc[2], kcp[2];
    for (ic[0] = 0; ic[0] < nc[0]; ic[0]++)
        for (ic[1] = 0; ic[1] < nc[1]; ic[1]++)
            for (ParticleList *i = grid[index(ic, nc)]; NULL != i; i = i->next)
            {
                //std::cout<<"compForces i "<<time<<" "<<i->p.getIndex()<<" "<<ic[0]<<" "<<ic[1]<<" "<<nc[0]<<" "<<nc[1]<<std::endl;
                if (!(i->p.isFixed()))
                {

                    i->p.setContactF(Vec2D(0., 0.));
                    i->p.setSocialF(Vec2D(0., 0.));
                    i->p.setDrivingForcePart(Vec2D(0., 0.));
                    i->p.setTorque(0);
                    bool outOfBounding;
                    for (kc[0] = ic[0] - 1; kc[0] <= ic[0] + 1; kc[0]++)
                        for (kc[1] = ic[1] - 1; kc[1] <= ic[1] + 1; kc[1]++)
                        {
                            outOfBounding = false;
                            
                            for (int d = 0; d < 2; d++)
                            {   
                                kcp[d] = kc[d];
                                if (kc[d] < 0)
                                {
                                    if (PBC[d] == true)
                                        kcp[d] = nc[d] + kc[d];
                                    else
                                    {
                                        outOfBounding = true;
                                    }
                                }
                                else if (kc[d] >= nc[d])
                                {
                                    if (PBC[d] == true)
                                        kcp[d] = kc[d] - nc[d];
                                    else
                                    {
                                        outOfBounding = true;
                                    }
                                }
                            }
                             
                            //if ((Vec2D::getDistanceSquared(i->p.getPosition() - minDom, Vec2D(kc[0] * l[0], kc[1] * l[1] )) < r_cut * r_cut) && !outOfBounding)
                            //std::cout <<outOfBounding<<std::endl;
                            if ( !outOfBounding)     
                                for (ParticleList *j = grid[index(kcp, nc)]; NULL != j; j = j->next)
                                {
                                    //if (nStep % frequence_print == 0)
                                    //std::cout <<"distance: "<<i->p.getIndex()<<" "<<j->p.getIndex()<<" "<<Vec2D::getDistanceSquared(i->p.getPosition() - minDom, Vec2D(ic[0] * l[0], ic[1] * l[1]))<<" "<<r_cut * r_cut <<" "<<Vec2D::getDistanceSquared(j->p.getPosition() - minDom, Vec2D(kc[0] * l[0], kc[1] * l[1]))<<std::endl;    
                                    if (i != j)
                                    {
                                        ///std::cout<<"compForces j "<<i->p.getIndex()<<" "<<ic[0]<<" "<<ic[1]<<" "<<nc[0]<<" "<<nc[1]<<std::endl;
                                        //std::cout<<"compForces j "<<i->p.getIndex()<<" "<<j->p.getIndex()<<std::endl;
                                        // std::cout<<"compForces j "<<j->p.getIndex()<<" "<<ic[0]<<" "<<ic[1]<<std::endl;
                                        Vec2D distance = i->p.getPosition() - j->p.getPosition();
                                        for (int d = 0; d < 2; d++)
                                        {
                                            if (PBC[d] == true)
                                            {
                                                if (distance.getComponent(d) > pbcLength[d] * 0.5)
                                                {
                                                    distance.setComponent(d, distance.getComponent(d) - pbcLength[d]);
                                                }
                                                if (distance.getComponent(d) < -pbcLength[d] * 0.5)
                                                {
                                                    distance.setComponent(d, distance.getComponent(d) + pbcLength[d]);
                                                }
                                            }
                                        }
                                        double distanceSqr = distance.getLengthSquared();
                                        //std::cout<<"compForces "<<distanceSqr<<std::endl;
                                        if (distanceSqr < r_cut * r_cut)
                                        {
                                            // std::cout << "long range\n";
                                            Vec2D socialForce;
                                            force.socialForceCS(&i->p, &j->p, &socialForce, distance);
                                            if (socialForce.getLengthSquared() > 0)
                                            {
                                                socialForceHandler sf;
                                                sf.indexI = i->p.getIndex();
                                                sf.indexJ = j->p.getIndex();
                                                sf.socialForce = socialForce;
                                                socialList.push_back(sf);
                                                // i->p.setSocialF(socialForce);
                                            }
                                            // computeSocialForce();
                                            //std::cout<<"compForces "<<distanceSqr<<" "<<long_axis * long_axis<<std::endl;
                                            if (distanceSqr < long_axis * long_axis)
                                            {
                                                // std::cout << "contact range\n";
                                                Vec2D contactForce, contactPoint;
                                                Vec2D shearForce;
                                                force.contactForceHellbing(&i->p, &j->p, &contactPoint, &contactForce, &shearForce, distance);
                                                if (contactForce.getLengthSquared() > 0)
                                                {
                                                    contactForceHandler cf;
                                                    cf.indexI = i->p.getIndex();
                                                    cf.indexJ = j->p.getIndex();
                                                    cf.contactForce = contactForce;
						    cf.shearForce = shearForce;
                                                    cf.contactPoint = contactPoint;
                                                    contactList.push_back(cf);
                                                    // i->p.setContactF(contactForce);
                                                }
                                            }
                                            

                                            // computeContactForce();
                                        }
                                    }
                                }
                             
                       
                        }
                    // driving force
                    /*if (nStep % frequence_print == 0){
                        std::cout << time<< " "<< socialList.size()<<std::endl;
                        for(auto &sf: socialList){
                        std::cout << sf.indexI <<" "<< sf.indexJ<<" "<< sf.socialForce <<std::endl;
                        }
                    }
                    */

                    force.drivingForce(&i->p);
                    
                    Vec2D socialForce;
                    // wall interaction
                    for (auto w : walls)
                    {
                        Vec2D distance = i->p.getPosition() - w.getNearestPoint(i->p.getPosition());
                        double distanceSqr = distance.getLengthSquared();
                        if (distanceSqr < r_cut * r_cut)
                        {

                            socialForceHandler sf;
                            force.socialForceCS(&i->p, &w, &socialForce);
                            /*if (socialForce.getLengthSquared() > 0)
                            {
                                sf.indexI = i->p.getIndex();
                                sf.indexJ = w.getIndex();
                                sf.socialForce = socialForce;
                                socialList.push_back(sf);
                            }*/
                            if (distanceSqr < long_axis * long_axis)
                            {
                                // std::cout << "contact range\n";
                                Vec2D contactForce, contactPoint;
                                force.contactForceHellbing(&i->p, &w, &contactPoint, &contactForce);
                                /*if (contactForce.getLengthSquared() > 0)
                                {
                                    contactForceHandler cf;
                                    cf.indexI = i->p.getIndex();
                                    cf.indexJ = w.getIndex();
                                    cf.contactForce = contactForce;
                                    cf.contactPoint = contactPoint;
                                    contactList.push_back(cf);
                                }*/
                            }
                        }
                    }
                }
            }
}

void List::moveParticlesLC(Cell *grid)
{
    int ic[2], kc[2];
    for (ic[0] = 0; ic[0] < nc[0]; ic[0]++)
        for (ic[1] = 0; ic[1] < nc[1]; ic[1]++)
        {

            ParticleList **q = &grid[index(ic, nc)]; // pointer to predecessor
            ParticleList *i = *q;
            while (NULL != i)
            {   

                if (i->p.getPosition().getComponent(0) >= 0)
                {

		    if(!(i->p.isDomainSel())){
                    	Nout++;
			double rand_num = ((double) rand() / (RAND_MAX));
			double radius_i = i->p.getRadius();
                    	i->p.setObjective(Vec2D(maxDom.getComponent(0), rand_num*(maxDom.Y - radius_i - (minDom.Y + radius_i)) + minDom.Y + radius_i));
			i->p.setDomainSel(true);
		    }
                    /*Particle p; //particle
                    double radius_i = i->p.getRadius()*2;
                    int insert_tries = 0;
                    Vec2D pos;
                    while((insert_tries < 2000)){
                        double rand_num = ((double) rand() / (RAND_MAX));
                        pos = Vec2D(minDom.X + radius_i,  rand_num*(maxDom.Y - radius_i - (minDom.Y + radius_i)) + minDom.Y + radius_i);
                        i->p.setPosition(pos);

                        //i->p.setContactF(Vec2D(0., 0.));
                        //i->p.setSocialF(Vec2D(0., 0.));
                        //i->p.setDrivingForcePart(Vec2D(0., 0.));
                        //i->p.setTorque(0);

                        if(checkIsInContact(&p)){
                            insert_tries++;  
                        }else{
                            insert_tries = 2000;
                        }
                    }
                    
                    i->p.setPosition(pos);*/
                }
                
                for (int d = 0; d < 2; d++)
                {
                    if (PBC[d] == true)
                    {
                        if (i->p.getPosition().getComponent(d) < minDom.getComponent(d))
                        {
                            i->p.setPositionComponent(d, i->p.getPosition().getComponent(d) + pbcLength[d]);
                        }
                        if (i->p.getPosition().getComponent(d) > maxDom.getComponent(d))
                        {
                            i->p.setPositionComponent(d, i->p.getPosition().getComponent(d) - pbcLength[d] + i->p.getRadius());
			    if((i->p.isDomainSel())){
				i->p.setObjective(Vec2D(0., 0.));
				i->p.setDomainSel(false);
			    }
                        }
                    }
                }
                //std::cout<< i->p.getIndex() <<" "<<i->p.getContactF()<<std::endl;
                

                if ((i->p.getPosition().getComponent(0) >= minDom.getComponent(0)) && 
                (i->p.getPosition().getComponent(1) >= minDom.getComponent(1)) && 
                (i->p.getPosition().getComponent(0) <= maxDom.getComponent(0)) && 
                (i->p.getPosition().getComponent(1) <= maxDom.getComponent(1)))
                {
                    for (int d = 0; d < 2; d++)
                        kc[d] = (int)floor((i->p.getPosition().getComponent(d) - minDom.getComponent(d)) / l[d]);
                    //std::cout << "move " << i->p.getIndex() << " " << time << " " << kc[0] << " " << kc[1] <<" " <<ic[0] << " " << ic[1] << std::endl;
                    
                      
                    if ((ic[0] != kc[0]) || (ic[1] != kc[1]))
                    {
                        deleteList(q);
                        insertList(&grid[index(kc, nc)], i);
                        //q = &i->next;
                        //i = *q;
                        //std::cout << "noooooo " << std::endl;
                        
                    }
                    else
                        q = &i->next;
                    
                   
                }
                else
                {
                   deleteList(q);
                   //q = &i->next;
                   //free(i);
                   NP--;
                   std::cout << "nooooo " <<std::endl;
                   //q = &i->next;
                   //i = *q;
                }
                i = *q;
                
                
            }
        }
}

void List::initLC(Cell *grid)
{
    int count = 0;
    for (auto &par : particles)
    {
        int kc[2];
        ParticleList *i;
        i = new ParticleList;
        //i = (ParticleList *)malloc(sizeof(*i));
        i->p = par;
        //std::cout<<i->p.getPosition()<<" "<<i->p.getRadius()<<std::endl;
        if ((i->p.getPosition().getComponent(0) > minDom.getComponent(0)) && 
                (i->p.getPosition().getComponent(1) > minDom.getComponent(1)) && 
                (i->p.getPosition().getComponent(0) < maxDom.getComponent(0)) && 
                (i->p.getPosition().getComponent(1) < maxDom.getComponent(1)))
        {
            for (int d = 0; d < 2; d++)
                kc[d] = (int)floor((i->p.getPosition().getComponent(d) - minDom.getComponent(d)) / l[d]);
            insertList(&grid[index(kc, nc)], i);
            //std::cout<<count<<std::endl;
            count++;
        }
        else{
            NP--;
        }
    }
}

bool List::checkIsInContact(Particle *part)
{

    Vec2D positionPart = part->getPosition();
    double radiusPart = part->getRadius();
    for (auto p : particles)
    {
        Vec2D distance = positionPart - p.getPosition();
        double Rij = radiusPart + p.getRadius();
        double distanceLength = distance.getLength();
        if (distanceLength < Rij)
            return true;
    }
    return false;
}



void List::timeIntegrationLC()
{
    fileSF.open(nameFileSF);
    fileCF.open(nameFileCF);
    filePos.open(nameFilePos);
    fileOut.open(nameFileOut);
    setupInitialConditions();
    initLC(grid_);
    // outputTerminal();
    
    compForceLC(grid_);
    outputResultsLC();
    outputSocialForce();
    outputContactForce();
    clearForceHandlers();
    while ((time < t_end) && (NP > 0))
    {
        actionsBeforeTimeStep();
        time += delta_t;

        compXLC(grid_, delta_t);

        compForceLC(grid_);
        compVLC(grid_, delta_t);
        actionsAfterTimeStep();
        if (nStep % frequence_print == 0)
        {
            outputResultsLC();
            outputSocialForce();
            outputContactForce();
            fileOut << time << " "<< Nout<<std::endl;
            std::cout << time << std::endl;
            
        }
        nStep++;
        
        clearForceHandlers();
        // compoutStatistic_basis(p, N, t);
        // outputResults_basis(p, N, t);
       
    }

    fileSF.close();
    fileCF.close();
    filePos.close();
    fileOut.close();
}
