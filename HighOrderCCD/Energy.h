#ifndef ENERGY_H
#define ENERGY_H

#include "Utils/CCDUtils.h"

#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"

#include <vector>

PRJ_BEGIN

class Energy
{
  public:
    typedef Eigen::MatrixXd Data;
    typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
    typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;
    

    static double plane_whole_energy(const Data& spline, const double& piece_time,
                                     const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                     const std::vector<std::vector<double>>& d_lists)
    { 
        //return dirichlet_energy(spline)+lambda*barrier_energy(spline,V,F,bvh);
        return ks*dynamic_energy(spline,piece_time)+
               lambda*plane_barrier_energy(spline,c_lists, d_lists)+ lambda*bound_energy(spline,piece_time)+
               kt*whole_weight*piece_time;
    }

    static double dynamic_energy(const Data& spline, const double& piece_time)
    {
        double energy=0;

        //dynamic energy
        Eigen::MatrixXd M;
        for(int sp_id=0;sp_id<piece_num;sp_id++)
        {
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
            
            M=convert_list[sp_id].transpose()*M_dynamic*convert_list[sp_id];
            for(int j=0;j<3;j++)
            {
                Eigen::VectorXd x=bz.col(j);
                energy+=1/std::pow(time_weight[sp_id]*piece_time,2*der_num-1)*0.5*x.transpose()*M*x;
            }
        }
        return energy;
    }

    static double plane_barrier_energy(const Data& spline, 
                                       const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                       const std::vector<std::vector<double>>& d_lists)
    {
        double energy=0;
       
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {        
            std::vector<Eigen::Vector3d> c_list=c_lists[tr_id];
            std::vector<double> d_list=d_lists[tr_id];
            if(c_list.size()==0)
                continue;

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        
            std::vector<Eigen::RowVector3d> P(order_num+1);
            for(int j=0;j<=order_num;j++)
            {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                    P[j]+=basis(j,j0)*bz.row(j0);
                } 
                
            }

            double d;

            for(unsigned int k=0;k<c_list.size();k++)
            {
                for(int j=0;j<=order_num;j++)
                {
                    d=P[j].dot(c_list[k])+d_list[k];

                    if(d<=0)
                        return INFINITY;
                    
                    if(d<margin)
                    { 
                       energy+=-weight*(d-margin)*(d-margin)*log(d/margin); 
                       //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
                    }
                }

            }
            
           
        }

        return energy;  
    }

    static double bound_energy(const Data& spline,const double& piece_time)
    {
        double energy=0;

        for(unsigned int tr_id=0;tr_id<vel_tree.size();tr_id++)
        {
        
            int sp_id=std::get<0>(vel_tree[tr_id]);
            double weight=std::get<1>(vel_tree[tr_id]).second-std::get<1>(vel_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(vel_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);

            std::vector<Eigen::RowVector3d> P(order_num+1);
            for(int j=0;j<=order_num;j++)
            {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                  P[j]+=basis(j,j0)*bz.row(j0);
                } 
            }
                
            double d;
           
            for(int j=0;j<order_num;j++)
            {
                Eigen::RowVector3d vel=order_num*(P[j+1]-P[j]);
                //d=vel_limit*piece_time-vel.norm()/weight;
                d=vel_limit-vel.norm()/(weight*time_weight[sp_id]*piece_time);
                if(d<=0)
                    return INFINITY;
                
                if(d<margin)
                { 
                   energy+=-weight*(d-margin)*(d-margin)*log(d/margin); 
                }
            }
        }


        

        for(unsigned int tr_id=0;tr_id<acc_tree.size();tr_id++)
        {
        
            int sp_id=std::get<0>(acc_tree[tr_id]);
            double weight=std::get<1>(acc_tree[tr_id]).second-std::get<1>(acc_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(acc_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);

            std::vector<Eigen::RowVector3d> P(order_num+1);
            for(int j=0;j<=order_num;j++)
            {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                  P[j]+=basis(j,j0)*bz.row(j0);
                } 
            }
                
            double d;
           
            for(int j=0;j<order_num-1;j++)
            {
                Eigen::RowVector3d acc=order_num*(order_num-1)*(P[j+2]-2*P[j+1]+P[j]);
                //d=acc_limit*piece_time*piece_time-acc.norm()/(weight*weight);
                d=acc_limit-acc.norm()/(weight*weight*time_weight[sp_id]*time_weight[sp_id]*piece_time*piece_time);
                if(d<=0)
                    return INFINITY;
                
                if(d<margin)
                { 
                    //std::cout<<"e-acc:"<<acc<<std::endl;
                   energy+=-weight*(d-margin)*(d-margin)*log(d/margin); 
                }
            }
        }

        return energy;  

    }

};

PRJ_END

#endif