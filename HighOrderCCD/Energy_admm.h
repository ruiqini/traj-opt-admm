#ifndef ENERGY_ADMM_H
#define ENERGY_ADMM_H

#include "Utils/CCDUtils.h"

//#include "HighOrderCCD/BVH/BVH.h"
//#include "HighOrderCCD/CCD/CCD.h"
//#include "HighOrderCCD/Distance.h"
//#include "HighOrderCCD/Distance_der.h"

#include <vector>

PRJ_BEGIN

class Energy_admm
{
  public:
    typedef Eigen::MatrixXd Data;
    typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
    typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;
    
    // whole energy
    /*ks*dynamic_energy(spline,piece_time)+kt*whole_weight*piece_time+
              lambda*plane_barrier_energy(spline,c_lists, d_lists)+ lambda*bound_energy(spline,piece_time)*/
    static double spline_energy(const Data& spline, const double& piece_time,
                                const Data& p_slack, const Eigen::VectorXd& t_slack, 
                                const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                                const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                const std::vector<std::vector<double>>& d_lists) //all piece
    { 
        double energy=lambda*plane_barrier_energy(spline,c_lists, d_lists)+ lambda*bound_energy(spline,piece_time);

        for(int sp_id=0;sp_id<piece_num;sp_id++)
        {
            
            Data p_delta = convert_list[sp_id]*spline.block<order_num+1,3>(sp_id*(order_num-2),0)
                            -p_slack.block<order_num+1,3>(sp_id*(order_num+1),0);
            Data lamb_part=p_lambda.block<order_num+1,3>(sp_id*(order_num+1),0);

            energy+=mu/2.0*p_delta.squaredNorm();
            energy+=mu/2.0*std::pow(piece_time-t_slack(sp_id),2);

            for(int j=0;j<3;j++)
            {
                Eigen::VectorXd x=p_delta.col(j);
                Eigen::VectorXd lamb=lamb_part.col(j);
                energy+=lamb.transpose()*x;
            }
            energy+=t_lambda(sp_id)*(piece_time-t_slack(sp_id));
        }

        return energy;
    }

    static double slack_energy(const Data& c_spline, const double& piece_time,
                               const Data& p_part, const double& t_part, 
                               const Data& p_lambda, const double& t_lambda)// one piece energy
    { 
        double energy=dynamic_energy(p_part,t_part);

        energy+=mu/2.0*(c_spline-p_part).squaredNorm();
        energy+=mu/2.0*(piece_time-t_part)*(piece_time-t_part);

        for(int j=0;j<3;j++)
        {
            Eigen::VectorXd x=(c_spline-p_part).col(j);
            Eigen::VectorXd lamb=p_lambda.col(j);
            energy+=lamb.transpose()*x;
        }
        energy+=t_lambda*(piece_time-t_part);

        return energy;
    }

    static double target_energy(const Eigen::Vector3d& endpoint, const Eigen::Vector3d& target)
    {        
        double energy=0.5*(endpoint-target).squaredNorm();
        
        return energy;
    }

    static double dynamic_energy(const Data& p_part, const double& t_part)
    {
        double energy=0;

        //dynamic energy
        
        for(int j=0;j<3;j++)
        {
            Eigen::VectorXd x=p_part.col(j);
            energy+=ks/std::pow(t_part,2*der_num-1)*0.5*x.transpose()*M_dynamic*x;
        }

        energy+=kt*t_part;
        //energy+=kt*t_part*t_part;

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
            
            Eigen::MatrixXd P=basis*bz;

            double d;

            for(unsigned int k=0;k<c_list.size();k++)
            {
                for(int j=0;j<=order_num;j++)
                {
                    //d=P[j].dot(c_list[k])+d_list[k];
                    d=P.row(j).dot(c_list[k])+d_list[k];

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
            
            Eigen::MatrixXd P=basis*bz;

            double d;
           
            for(int j=0;j<order_num;j++)
            {
                //Eigen::RowVector3d vel=order_num*(P[j+1]-P[j]);
                Eigen::RowVector3d vel=order_num*(P.row(j+1)-P.row(j));
                  //d=vel_limit*piece_time-vel.norm()/weight;

                //d=vel_limit-vel.norm()/(weight*piece_time);
                d=vel_limit*vel_limit-vel.squaredNorm()/std::pow(weight*piece_time,2);
                if(d<=0)
                    return INFINITY;
                double vel_margin=2*vel_limit*margin;
                
                if(d<vel_margin)
                { 
                   energy+=-weight*(d-vel_margin)*(d-vel_margin)*log(d/vel_margin); 
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
            
            Eigen::MatrixXd P=basis*bz;

            double d;
           
            for(int j=0;j<order_num-1;j++)
            {
                //Eigen::RowVector3d acc=order_num*(order_num-1)*(P[j+2]-2*P[j+1]+P[j]);
                Eigen::RowVector3d acc=order_num*(order_num-1)*(P.row(j+2)-2*P.row(j+1)+P.row(j));
                //d=acc_limit*piece_time*piece_time-acc.norm()/(weight*weight);
                
                //d=acc_limit-acc.norm()/(weight*weight*piece_time*piece_time);
                d=acc_limit*acc_limit-acc.squaredNorm()/std::pow(weight*weight*piece_time*piece_time,2);
                if(d<=0)
                    return INFINITY;
                                
                double acc_margin=2*acc_limit*margin;

                if(d<acc_margin)
                { 
                    //std::cout<<"e-acc:"<<acc<<std::endl;
                   energy+=-weight*(d-acc_margin)*(d-acc_margin)*log(d/acc_margin); 
                }
            }
        }

        return energy;  

    }
    

};

PRJ_END

#endif