#ifndef STEP_H
#define STEP_H

#include "Utils/CCDUtils.h"

#include "BVH/BVH.h"
#include "CCD/CCD.h"
#include <vector>


PRJ_BEGIN

class Step
{
  public:
    typedef Eigen::MatrixXd Data;
    typedef std::pair<unsigned int, unsigned int> id_pair;

    

  
    static double position_step(const Data& spline, const Data& direction, 
                                const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                                BVH& bvh)
    {
        double step;
        if(step_choose)
        {
            step=1.0;
        }
        else
        {
            step=std::min(1.5*max_step,1.0);
        }

        //clock_t time1 = clock();
        //std::cout<<std::endl<<"bvhtime:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
        //bool is_collided=true;
        //std::cout<<"bvh: "<<collision_pair.size()<<std::endl;
        std::vector<std::vector<unsigned int>> collision_pairs;
        //bvh.CheckCollision(collision_pair,margin);
        bvh.CCDCollision(spline, direction, collision_pairs,offset);
        //bvh.CCDCollision(spline, direction, collision_pairs,offset);
        
        //double temp_step=step;
        //double step0=0.0;
        
        //std::cout<<"start:"<<collision_pair.size()<<std::endl;
        //time0 = clock();
        unsigned int whole_size=0;
        for(unsigned int tr_id=0;tr_id<collision_pairs.size();tr_id++)
        {
            whole_size+=collision_pairs[tr_id].size();
        }
        if(whole_size==0)
          return step;
        
        std::cout<<"bvh: "<<whole_size<<std::endl;

            
        for(unsigned int tr_id=0;tr_id<collision_pairs.size();tr_id++)
        {

            //std::vector<unsigned int> collision_pair=temp_pairs[tr_id];
            int collision_size =collision_pairs[tr_id].size();
            if(collision_size==0)
                continue;
            int sp_id=std::get<0>(subdivide_tree[tr_id]);

            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            Eigen::MatrixXd bz, bz_d;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
            bz_d=direction.block<order_num+1,3>(sp_id*(order_num-2),0);
            
            Eigen::MatrixXd P(order_num+1,3),D(order_num+1,3);
            P=basis*bz;
            D=basis*bz_d;
            //for(auto it = temp_pairs[tr_id].begin(); it != temp_pairs[tr_id].end();)
            for(int i=0;i<collision_size;i++)
            {
                //int ob_id=*it;
                int ob_id=collision_pairs[tr_id][i];

                int f0=F(ob_id,0); int f1=F(ob_id,1); int f2=F(ob_id,2);
                
                Eigen::Matrix3d _position;
                _position<<V.row(f0),V.row(f1),V.row(f2);

                                                    
                //std::cout<<"is_collided:"<<is_collided<<" \n";
                
                //is_collided=CCD::ExactCCD(P,D,_position);
                bool is_collided = CCD::KDOPCCD(P,D,_position,offset,0,step);
                while(is_collided)
                {  
                    //is_collided= CCD::KDOPCCD(P,D,_position,offset,step0,step0+temp_step);
                    //if(is_collided)
                    //{
                        is_collided= CCD::GJKCCD(P,D,_position, offset,0,step);  //cgal
                        if(is_collided)
                        {
                            step*=0.8;
                        }     
                    //}
                }
            }
            
        }

        
        

        return step;

    }
    
    static std::vector<double> self_step( const std::vector<Data>& spline_list,const std::vector<Data>& direction_list, BVH& bvh)
    {
        std::vector<double> step_list; step_list.resize(uav_num,1.0);

        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            //double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            std::vector<Eigen::MatrixXd> P_list, D_list;
            for(int i=0;i<uav_num;i++)
            {
              Eigen::MatrixXd bz, bz_d;
              bz=spline_list[i].block<order_num+1,3>(sp_id*(order_num-2),0);
              bz_d=direction_list[i].block<order_num+1,3>(sp_id*(order_num-2),0);
           
              Eigen::MatrixXd P = basis*bz;
              Eigen::MatrixXd D = basis*bz_d;
              P_list.push_back(P);
              D_list.push_back(D);

            }
            std::vector<std::pair<unsigned int, unsigned int>> self_collision_pair;
            
            bvh.SelfCCDCollision(P_list, D_list, self_collision_pair, offset);

            for(unsigned int i=0;i<self_collision_pair.size();i++)
            {
                int p0=self_collision_pair[i].first;
                int p1=self_collision_pair[i].second;
                
                Eigen::MatrixXd P0=P_list[p0];
                Eigen::MatrixXd P1=P_list[p1];

                Eigen::MatrixXd D0=D_list[p0];
                Eigen::MatrixXd D1=D_list[p1];

                double temp_step0=step_list[p0];
                double temp_step1=step_list[p1];

                bool is_collided=CCD::SelfKDOPCCD(P0,D0,P1,D1,offset,0,temp_step0, 0, temp_step1);
                while(is_collided)
                {  
                    is_collided= CCD::SelfGJKCCD(P0,D0,P1,D1, offset,0,temp_step0, 0, temp_step1);  //cgal
                    if(is_collided)
                    {
                        temp_step0*=0.8;
                        temp_step1*=0.8;
                    }             
                }
                step_list[p0]=temp_step0;
                step_list[p1]=temp_step1;
            }
        }

        return step_list;

    }
    
    static double plane_step(const Data& spline, const Data& direction, 
                             const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                             const std::vector<std::vector<double>>& d_lists)
    {
        double step;
        if(step_choose)
        {
            step=1.0;
        }
        else
        {
            step=std::min(1.5*max_step,1.0);
        }
        double temp_step=step;
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {
            std::vector<Eigen::Vector3d> c_list=c_lists[tr_id];
            std::vector<double> d_list=d_lists[tr_id];
            if(c_list.size()==0)
                continue;

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            //double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            Eigen::MatrixXd bz, bz_d;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
            bz_d=direction.block<order_num+1,3>(sp_id*(order_num-2),0);
            
            Eigen::MatrixXd P(order_num+1,3),D(order_num+1,3);
            P=basis*bz;
            D=basis*bz_d;

            double d;

            for(unsigned int k=0;k<c_list.size();k++)
            {
                for(int j=0;j<=order_num;j++)
                {
                    //d=P[j].dot(c_list[k])+d_list[k];
                    d=(P+temp_step*D).row(j).dot(c_list[k])+d_list[k];

                    while(d<=0)
                    {
                        temp_step*=0.8;
                        d=(P+temp_step*D).row(j).dot(c_list[k])+d_list[k];
                    }
                        
                    
                }

            }
        }
        step=temp_step;

        return step;
        
    }


    static double mix_step(const Data& spline, const Data& direction, 
                           const Eigen::MatrixXd & V,const Eigen::MatrixXi& F, BVH& bvh,
                            const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                            const std::vector<std::vector<double>>& d_lists)
    {
        double step;
        if(step_choose)
        {
            step=1.0;
        }
        else
        {
            step=std::min(1.5*max_step,1.0);
        }

        bool is_collided=true;
        //std::cout<<"bvh: "<<collision_pair.size()<<std::endl;
        std::vector<std::vector<unsigned int>> collision_pairs;
        //bvh.CheckCollision(collision_pair,margin);
        bvh.CCDCollision(spline, direction, collision_pairs, offset);
    
        //std::cout<<"start:"<<collision_pair.size()<<std::endl;
        //time0 = clock();
        unsigned int whole_size=0;
        for(unsigned int tr_id=0;tr_id<collision_pairs.size();tr_id++)
        {
            whole_size+=collision_pairs[tr_id].size();
        }
        if(whole_size==0)
          return step;
        std::cout<<"bvh: "<<whole_size<<std::endl;

        double temp_step=step;
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {
            std::vector<Eigen::Vector3d> c_list=c_lists[tr_id];
            std::vector<double> d_list=d_lists[tr_id];

            int collision_size =collision_pairs[tr_id].size();
            if(c_list.size()==0 && collision_size==0)
                continue;

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            //double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            Eigen::MatrixXd bz, bz_d;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
            bz_d=direction.block<order_num+1,3>(sp_id*(order_num-2),0);
            
            Eigen::MatrixXd P(order_num+1,3),D(order_num+1,3);
            P=basis*bz;
            D=basis*bz_d;

            double d;

            for(unsigned int k=0;k<c_list.size();k++)
            {
                for(int j=0;j<=order_num;j++)
                {
                    //d=P[j].dot(c_list[k])+d_list[k];
                    d=(P+temp_step*D).row(j).dot(c_list[k])+d_list[k];

                    while(d<=0)
                    {
                        temp_step*=0.8;
                        d=(P+temp_step*D).row(j).dot(c_list[k])+d_list[k];
                    }
                        
                    
                }

            }
            if(c_list.size()==0 && collision_size!=0)
            for(int i=0;i<collision_size;i++)
            {
                int ob_id=collision_pairs[tr_id][i];

                int f0=F(ob_id,0); int f1=F(ob_id,1); int f2=F(ob_id,2);
                
                Eigen::Matrix3d _position;
                _position<<V.row(f0),V.row(f1),V.row(f2);

                is_collided=CCD::KDOPCCD(P,D,_position,offset,0,temp_step);
                while(is_collided)
                {  
                    //is_collided= CCD::KDOPCCD(P,D,_position,offset,step0,step0+temp_step);
                    //if(is_collided)
                    //{
                        is_collided= CCD::GJKCCD(P,D,_position, offset,0,temp_step);  //cgal
                        if(is_collided)
                        {
                            temp_step*=0.8;
                        }     
                    //}
                }
            }
        }

        step=temp_step;

        return step;
        
    }

};

PRJ_END

#endif