#ifndef OPTIMIZATION3D_MULTI_H
#define OPTIMIZATION3D_MULTI_H

#include "HighOrderCCD/Utils/CCDUtils.h"

#include "HighOrderCCD/Energy.h"

#include "HighOrderCCD/Energy_admm.h"
#include "HighOrderCCD/Gradient_admm.h"

#include "HighOrderCCD/Step.h"
#include "HighOrderCCD/Separate.h"

#include <vector>
#include <ctime>
#include <Eigen/SparseCholesky>


PRJ_BEGIN

class Optimization3D_multi
{
public:

  typedef Eigen::MatrixXd Data;
 
  typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

  static void optimization(std::vector<Data>& spline_list, std::vector<double>& piece_time_list,
                           std::vector<Data>& p_slack_list, std::vector<Eigen::VectorXd>& t_slack_list, 
                           std::vector<Data>& p_lambda_list, std::vector<Eigen::VectorXd>& t_lambda_list, 
                           const std::vector<Eigen::Matrix3d>& face_list,
                           BVH& bvh)
  {
   
      clock_t time_0 = clock();
      std::vector<std::vector<std::vector<Eigen::Vector3d>>> c_lists;  c_lists.resize(uav_num);
      std::vector<std::vector<std::vector<double>>> d_lists;  d_lists.resize(uav_num);

      for(int i=0;i<uav_num;i++)
      {
        std::vector<std::vector<Eigen::Vector3d>> c_list;
        std::vector<std::vector<double>> d_list;

        separate_plane(spline_list[i], face_list, c_list, d_list, bvh);

        c_lists[i]=c_list;
        d_lists[i]=d_list;
      }

      separate_self( spline_list, 
                      c_lists, d_lists, bvh);

      std::vector<Data> direction_list; direction_list.resize(uav_num);
      std::vector<double> t_direction_list; t_direction_list.resize(uav_num);
      
      clock_t time0 = clock();
      for(int i=0;i<uav_num;i++)
      {
        Data direction;
        double t_direction;
                
        spline_descent_direction( spline_list[i], direction,  piece_time_list[i],  t_direction,
                                  p_slack_list[i],  t_slack_list[i],
                                  p_lambda_list[i],  t_lambda_list[i],
                                  c_lists[i], d_lists[i]);
        
        direction_list[i]=direction;
        t_direction_list[i]=t_direction;
      }

      std::vector<double> step_list; 

      Step::self_step(spline_list, direction_list,step_list, bvh);

      for(int i=0;i<uav_num;i++)
      {
        double step=Step::position_step(spline_list[i], direction_list[i],face_list, bvh);  
        //double step=Step::mix_step(spline_list[i], direction_list[i],V,F, bvh, c_lists[i], d_lists[i]);  
        if(step<step_list[i])
           step_list[i]=step;
        spline_line_search( spline_list[i], direction_list[i],  piece_time_list[i], t_direction_list[i],
                            p_slack_list[i], t_slack_list[i],
                            p_lambda_list[i],  t_lambda_list[i],
                            c_lists[i], d_lists[i],
                            step_list[i]);
      }
      
      clock_t time1 = clock();   

      //gn_file <<gnorm<<std::endl;
      //step_file<<max_step<<std::endl;
      for(int i=0;i<uav_num;i++)
      {

        update_slack_lambda(spline_list[i],  piece_time_list[i],
                            p_slack_list[i], t_slack_list[i],
                            p_lambda_list[i],  t_lambda_list[i]);
      }
      
      for(int i=0;i<uav_num;i++)
      {
        std::cout<<"piece_time:"<<piece_time_list[i]<<std::endl;
      }
      //energy_file <<Energy::fast_whole_energy( spline, V, F, bvh)<<std::endl;
      
      clock_t time2 = clock();
      
      std::cout<<std::endl<<"time10:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
      std::cout<<"time21:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
      std::cout<<"separate:"<<(time0-time_0)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;

  }

  static void separate_plane(const Data& spline, 
                             const std::vector<Eigen::Matrix3d>& face_list,
                             std::vector<std::vector<Eigen::Vector3d>>& c_list,
                             std::vector<std::vector<double>>& d_list,
                             BVH& bvh)
  {
     
        std::vector<std::vector<unsigned int>> collision_pairs;

        bvh.DCDCollision(spline, collision_pairs,offset+margin);

        c_list.resize(subdivide_tree.size());
        d_list.resize(subdivide_tree.size());

        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {
            c_list[tr_id].clear();
            d_list[tr_id].clear();
        }

   
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {
            std::vector<unsigned int> collision_pair=collision_pairs[tr_id];
            int collision_size =collision_pair.size();
            if(collision_size==0)
              continue;

            int sp_id=std::get<0>(subdivide_tree[tr_id]);

            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
            
            Eigen::MatrixXd P =basis*bz;
            for(int i=0;i<collision_size;i++)
            {
              int ob_id=collision_pair[i];

              //int f0=F(ob_id,0); int f1=F(ob_id,1); int f2=F(ob_id,2);
              
              Eigen::Matrix3d _position=face_list[ob_id];
              //_position<<V.row(f0),V.row(f1),V.row(f2);

              Eigen::Vector3d c;
              double d;
              if(CCD::KDOPDCD(P, _position,offset+margin))
              {
                  if(Separate::opengjk(P, _position, offset+margin,
                                    c, d))//cgal
                  {
                    c_list[tr_id].push_back(c);
                    d_list[tr_id].push_back(d);
                  }
              }

            }                      
        }
  }

  static void separate_self( const std::vector<Data>& spline_list, 
                             std::vector<std::vector<std::vector<Eigen::Vector3d>>>& self_c_lists,
                             std::vector<std::vector<std::vector<double>>>& self_d_lists,
                             BVH& bvh)
  {

        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            //double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            std::vector<Eigen::MatrixXd> P_list;
            for(int i=0;i<uav_num;i++)
            {
              Eigen::MatrixXd bz;
              bz=spline_list[i].block<order_num+1,3>(sp_id*(order_num-2),0);
           
              Eigen::MatrixXd P = basis*bz;
              P_list.push_back(P);

            }
            
            std::vector<std::pair<unsigned int, unsigned int>> self_collision_pair;
            
            bvh.SelfDCDCollision(P_list, self_collision_pair, offset+2*margin);

            for(unsigned int i=0;i<self_collision_pair.size();i++)
            {
                int p0=self_collision_pair[i].first;
                int p1=self_collision_pair[i].second;
                Eigen::Vector3d c;
                double d;
                Eigen::MatrixXd P0=P_list[p0];
                Eigen::MatrixXd P1=P_list[p1];
               
                if(CCD::SelfKDOPDCD(P0, P1,offset+2*margin))
                {
                  if(is_optimal_plane)
                  {
                    if(is_self_seperate[tr_id][p0][p1]==false)
                    {
                      if(Separate::selfgjk(P0, P1, offset+2*margin,
                                    c, d))
                      {
                        self_seperate_c[tr_id][p0][p1]=c; 
                        self_seperate_d[tr_id][p0][p1]=d;
                        is_self_seperate[tr_id][p0][p1]=true;
                        
                      }

                    }
                    
                  }
                  else
                  {
                      if(Separate::selfgjk(P0, P1, offset+2*margin,
                                          c, d))//cgal
                      {
                        Optimal_plane::optimal_d(P0, P1, c, d);
                        //Optimal_plane::self_optimal_cd(P0,  P1, 
                        //                               c,  d);
                        self_c_lists[p0][tr_id].push_back(c);
                        self_d_lists[p0][tr_id].push_back(d-0.5*offset);

                        self_c_lists[p1][tr_id].push_back(-c);
                        self_d_lists[p1][tr_id].push_back(-d-0.5*offset);

                      }
                  }
                }
            }
            if(is_optimal_plane)
                for(int p0=0;p0<uav_num;p0++)
                {
                  for(int p1=p0+1;p1<uav_num;p1++)
                  {
                      if(is_self_seperate[tr_id][p0][p1]==true)
                      {
                          Eigen::MatrixXd P0=P_list[p0];
                          Eigen::MatrixXd P1=P_list[p1]; 
                          Eigen::Vector3d c;
                          double d;
                          //std::cout<<"already\n";
                          c=self_seperate_c[tr_id][p0][p1];
                          d=self_seperate_d[tr_id][p0][p1];

                          Optimal_plane::self_optimal_cd(P0,  P1, 
                                                      c,  d);
                          self_seperate_c[tr_id][p0][p1]=c;
                          self_seperate_d[tr_id][p0][p1]=d;

                          self_c_lists[p0][tr_id].push_back(c);
                          self_d_lists[p0][tr_id].push_back(d-0.5*offset);

                          self_c_lists[p1][tr_id].push_back(-c);
                          self_d_lists[p1][tr_id].push_back(-d-0.5*offset);
                      }
                  
                  }
                }
        }

       
  }

  static void update_slack_lambda(const Data& spline, const double& piece_time,
                                  Data& p_slack,  Eigen::VectorXd& t_slack,
                                  Data& p_lambda,  Eigen::VectorXd& t_lambda)
  {


    for(int sp_id=0;sp_id<piece_num;sp_id++)
    {
        int init=sp_id*(order_num-2);

        Data c_spline = convert_list[sp_id]*spline.block<order_num+1,3>(init,0);
        Data p_part = p_slack.block<order_num+1,3>(sp_id*(order_num+1),0);

        Data p_lambda_part=p_lambda.block<order_num+1,3>(sp_id*(order_num+1),0);

        double t_part=t_slack(sp_id);
        double t_lambda_part=t_lambda(sp_id);

        Eigen::VectorXd grad, partgrad;
        Eigen::MatrixXd hessian;

        Gradient_admm::slack_gradient( c_spline,  piece_time,
                                        p_part,  t_part, 
                                        p_lambda_part, t_lambda_part, 
                                        grad, hessian);
        Eigen::VectorXd ng0;
        Eigen::MatrixXd h0;
        int t_n;

        if(sp_id==0)
        {
          t_n=order_num+1-2;

          ng0.resize(t_n*3+1);
          h0.resize(t_n*3+1,t_n*3+1);

          ng0.head(t_n*3)=-grad.segment(6,t_n*3);
          ng0(t_n*3)=-grad(3*(order_num+1));

          h0.block(0,0,t_n*3,t_n*3)=hessian.block(6,6,t_n*3,t_n*3);
          h0(t_n*3,t_n*3)=hessian(3*(order_num+1),3*(order_num+1));//h_t;

          partgrad=hessian.block(6,3*(order_num+1),3*t_n,1);

          h0.block(0,t_n*3,t_n*3,1)=partgrad;
          h0.block(t_n*3,0,1,t_n*3)=partgrad.transpose();

        }
        else if(sp_id==piece_num-1)
        {
          t_n=order_num+1-2;

          ng0.resize(t_n*3+1);
          h0.resize(t_n*3+1,t_n*3+1);

          ng0.head(t_n*3)=-grad.segment(0,t_n*3);
          ng0(t_n*3)=-grad(3*(order_num+1));

          h0.block(0,0,t_n*3,t_n*3)=hessian.block(0,0,t_n*3,t_n*3);
          h0(t_n*3,t_n*3)=hessian(3*(order_num+1),3*(order_num+1));//h_t;

          partgrad=hessian.block(0,3*(order_num+1),3*t_n,1);

          h0.block(0,t_n*3,t_n*3,1)=partgrad;
          h0.block(t_n*3,0,1,t_n*3)=partgrad.transpose();

        }
        else
        {
          t_n=order_num+1;
          h0=hessian;
          ng0=-grad;

        }

          Eigen::VectorXd x, x0;
    
          Eigen::MatrixXd I=h0; I.setIdentity();
          Eigen::LLT<Eigen::MatrixXd> solver; 
          
          solver.compute(h0);
          
          if(solver.info() == Eigen::NumericalIssue)//50
          {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(h0);
            Eigen::MatrixXd eigenvalue=eigensolver.eigenvalues();
            if(eigenvalue(0)<0)
            {
              //std::cout<<"eigenvalue:"<<eigenvalue(0)<<std::endl;
              h0=h0-eigenvalue(0)*I+0.01*I;
            }
            solver.compute(h0);    

          }
         
          
          x0 = solver.solve(ng0);
          //x0=ng0;
          wolfe=x0.dot(ng0);

          Data direction;
          double t_direction;

          x=x0.head(t_n*3);
          
          t_direction=x0(t_n*3);

          Eigen::MatrixXd d(Eigen::Map<Eigen::MatrixXd>(x.data(), 3,t_n));

          Eigen::MatrixXd d_=d.transpose();
          

          direction.resize(order_num+1,3);
          direction.setZero();

        if(sp_id==0)
        {
          direction.block(2,0,t_n,3)=d_;
        }
        else if(sp_id==piece_num-1)
        {
          direction.block(0,0,t_n,3)=d_;
        }
        else
        {
          direction.block(0,0,t_n,3)=d_;
        }


        double step=1.0;

        if(t_part+step*t_direction<=0)
        {
          step=-0.95*t_part/t_direction;
        }

        double e=Energy_admm::slack_energy(c_spline, piece_time,
                                           p_part, t_part, 
                                           p_lambda_part, t_lambda_part);
        double init_time=t_part;
        t_part=init_time+step*t_direction;

        while(e-1e-4*wolfe*step<Energy_admm::slack_energy(c_spline, piece_time,
                                                          p_part+step*direction, t_part, 
                                                          p_lambda_part, t_lambda_part))
        {
          step*=0.8;
          t_part=init_time+step*t_direction;
        }

        p_part=p_part+step*direction;

        p_slack.block<order_num+1,3>(sp_id*(order_num+1),0)=p_part;
        t_slack(sp_id)=t_part;

        p_lambda.block<order_num+1,3>(sp_id*(order_num+1),0)+=mu*(c_spline-p_part);

        t_lambda(sp_id)+=mu*(piece_time-t_part);
        
    }
    
    
  }

  static void spline_descent_direction(const Data& spline, Data& direction, const double& piece_time, double& t_direction,
                                       const Data& p_slack, const Eigen::VectorXd& t_slack,
                                       const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                                       const std::vector<std::vector<Eigen::Vector3d>>& c_list,
                                       const std::vector<std::vector<double>>& d_list)
  {
    int t_n=trajectory_num;

    Eigen::VectorXd grad, partgrad;
    Eigen::MatrixXd hessian;
    double g_t,h_t;
    /*
    Gradient_admm::spline_gradient(spline,  piece_time,
                                   p_slack,  t_slack, 
                                   p_lambda,  t_lambda,
                                   c_list, d_list,
                                   grad, hessian);
    */
    Gradient_admm::global_spline_gradient(spline,  piece_time,
                                          p_slack,  t_slack, 
                                          p_lambda,  t_lambda,
                                          c_list, d_list,
                                          grad, hessian);

    g_t=grad(3*t_n);
    h_t=hessian(3*t_n,3*t_n);
    partgrad=hessian.block(6,3*t_n,3*(t_n-4),1);

    Eigen::VectorXd ng = -grad.segment(6,(t_n-4)*3);
    Eigen::MatrixXd h = hessian.block(6,6,(t_n-4)*3,(t_n-4)*3);
    
    Eigen::VectorXd ng0((t_n-4)*3+1);
    Eigen::MatrixXd h0((t_n-4)*3+1,(t_n-4)*3+1);
    ng0.head((t_n-4)*3)=ng;
    ng0((t_n-4)*3)=-g_t;

    h0.block(0,0,(t_n-4)*3,(t_n-4)*3)=h;
    h0((t_n-4)*3,(t_n-4)*3)=h_t;//h_t;

    h0.block(0,(t_n-4)*3,(t_n-4)*3,1)=partgrad;
    h0.block((t_n-4)*3,0,1,(t_n-4)*3)=partgrad.transpose();
    //std::cout<<h0<<"\n";
    Eigen::VectorXd x, x0;
    
    Eigen::MatrixXd I=h0; I.setIdentity();
    Eigen::LLT<Eigen::MatrixXd> solver; 
    
    solver.compute(h0);
    
    if(solver.info() == Eigen::NumericalIssue)//50
    {
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(h0);
      Eigen::MatrixXd eigenvalue=eigensolver.eigenvalues();
      if(eigenvalue(0)<0)
      {
        //std::cout<<"eigenvalue:"<<eigenvalue(0)<<std::endl;
        h0=h0-eigenvalue(0)*I+0.01*I;
      }
      solver.compute(h0);    

    }

    /*
    SpMat H=h0.sparseView();
    Eigen::SimplicialLLT<SpMat> solver;  // performs a Cholesky factorization of A
    solver.compute(H);
    */
   
    
    x0 = solver.solve(ng0);
    //x0=ng0;
    wolfe=x0.dot(ng0);

    x=x0.head((t_n-4)*3);
    
    t_direction=x0((t_n-4)*3);
    

    Eigen::MatrixXd d(Eigen::Map<Eigen::MatrixXd>(x.data(), 3,t_n-4));

    Eigen::MatrixXd d_=d.transpose();
    
    direction.resize(t_n,3);
   
    direction.setZero();
    direction.block(2,0,t_n-4,3)=d_;
    
    std::cout<<"gn:"<<ng0.norm()<<std::endl;
    std::cout<<"dn:"<<x0.norm()<<std::endl;
    std::cout<<"t_direction:"<<t_direction<<std::endl;
    
    gnorm=ng0.norm();
    //std::cout<<"gnorm:"<<gnorm<<std::endl;
  }

  static void spline_line_search(Data& spline, const Data& direction, double& piece_time, const double& t_direction,
                                 const Data& p_slack, const Eigen::VectorXd& t_slack,
                                 const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                                 const std::vector<std::vector<Eigen::Vector3d>>& c_list,
                                 const std::vector<std::vector<double>>& d_list,
                                 double& step)
  {
    clock_t time0 = clock();

    //double step=Step::position_step(spline, direction,V,F, bvh);
    //double step=Step::plane_step(spline, direction,c_list, d_list);
    //if(step0<step)
    //  step=step0;
    
    clock_t time1 = clock();

    std::cout<<"time ccd:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
    //double time_step=Step::time_step(spline, direction, bvh);
    //if(time_step<step)
      //step=time_step;
    
    std::cout<<"highcdd:"<<step<<std::endl<<std::endl;
    if(piece_time+step*t_direction<=0)
    {
      step=-0.95*piece_time/t_direction;
    }

    std::cout.precision(10);
   
    std::cout<<"wolfe:"<<wolfe<<std::endl;
    

    double e=Energy_admm::spline_energy(spline, piece_time,
                                        p_slack,  t_slack, 
                                        p_lambda,  t_lambda,
                                        c_list, d_list);
    double init_time=piece_time;
    piece_time=init_time+step*t_direction;
    while(e-1e-4*wolfe*step<Energy_admm::spline_energy(spline+step*direction, piece_time,
                                                        p_slack,  t_slack, 
                                                        p_lambda,  t_lambda,
                                                        c_list, d_list))
    {
      step*=0.8;
      piece_time=init_time+step*t_direction;
    }
    

    max_step=step;
    
    std::cout<<"step:"<<step<<std::endl;
    std::cout<<"result:"<<Energy::dynamic_energy(spline+step*direction,piece_time)<<
               " "<<lambda*Energy::plane_barrier_energy(spline+step*direction,c_list,d_list)<<
               " "<<lambda*Energy::bound_energy(spline+step*direction,piece_time)<<
               " "<<kt*whole_weight*piece_time<<std::endl<<std::endl;
    //<<" "<<limit_energy(spline+step*direction)
    spline=spline+step*direction;

  }
};


PRJ_END

#endif
