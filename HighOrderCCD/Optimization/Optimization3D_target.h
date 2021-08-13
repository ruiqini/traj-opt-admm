#ifndef OPTIMIZATION3D_TARGET_H
#define OPTIMIZATION3D_TARGET_H

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

class Optimization3D_target
{
public:

  typedef Eigen::MatrixXd Data;
 
  typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double


  static void optimization(const std::vector<Eigen::Vector3d>& target_list,
                           std::vector<Data>& spline_list, std::vector<double>& piece_time_list,
                           std::vector<Data>& p_slack_list, std::vector<Eigen::VectorXd>& t_slack_list, 
                           std::vector<Data>& p_lambda_list, std::vector<Eigen::VectorXd>& t_lambda_list, 
                           const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                           BVH& bvh)
  {
    
      bool all_target=true;
      for(int i=0;i<uav_num;i++)
      {
        if((spline_list[i].row(trajectory_num-1).transpose() - target_list[i]).squaredNorm()>1e-2)
        {
          all_target=false;
        }
        else
        {
          reach_target[i]=true;
        }
      }    
      if(all_target)
         tnorm=1e-5;
      
         
      clock_t time_0 = clock();
      std::vector<std::vector<std::vector<Eigen::Vector3d>>> c_lists;  c_lists.resize(uav_num);
      std::vector<std::vector<std::vector<double>>> d_lists;  d_lists.resize(uav_num);

      for(int i=0;i<uav_num;i++)
      {

        std::vector<std::vector<Eigen::Vector3d>> c_list;
        std::vector<std::vector<double>> d_list;

        separate_plane(spline_list[i], V, F, c_list, d_list, bvh);

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
        Data direction(trajectory_num,3);
        double t_direction;

        if(reach_target[i])
        {
          direction.setZero();
          direction_list[i]=direction;
          t_direction_list[i]=0;
          continue;
        }
        
        spline_descent_direction( spline_list[i], direction,  piece_time_list[i],  t_direction,
                                  p_slack_list[i],  t_slack_list[i],
                                  p_lambda_list[i],  t_lambda_list[i],
                                  c_lists[i], d_lists[i]);
        
        direction_list[i]=direction;
        t_direction_list[i]=t_direction;
      }

      std::vector<double> step_list; 

      step_list=Step::self_step(spline_list, direction_list, bvh);

      for(int i=0;i<uav_num;i++)
      {
        if(reach_target[i])
        {
          continue;
        }

        double step=Step::position_step(spline_list[i], direction_list[i],V,F, bvh);  
        //double step=Step::mix_step(spline_list[i], direction_list[i],V,F, bvh, c_lists[i], d_lists[i]);  
        if(step<step_list[i])
           step_list[i]=step;
        spline_line_search( spline_list[i], direction_list[i],  piece_time_list[i], t_direction_list[i],
                            p_slack_list[i], t_slack_list[i],
                            p_lambda_list[i],  t_lambda_list[i],
                            V,F, bvh, 
                            c_lists[i], d_lists[i],
                            step_list[i]);
      }
      
      clock_t time1 = clock();   

      //gn_file <<gnorm<<std::endl;
      //step_file<<max_step<<std::endl;
      for(int i=0;i<uav_num;i++)
      {
        if(reach_target[i])
        {
          continue;
        }

        update_slack_lambda(target_list[i],
                           spline_list[i],  piece_time_list[i],
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
                             const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
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

              int f0=F(ob_id,0); int f1=F(ob_id,1); int f2=F(ob_id,2);
              
              Eigen::Matrix3d _position;
              _position<<V.row(f0),V.row(f1),V.row(f2);

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
                  if(Separate::selfgjk(P0, P1, offset+2*margin,
                                      c, d))//cgal
                  {
                    
                    self_c_lists[p0][tr_id].push_back(c);
                    self_d_lists[p0][tr_id].push_back(d-0.5*offset);

                    self_c_lists[p1][tr_id].push_back(-c);
                    self_d_lists[p1][tr_id].push_back(-d-0.5*offset);
                  }
                }
            }
        }
  }

  static void update_slack_lambda(const Eigen::Vector3d& target,
                                  const Data& spline, const double& piece_time,
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

        double target_weight=1;

        Gradient_admm::slack_gradient( c_spline,  piece_time,
                                        p_part,  t_part, 
                                        p_lambda_part, t_lambda_part, 
                                        grad, hessian);
        if(sp_id == piece_num-1)
        {
          Eigen::VectorXd g0;
          Eigen::MatrixXd h0;
          Gradient_admm::target_gradient( p_part.row(order_num-2).transpose(),  target,
                                          g0,  h0);
          grad.segment<3>(3*(order_num-2))+=target_weight*g0;
          hessian.block<3,3>(3*(order_num-2),3*(order_num-2))+=target_weight*h0;
        }

        Eigen::VectorXd ng0;
        Eigen::MatrixXd h0;
        int t_n;
          
        Eigen::RowVector3d vel=1/t_part*(p_part.row(1)-p_part.row(0));
        
        if(sp_id==0)
        {
          
          t_n=order_num+1-2;

          for(int i=0;i<3;i++)
          {
            grad(3*(order_num+1))+=vel(i)*grad(3+i);
          }
          for(int i=0;i<3;i++)
          {
            hessian.col(3*(order_num+1))+=vel(i)*hessian.col(3+i);
          }
          for(int i=0;i<3;i++)
          {
            hessian.row(3*(order_num+1))+=vel(i)*hessian.row(3+i);
          }

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
          wolfe = x0.dot(ng0);

          Data direction;
          double t_direction;

          x=x0.head(t_n*3);
          
          t_direction=x0(t_n*3);

          t_direction=0;

          Eigen::MatrixXd d(Eigen::Map<Eigen::MatrixXd>(x.data(), 3,t_n));

          Eigen::MatrixXd d_=d.transpose();
          
          direction.resize(order_num+1,3);
          direction.setZero();

        if(sp_id==0)
        {
          direction.row(1)=t_direction*vel;
          direction.block(2,0,t_n,3)=d_;
        }
        else if(sp_id==piece_num-1)
        {
          direction.block(0,0,t_n,3)=d_;
          direction.row(order_num-1)=direction.row(order_num-2);
          direction.row(order_num)=direction.row(order_num-2);
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

        double e1=Energy_admm::slack_energy(c_spline, piece_time,
                                           p_part+step*direction, t_part, 
                                           p_lambda_part, t_lambda_part);
        double e=Energy_admm::slack_energy(c_spline, piece_time,
                                           p_part, t_part, 
                                           p_lambda_part, t_lambda_part);
        if(sp_id == piece_num-1)
        {
          e+=target_weight*Energy_admm::target_energy( p_part.row(order_num).transpose(),  target);
          e1+=target_weight*Energy_admm::target_energy( (p_part+step*direction).row(order_num).transpose(), target);
        }

        double init_time=t_part;
        t_part=init_time+step*t_direction;

        while(e-1e-4*wolfe*step<e1)
        {
          step*=0.8;
          t_part=init_time+step*t_direction;
          e1=Energy_admm::slack_energy(c_spline, piece_time,
                                        p_part+step*direction, t_part, 
                                        p_lambda_part, t_lambda_part);
          if(sp_id == piece_num-1)
          {
            e1+=target_weight*Energy_admm::target_energy( (p_part+step*direction).row(order_num).transpose(), target);
          }
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

    Gradient_admm::spline_gradient(spline,  piece_time,
                                   p_slack,  t_slack, 
                                   p_lambda,  t_lambda,
                                   c_list, d_list,
                                   grad, hessian);
    /*                               
    Gradient_admm::global_spline_gradient(spline,  piece_time,
                                          p_slack,  t_slack, 
                                          p_lambda,  t_lambda,
                                          c_list, d_list,
                                          grad, hessian);
    */
    Eigen::RowVector3d vel=1/piece_time*(spline.row(1)-spline.row(0));

    for(int i=0;i<3;i++)
    {
      grad(3*t_n)+=vel(i)*grad(3+i);
    }
    for(int i=0;i<3;i++)
    {
      hessian.col(3*t_n)+=vel(i)*hessian.col(3+i);
    }
    for(int i=0;i<3;i++)
    {
      hessian.row(3*t_n)+=vel(i)*hessian.row(3+i);
    }

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

    t_direction=0;


    Eigen::MatrixXd d(Eigen::Map<Eigen::MatrixXd>(x.data(), 3, t_n-4));

    Eigen::MatrixXd d_=d.transpose();
    
    direction.resize(t_n,3);
   
    direction.setZero();
    direction.block(2,0,t_n-4,3)=d_;

    direction.row(1)=t_direction*vel;
    
    direction.row(t_n-2)=direction.row(t_n-3);
    direction.row(t_n-1)=direction.row(t_n-3);

    
    std::cout<<"gn:"<<ng0.norm()<<std::endl;
    std::cout<<"dn:"<<x0.norm()<<std::endl;
    std::cout<<"t_direction:"<<t_direction<<std::endl;
    
    //gnorm=ng0.norm();
    //std::cout<<"gnorm:"<<gnorm<<std::endl;
  }

  static void spline_line_search(Data& spline, const Data& direction, double& piece_time, const double& t_direction,
                                 const Data& p_slack, const Eigen::VectorXd& t_slack,
                                 const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                                 const Eigen::MatrixXd & V,const Eigen::MatrixXi& F, BVH& bvh, 
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
