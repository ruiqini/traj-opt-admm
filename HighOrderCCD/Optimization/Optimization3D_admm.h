#ifndef OPTIMIZATION3D_ADMM_H
#define OPTIMIZATION3D_ADMM_H



#include "HighOrderCCD/Utils/CCDUtils.h"

#include "HighOrderCCD/CCD/CCD.h"

#include "HighOrderCCD/Energy.h"

#include "HighOrderCCD/Energy_admm.h"
#include "HighOrderCCD/Gradient_admm.h"

#include "HighOrderCCD/Step.h"
#include "HighOrderCCD/Separate.h"


PRJ_BEGIN

class Optimization3D_admm 
{
public:

  typedef Eigen::MatrixXd Data;
 
  typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

  static void optimization(Data& spline, double& piece_time,
                           Data& p_slack, Eigen::VectorXd& t_slack, 
                           Data& p_lambda, Eigen::VectorXd& t_lambda, 
                           const std::vector<Eigen::Matrix3d>& face_list,
                           BVH& bvh)
  {

      std::vector<std::vector<Eigen::Vector3d>> c_lists;
      std::vector<std::vector<double>> d_lists;
      clock_t time_0 = clock();

      separate_plane(spline, face_list, c_lists, d_lists, bvh);
      
      clock_t time0 = clock();
      
      std::cout<<"\nupdate spline\n";
      update_spline(spline,  piece_time,
                    p_slack, t_slack,
                    p_lambda,  t_lambda,
                    face_list, bvh,
                    c_lists, d_lists);
      
      clock_t time1 = clock();     
      
      std::cout<<"\nupdate slack&lambda\n";
      update_slack_lambda(spline,  piece_time,
                          p_slack, t_slack,
                          p_lambda,  t_lambda);

      std::cout<<"piece_time:"<<piece_time<<std::endl;

      //energy_file <<Energy::fast_whole_energy( spline, V, F, bvh)<<std::endl;
      
      clock_t time2 = clock();
      std::cout<<std::endl<<"separate:"<<(time0-time_0)/(CLOCKS_PER_SEC/1000)<<std::endl;
      std::cout<<"spline:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
      std::cout<<"slack&lambda:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl;

  }

  static void separate_plane(const Data& spline, 
                             const std::vector<Eigen::Matrix3d>& face_list,
                             std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                             std::vector<std::vector<double>>& d_lists,
                             BVH& bvh)
  {
        std::vector<std::vector<unsigned int>> collision_pairs;

        bvh.DCDCollision(spline, collision_pairs,offset+margin);

        c_lists.resize(subdivide_tree.size());
        d_lists.resize(subdivide_tree.size());
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {
            c_lists[tr_id].clear();
            d_lists[tr_id].clear();
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
            
            Eigen::MatrixXd P; P.noalias()=basis*bz;
            
            
            for(int i=0;i<collision_size;i++)
            {
              int ob_id=collision_pair[i];

              //int f0=F(ob_id,0); int f1=F(ob_id,1); int f2=F(ob_id,2);
              
              Eigen::Matrix3d _position=face_list[ob_id];

              Eigen::Vector3d c;
              double d;
              if(CCD::KDOPDCD(P, _position,offset+margin))
              {
                if(Separate::opengjk(P, _position, offset+margin,
                                      c, d))//cgal
                {
          
                  c_lists[tr_id].push_back(c);
                  d_lists[tr_id].push_back(d);
                }
              }


            }     

        }

  }

  static void update_spline(Data& spline, double& piece_time,
                            const Data& p_slack, const Eigen::VectorXd& t_slack,
                            const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                            const std::vector<Eigen::Matrix3d>& face_list, BVH& bvh,
                            const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                            const std::vector<std::vector<double>>& d_lists)
  {
    Data direction;
    double t_direction;
            
    //clock_t time0 = clock();
    
    spline_descent_direction( spline, direction,  piece_time,  t_direction,
                              p_slack,  t_slack,
                              p_lambda,  t_lambda,
                              c_lists, d_lists);

    //clock_t time1 = clock();
                              
    spline_line_search( spline, direction,  piece_time, t_direction,
                        p_slack, t_slack,
                        p_lambda,  t_lambda,
                        face_list, bvh, 
                        c_lists, d_lists);
    //clock_t time2 = clock();
    
   
    //std::cout<<"descent_direction:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
    //std::cout<<"line_search:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl;

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
    
          Eigen::LLT<Eigen::MatrixXd> solver; 
          Eigen::MatrixXd I=h0; I.setIdentity();
          
          solver.compute(h0);
          
          if(solver.info() == Eigen::NumericalIssue)//50
          {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(h0);
            Eigen::MatrixXd eigenvalue=eigensolver.eigenvalues();
            if(eigenvalue(0)<0)
            {
              std::cout<<"eigenvalue:"<<eigenvalue(0)<<"-----------------------------"<<std::endl;
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
        //std::cout<<step<<"\n";
        p_part=p_part+step*direction;

        p_slack.block<order_num+1,3>(sp_id*(order_num+1),0)=p_part;
        t_slack(sp_id)=t_part;

        p_lambda.block<order_num+1,3>(sp_id*(order_num+1),0)+=mu*(c_spline-p_part);

        t_lambda(sp_id)+=mu*(piece_time-t_part);
        
    }
  }

  static int spline_descent_direction(const Data& spline, Data& direction, const double& piece_time, double& t_direction,
                                       const Data& p_slack, const Eigen::VectorXd& t_slack,
                                       const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                                       const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                       const std::vector<std::vector<double>>& d_lists)
  {
    int t_n=trajectory_num;

    Eigen::VectorXd grad, partgrad;
    Eigen::MatrixXd hessian;
    double g_t,h_t;
    /*
    Gradient_admm::spline_gradient(spline,  piece_time,
                                   p_slack,  t_slack, 
                                   p_lambda,  t_lambda,
                                   c_lists, d_lists,
                                   grad, hessian);
                                   */
                           
    Gradient_admm::global_spline_gradient(spline,  piece_time,
                                          p_slack,  t_slack, 
                                          p_lambda,  t_lambda,
                                          c_lists, d_lists,
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
    clock_t time0 = clock();

    
    /*
    Eigen::LLT<Eigen::MatrixXd> solver; 
    
    Eigen::MatrixXd I=h0; I.setIdentity();

    solver.compute(h0);

    if(solver.info() == Eigen::NumericalIssue)//50
    {
      
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(h0);
      Eigen::MatrixXd eigenvalue=eigensolver.eigenvalues();
      if(eigenvalue(0)<0)
      {
        std::cout<<"eigenvalue:"<<eigenvalue(0)<<"-----------------------------"<<std::endl;
        h0=h0-eigenvalue(0)*I+0.01*I;
      }
      solver.compute(h0);    
      
      //return 0;
    }
    */
    
    SpMat H=h0.sparseView();
    Eigen::SimplicialLLT<SpMat> solver;  // performs a Cholesky factorization of A
    solver.compute(H);
    

    x0 = solver.solve(ng0);
    //x0=ng0;
    wolfe=x0.dot(ng0);
    
    clock_t time1 = clock();
    std::cout<<"\ntime solve:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;

    x=x0.head((t_n-4)*3);
    
    t_direction=x0((t_n-4)*3);

    Eigen::MatrixXd d(Eigen::Map<Eigen::MatrixXd>(x.data(), 3,t_n-4));

    Eigen::MatrixXd d_=d.transpose();
    
    direction.resize(t_n,3);
   
    direction.setZero();
    direction.block(2,0,t_n-4,3)=d_;
    
    std::cout<<"gn:"<<ng0.norm()<<std::endl;
    std::cout<<"dn:"<<x0.norm()<<std::endl;
    std::cout<<"t_direction:"<<t_direction<<std::endl<<std::endl;
    
    gnorm=ng0.norm();

    return 1;
    //std::cout<<"gnorm:"<<gnorm<<std::endl;
  }

  static void spline_line_search(Data& spline, const Data& direction, double& piece_time, const double& t_direction,
                                 const Data& p_slack, const Eigen::VectorXd& t_slack,
                                 const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                                 const std::vector<Eigen::Matrix3d>& face_list, BVH& bvh, 
                                 const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                 const std::vector<std::vector<double>>& d_lists)
  {
    clock_t time0 = clock();

    double step=Step::position_step(spline, direction,face_list, bvh);

    clock_t time1 = clock();

    std::cout<<"\ntime ccd:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
    
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
                                        c_lists, d_lists);
    double init_time=piece_time;
    piece_time=init_time+step*t_direction;
    while(e-1e-4*wolfe*step<Energy_admm::spline_energy(spline+step*direction, piece_time,
                                                        p_slack,  t_slack, 
                                                        p_lambda,  t_lambda,
                                                        c_lists, d_lists))
    {
      step*=0.8;
      piece_time=init_time+step*t_direction;
    }
    

    max_step=step;
    
    std::cout<<"step:"<<step<<std::endl;
    std::cout<<"result:"<<Energy::dynamic_energy(spline+step*direction,piece_time)<<
               " "<<lambda*Energy::plane_barrier_energy(spline+step*direction,c_lists,d_lists)<<
               " "<<lambda*Energy::bound_energy(spline+step*direction,piece_time)<<
               " "<<kt*whole_weight*piece_time<<std::endl<<std::endl;
    //<<" "<<limit_energy(spline+step*direction)
    spline=spline+step*direction;

  }

};


PRJ_END

#endif
