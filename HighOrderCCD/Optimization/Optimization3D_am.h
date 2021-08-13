#ifndef OPTIMIZATION3D_AM_H
#define OPTIMIZATION3D_AM_H



#include "HighOrderCCD/Utils/CCDUtils.h"

#include "HighOrderCCD/CCD/CCD.h"

#include "HighOrderCCD/Energy.h"

#include "HighOrderCCD/Energy_admm.h"
#include "HighOrderCCD/Gradient_admm.h"

#include "HighOrderCCD/Step.h"
#include "HighOrderCCD/Separate.h"



#include <vector>
#include <ctime>
#include <Eigen/SparseCholesky>


PRJ_BEGIN

class Optimization3D_am 
{
public:

  typedef Eigen::MatrixXd Data;
 
  typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

  static void optimization(Data& spline, double& piece_time,
                           const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
                           BVH& bvh)
  {

      std::vector<std::vector<Eigen::Vector3d>> c_lists;
      std::vector<std::vector<double>> d_lists;
      clock_t time_0 = clock();

      separate_plane(spline, V, F, c_lists, d_lists, bvh);
      
      clock_t time0 = clock();
      
      std::cout<<"\nupdate spline\n";
      update_spline(spline,  piece_time,
                    V, F, bvh,
                    c_lists, d_lists);
      
      clock_t time1 = clock();     
      
      

      std::cout<<"piece_time:"<<piece_time<<std::endl;

      //energy_file <<Energy::fast_whole_energy( spline, V, F, bvh)<<std::endl;
      
      clock_t time2 = clock();
      std::cout<<std::endl<<"separate:"<<(time0-time_0)/(CLOCKS_PER_SEC/1000)<<std::endl;
      std::cout<<"spline:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
  }

  static void separate_plane(const Data& spline, 
                             const Eigen::MatrixXd & V,const Eigen::MatrixXi& F,
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
                  c_lists[tr_id].push_back(c);
                  d_lists[tr_id].push_back(d);
                }
              }


            }     

        }

  }

  static void update_spline(Data& spline, double& piece_time,
                            const Eigen::MatrixXd & V,const Eigen::MatrixXi& F, BVH& bvh,
                            const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                            const std::vector<std::vector<double>>& d_lists)
  {
    Data direction;
    double t_direction;
            
    //clock_t time0 = clock();
    
    spline_descent_direction( spline, direction,  piece_time,  t_direction,
                              c_lists, d_lists);

    //clock_t time1 = clock();
                              
    spline_line_search( spline, direction,  piece_time, t_direction,
                        V,F, bvh, 
                        c_lists, d_lists);
    //clock_t time2 = clock();
 
  }

  static int spline_descent_direction(const Data& spline, Data& direction, const double& piece_time, double& t_direction,
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
    global_spline_gradient(spline,  piece_time,
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
                                 const Eigen::MatrixXd & V,const Eigen::MatrixXi& F, BVH& bvh, 
                                 const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                 const std::vector<std::vector<double>>& d_lists)
  {
    clock_t time0 = clock();

    double step=Step::position_step(spline, direction,V,F, bvh);
    //double step=Step::plane_step(spline, direction,c_lists, d_lists);
    //double step=Step::mix_step(spline, direction,V,F, bvh,c_lists, d_lists);
    
    clock_t time1 = clock();

    std::cout<<"\ntime ccd:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
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
    

    double e=spline_energy(spline, piece_time,
                           c_lists, d_lists);

    double init_time=piece_time;
    piece_time=init_time+step*t_direction;
    while(e-1e-4*wolfe*step<spline_energy(spline+step*direction, piece_time,
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

  static double spline_energy(  const Data& spline, const double& piece_time,
                                const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                const std::vector<std::vector<double>>& d_lists)
  {
    double energy=lambda*Energy_admm::plane_barrier_energy(spline,c_lists, d_lists) + 
                  lambda*Energy_admm::bound_energy(spline,piece_time);
    
    for(int sp_id=0;sp_id<piece_num;sp_id++)
    {
        int init=sp_id*(order_num-2);

        Data c_spline = convert_list[sp_id]*spline.block<order_num+1,3>(init,0);

        energy+=Energy_admm::dynamic_energy(c_spline, piece_time);
    }

    return energy;


  }

  static void global_spline_gradient( const Data& spline, const double& piece_time,
                                      const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                      const std::vector<std::vector<double>>& d_lists,
                                      Eigen::VectorXd& grad, Eigen::MatrixXd& hessian)
  {
    int n=3*trajectory_num;
        
    grad.resize(n+1); grad.setZero();
    hessian.resize(n+1,n+1); hessian.setZero();
    
    int num=3*(order_num+1);

    Eigen::MatrixXd I; I.resize(num+1,num+1); I.setIdentity();

    for(int sp_id=0;sp_id<piece_num;sp_id++)
    {
        int init=sp_id*(order_num-2);
        Eigen::VectorXd g0;
        Eigen::MatrixXd h0;
        
        local_spline_gradient( spline,  piece_time,
                               c_lists, d_lists,
                               g0,  h0, sp_id);
              
        Eigen::LLT<Eigen::MatrixXd> solver; 
    
        solver.compute(h0);
        
        if(solver.info() == Eigen::NumericalIssue)//50
        {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(h0);
            Eigen::MatrixXd eigenvalue=eigensolver.eigenvalues();
            if(eigenvalue(0)<0)
            {
              //std::cout<<"eigentheta:"<<eigenvalue(0)<<"-----------------------------"<<std::endl;
              h0=h0-eigenvalue(0)*I+0.01*I;
            }
        }

        grad.segment(3*init,num) += g0.segment(0,num);                
        hessian.block(3*init,3*init,num,num)+=  h0.block(0,0,num,num);

        grad(n)+=g0(num);
        hessian(n,n)+=h0(num,num);

        hessian.block(3*init,n,num,1)+=  h0.block(0,num,num,1);
        hessian.block(n,3*init,1,num)+=  h0.block(num,0,1,num);
          
    }
    

  }

  static void local_spline_gradient( const Data& spline, const double& piece_time,
                                      const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                      const std::vector<std::vector<double>>& d_lists, 
                                      Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                                      int sp_id)
  {
        int num=3*(order_num+1);
        grad.resize(num+1); grad.setZero();
        hessian.resize(num+1,num+1); hessian.setZero();
        
        int init=sp_id*(order_num-2);
        
        Eigen::MatrixXd bz;
        bz=spline.block<order_num+1,3>(init,0);

        Eigen::VectorXd g;
        Eigen::MatrixXd h;
        //plane barrier
        for(int i=0;i<res;i++)
        {
            int tr_id=sp_id*res+i;

            std::vector<Eigen::Vector3d> c_list=c_lists[tr_id];
            std::vector<double> d_list=d_lists[tr_id];
            if(c_list.size()==0)
                continue;

            Gradient_admm::local_plane_barrier_gradient( tr_id,  spline, 
                                                        c_list, d_list,
                                                        g,  h);
            grad.segment(0,num) += g;

            hessian.block(0,0,num,num) += h;
             
        }
        //bound
        double g_t,  h_t;
        Eigen::VectorXd partgrad;
        for(int i=0;i<res;i++)
        {
            int tr_id=sp_id*res+i;
            
            Gradient_admm::local_bound_gradient( tr_id, spline,  piece_time,
                                                g,  h,
                                                g_t,  h_t, 
                                                partgrad);
                   
            grad.segment(0,num) += g;
            
            grad(num)+=g_t;

                    
            hessian.block(0,0,num,num)+=h;
                                        
            hessian.block(0,num,num,1)+=partgrad;
            hessian.block(num,0,1,num)+=partgrad.transpose();

            hessian(num,num)+=h_t;
            
        }
        
        grad*=lambda;
        hessian*=lambda;
        
        double dynamic_energy=0;

        //dynamic energy
        Eigen::Matrix3d I;
        I.setIdentity();

        Data c_spline = convert_list[sp_id]*spline.block<order_num+1,3>(init,0);

        Eigen::MatrixXd M=convert_list[sp_id].transpose()*M_dynamic*convert_list[sp_id];

        Eigen::MatrixXd B = Eigen::kroneckerProduct(M,I);
        Eigen::MatrixXd x = convert_list[sp_id].transpose()*M_dynamic*c_spline;

        x.transposeInPlace();
        Eigen::VectorXd v1(Eigen::Map<Eigen::VectorXd>(x.data(), 3*(order_num+1)));

        
        
        for(int j=0;j<3;j++)
        {
            Eigen::VectorXd x1=c_spline.col(j);
            dynamic_energy+=ks/std::pow(piece_time,2*der_num-1)*0.5*x1.transpose()*M_dynamic*x1;
        }

        //dynamic_energy+=kt*t_part;
        g=ks/std::pow(piece_time,2*der_num-1)*v1;
        
        grad.segment(0,num) += g;
        hessian.block(0,0,num,num) += ks/std::pow(piece_time,2*der_num-1)*B;

        grad(num)+=-(2*der_num-1)*dynamic_energy/piece_time;
        //g_t+=kt;
        grad(num)+=kt*1.1*std::pow(piece_time,0.1);

        hessian(num,num)+=(2*der_num-1)*(2*der_num)*dynamic_energy/(piece_time*piece_time);  
        
        hessian(num,num)+=kt*0.11*std::pow(piece_time,-0.9);

        hessian.block(0,num,num,1)+=-(2*der_num-1)*g/piece_time;
        hessian.block(num,0,1,num)+=-(2*der_num-1)*g.transpose()/piece_time;

  }

  
};


PRJ_END

#endif
