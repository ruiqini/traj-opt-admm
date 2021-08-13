#ifndef GRADIENT_ADMM_H
#define GRADIENT_ADMM_H

#include "Utils/CCDUtils.h"

PRJ_BEGIN

class Gradient_admm
{
  public:
    typedef Eigen::MatrixXd Data;
    
    static void global_spline_gradient(const Data& spline, const double& piece_time,
                                        const Data& p_slack, const Eigen::VectorXd& t_slack, 
                                        const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                                        const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                        const std::vector<std::vector<double>>& d_lists,
                                        Eigen::VectorXd& grad, Eigen::MatrixXd& hessian)//all piece no fix
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
            
            local_spline_gradient(spline,  piece_time,
                                    p_slack,  t_slack, 
                                    p_lambda,  t_lambda,
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

    static void local_spline_gradient(const Data& spline, const double& piece_time,
                                        const Data& p_slack, const Eigen::VectorXd& t_slack, 
                                        const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                                        const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                        const std::vector<std::vector<double>>& d_lists,
                                        Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                                        int sp_id)//all piece no fix
    { 
        //std::cout<<"begin\n";
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
            
            local_plane_barrier_gradient( tr_id,  spline, 
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

            local_bound_gradient( tr_id, spline,  piece_time,
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

        Eigen::Matrix3d I;
        I.setIdentity();

        Eigen::MatrixXd B,M; 

       
        Data p_delta = convert_list[sp_id]*spline.block<order_num+1,3>(init,0)
                        -p_slack.block<order_num+1,3>(sp_id*(order_num+1),0);

        Data lamb_part=p_lambda.block<order_num+1,3>(sp_id*(order_num+1),0);
        
        Eigen::MatrixXd x1=convert_list[sp_id].transpose()*p_delta;
        Eigen::MatrixXd x2=convert_list[sp_id].transpose()*lamb_part;

        M=convert_list[sp_id].transpose()*convert_list[sp_id];
        B = Eigen::kroneckerProduct(M,I);
        
        x1.transposeInPlace();
        x2.transposeInPlace();
        
        Eigen::VectorXd v1(Eigen::Map<Eigen::VectorXd>(x1.data(), num));
        Eigen::VectorXd v2(Eigen::Map<Eigen::VectorXd>(x2.data(), num));
        
        grad.segment(0,num)+=mu*v1+v2;
        
        hessian.block(0,0,num,num)+=mu*B;

        grad(num)+=mu*(piece_time-t_slack(sp_id)) + t_lambda(sp_id);
        hessian(num,num)+=mu;
    }

    static void spline_gradient(const Data& spline, const double& piece_time,
                                const Data& p_slack, const Eigen::VectorXd& t_slack, 
                                const Data& p_lambda, const Eigen::VectorXd& t_lambda,
                                const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                const std::vector<std::vector<double>>& d_lists,
                                Eigen::VectorXd& grad, Eigen::MatrixXd& hessian)//all piece no fix
    { 
        Eigen::VectorXd g0,g1,g2, partgrad;
        Eigen::MatrixXd h0,h1,h2;
        double g_t,h_t;
        
        plane_barrier_gradient(spline, c_lists, d_lists,g1, h1);

        bound_gradient(spline,piece_time, g2,h2, g_t, h_t, partgrad);
        
        g0=lambda*g1+lambda*g2;
        h0=lambda*h1+lambda*h2;

        g_t*=lambda;
        h_t*=lambda;

        Eigen::Matrix3d I;
        I.setIdentity();

        Eigen::MatrixXd B,M; 

        for(int sp_id=0;sp_id<piece_num;sp_id++)
        {
            int init=sp_id*(order_num-2);

            Data p_delta = convert_list[sp_id]*spline.block<order_num+1,3>(init,0)
                            -p_slack.block<order_num+1,3>(sp_id*(order_num+1),0);

            Data lamb_part=p_lambda.block<order_num+1,3>(sp_id*(order_num+1),0);
            
            Eigen::MatrixXd x1=convert_list[sp_id].transpose()*p_delta;
            Eigen::MatrixXd x2=convert_list[sp_id].transpose()*lamb_part;

            M=convert_list[sp_id].transpose()*convert_list[sp_id];
            B = Eigen::kroneckerProduct(M,I);
            
            x1.transposeInPlace();
            x2.transposeInPlace();
            
            Eigen::VectorXd v1(Eigen::Map<Eigen::VectorXd>(x1.data(), 3*(order_num+1)));
            Eigen::VectorXd v2(Eigen::Map<Eigen::VectorXd>(x2.data(), 3*(order_num+1)));
            
            g0.segment<3*(order_num+1)>(3*init)+=mu*v1+v2;
            
            h0.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init)+=mu*B;

            g_t+=mu*(piece_time-t_slack(sp_id)) + t_lambda(sp_id);
            h_t+=mu;
        }


        int n=3*trajectory_num;
        
        grad.resize(n+1); grad.setZero();
        grad.head(n)=g0;
        grad(n)=g_t;

        hessian.resize(n+1,n+1); hessian.setZero();
        hessian.block(0,0,n,n)=h0;

        hessian.block(0,n,n,1)=lambda*partgrad;
        hessian.block(n,0,1,n)=lambda*partgrad.transpose();
        hessian(n,n)=h_t;
    }

    static void plane_barrier_gradient(const Data& spline, 
                                       const std::vector<std::vector<Eigen::Vector3d>>& c_lists,
                                       const std::vector<std::vector<double>>& d_lists,
                                       Eigen::VectorXd& grad, Eigen::MatrixXd& hessian)
    {
        
        int num=3*trajectory_num;
        
        grad.resize(num);
        grad.setZero();

        hessian.resize(num,num);
        hessian.setZero();

        //Eigen::VectorXd auto_grad=grad;
        //Eigen::MatrixXd auto_hessian=hessian;
        //double dmin=margin;

        Eigen::VectorXd g;
        Eigen::MatrixXd h;
        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {        
            std::vector<Eigen::Vector3d> c_list=c_lists[tr_id];
            std::vector<double> d_list=d_lists[tr_id];
            if(c_list.size()==0)
                continue;

            local_plane_barrier_gradient( tr_id,  spline, 
                                          c_list, d_list,
                                          g,  h);
            
            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            int init=sp_id*(order_num-2);
            
            grad.segment(3*init,3*(order_num+1)) += g;
                                        
            hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init) += h;
            
        }
        //std::cout<<std::endl<<"dmin:"<<dmin<<std::endl;
    }
    
  
    static void bound_gradient(const Data& spline, const double& piece_time,
                               Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                               double& g_t, double& h_t, 
                               Eigen::VectorXd& partgrad)
    {
        int num=3*trajectory_num;
        
        grad.resize(num);
        grad.setZero();

        hessian.resize(num,num);
        hessian.setZero();

        g_t=0;
        h_t=0;

        partgrad.resize(num);
        partgrad.setZero();
    
        //double max_vel=0;
        //double max_acc=0;

        Eigen::VectorXd g; Eigen::MatrixXd h;
        double local_g_t,  local_h_t;
        Eigen::VectorXd local_partgrad;

        for(unsigned int tr_id=0;tr_id<subdivide_tree.size();tr_id++)
        {
            int sp_id=std::get<0>(subdivide_tree[tr_id]);            
            int init=sp_id*(order_num-2);

            local_bound_gradient( tr_id, spline, piece_time,
                                   g,  h,
                                   local_g_t,  local_h_t, 
                                   local_partgrad);
            g_t+=local_g_t;
                   
            h_t+=local_h_t;

            grad.segment(3*init,3*(order_num+1)) += g;
                    
            hessian.block<3*(order_num+1),3*(order_num+1)>(3*init,3*init) += h;
 
            partgrad.segment(3*init,3*(order_num+1)) += local_partgrad;
          
        }
        
        
    }
    
    static void local_plane_barrier_gradient( int tr_id, const Data& spline, 
                                              const std::vector<Eigen::Vector3d>& c_list,
                                              const std::vector<double>& d_list,
                                              Eigen::VectorXd& grad, Eigen::MatrixXd& hessian)
    {
        
        grad.resize(3*(order_num+1));
        grad.setZero();

        hessian.resize(3*(order_num+1),3*(order_num+1));
        hessian.setZero();
            

            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            int init=sp_id*(order_num-2);

            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(init,0);
            
            Eigen::MatrixXd P=basis*bz;
            double d;
            
            Eigen::Matrix3d I; I.setIdentity();
            
            for(int j=0;j<=order_num;j++)
            {
                Eigen::MatrixXd A=Eigen::kroneckerProduct(basis.row(j),I);
                A.transposeInPlace();
                for(unsigned int k=0;k<c_list.size();k++)
                {
                    //d=P[j].dot(c_list[k])+d_list[k];
                    d=P.row(j).dot(c_list[k])+d_list[k];
                    
                    if(d<margin)
                    { 
                        //std::cout<<A<<"\n";
                       Eigen::MatrixXd d_x=A * c_list[k];

                       double e1=-weight*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);

                       grad += e1*d_x;

                       double e2=-weight*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));
                                        
                       hessian += e2*d_x*d_x.transpose();

                    }
                }

            }
    }

    static void local_bound_gradient(int tr_id, const Data& spline, const double& piece_time,
                                     Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                                     double& g_t, double& h_t, 
                                     Eigen::VectorXd& partgrad)
    {
        grad.resize(3*(order_num+1));
        grad.setZero();

        hessian.resize(3*(order_num+1),3*(order_num+1));
        hessian.setZero();

        g_t=0;
        h_t=0;

        partgrad.resize(3*(order_num+1));
        partgrad.setZero();
    
        //double max_vel=0;
        //double max_acc=0;
        
            int sp_id=std::get<0>(subdivide_tree[tr_id]);
            double weight=std::get<1>(subdivide_tree[tr_id]).second-std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[tr_id]);
            
            int init=sp_id*(order_num-2);
            
            Eigen::MatrixXd bz;
            bz=spline.block<order_num+1,3>(init,0);
            
            Eigen::MatrixXd P=basis*bz;

            double d;
            
            Eigen::Matrix3d I; I.setIdentity();

           
            for(int j=0;j<order_num;j++)
            {
                //Eigen::RowVector3d P_=P[j+1]-P[j];
                Eigen::RowVector3d P_=P.row(j+1)-P.row(j);
                //Eigen::RowVector3d vel=order_num*P_;
                double d_=P_.norm();

                double v=order_num*d_/(weight);
                //d=vel_limit*piece_time-vel.norm()/weight;
                d=vel_limit-v/piece_time;

                //if(v/piece_time>max_vel)
                //  max_vel=v/piece_time;
                
                if(d<margin)
                { 
                   double e1=-weight*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);
                   
                   double e2=-weight*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));

                   //g_t, h_t
                   
                   g_t+=e1*v/std::pow(piece_time,2);
                   
                   h_t+=-2*e1*v/std::pow(piece_time,3)+
                        +e2*v*v/std::pow(piece_time,4);
                   
                    Eigen::MatrixXd A=Eigen::kroneckerProduct(basis.row(j+1),I)-
                                      Eigen::kroneckerProduct(basis.row(j),I);
                    //std::cout<<A<<"\n";
                    Eigen::RowVector3d d_p;
                    Eigen::Matrix3d h_p;

                    
                    d_p=-order_num/(weight*piece_time)*P_/d_;
        
                    h_p=-order_num/(weight*piece_time)*(I/d_-P_.transpose()*P_/std::pow(d_,3));

                    Eigen::MatrixXd d_x=d_p*A;
                    
                    grad += e1*d_x.transpose();
                    
                    hessian += e2*d_x.transpose()*d_x+e1*A.transpose()*h_p*A;
                    
                    double e3=-e1/piece_time
                              +e2*(vel_limit-d)/piece_time;
                    
                    partgrad += e3*d_x.transpose();
                }
            }
           
            for(int j=0;j<order_num-1;j++)
            {
                //Eigen::RowVector3d P_=P[j+2]-2*P[j+1]+P[j];
                Eigen::RowVector3d P_=P.row(j+2)-2*P.row(j+1)+P.row(j);
                //Eigen::RowVector3d acc=order_num*(order_num-1)*P_;
                double d_ = P_.norm();
                double a = order_num*(order_num-1)*d_/(weight*weight);
                
                //d=acc_limit*piece_time*piece_time-acc.norm()/(weight*weight);
                d=acc_limit-a/(piece_time*piece_time);

                //if(a/(piece_time*piece_time)>max_acc)
                //  max_acc=a/(piece_time*piece_time);
                
                if(d<margin)
                { 
                   double e1=-weight*(2*(d-margin)*log(d/margin)+(d-margin)*(d-margin)/d);
                   
                   double e2=-weight*(2*log(d/margin)+4*(d-margin)/d-(d-margin)*(d-margin)/(d*d));
                   
                   //g_t, h_t
                   g_t+=2*e1*a/std::pow(piece_time,3);
                   
                   h_t+=-6*e1*a/std::pow(piece_time,4)+
                        +4*e2*a*a/std::pow(piece_time,6);
                                           
                    //Eigen::Matrix3d I; I.setIdentity();
                    Eigen::MatrixXd A=Eigen::kroneckerProduct(basis.row(j+2),I)-
                                      2 * Eigen::kroneckerProduct(basis.row(j+1),I)+
                                      Eigen::kroneckerProduct(basis.row(j),I);
                    //std::cout<<A<<"\n";
                    Eigen::RowVector3d d_p;
                    Eigen::Matrix3d h_p;
                    
                    d_p=-order_num*(order_num-1)/std::pow(weight*piece_time,2)*P_/d_;
        
                    h_p=-order_num*(order_num-1)/std::pow(weight*piece_time,2)*(I/d_-P_.transpose()*P_/std::pow(d_,3));

                    Eigen::MatrixXd d_x=d_p*A;
                    
                    grad += e1*d_x.transpose();
                                        
                    hessian+=e2*d_x.transpose()*d_x+e1*A.transpose()*h_p*A;
                   
                    double e3=-2*e1/piece_time
                              +2*e2*(acc_limit-d)/piece_time;
                   
                    partgrad += e3*d_x.transpose();
                }
            }
        

      

        //std::cout<<"max_vel:"<<max_vel<<std::endl;
        //std::cout<<"max_acc:"<<max_acc<<std::endl;
    }

    static void slack_gradient(const Data& c_spline, const double& piece_time,
                               const Data& p_part, const double& t_part, 
                               const Data& p_lambda, const double& t_lambda,
                               Eigen::VectorXd& grad, Eigen::MatrixXd& hessian)// one piece no fix
    { 
        Eigen::VectorXd g0,g1,g2, partgrad;
        Eigen::MatrixXd h0,h1,h2;
        double g_t,h_t;

        dynamic_gradient(p_part, t_part, 
                         g1, h1,
                         g_t, h_t, 
                         partgrad);
    
        int n=3*(order_num+1);

        Eigen::MatrixXd x1=p_part-c_spline;
        Eigen::MatrixXd x2=p_lambda;
        
        x1.transposeInPlace();
        x2.transposeInPlace();
            
        Eigen::VectorXd v1(Eigen::Map<Eigen::VectorXd>(x1.data(), n));
        Eigen::VectorXd v2(Eigen::Map<Eigen::VectorXd>(x2.data(), n));
        
        g2=mu*v1-v2;

        h2.resize(n,n); 
        h2.setIdentity();
        h2*=mu;

        g0=g1+g2;
        h0=h1+h2;

        g_t+=mu*(t_part-piece_time)-t_lambda;
        h_t+=mu;

        grad.resize(n+1); grad.setZero();
        grad.head(n)=g0;
        grad(n)=g_t;

        hessian.resize(n+1,n+1); hessian.setZero();
        hessian.block(0,0,n,n)=h0;

        hessian.block(0,n,n,1)=partgrad;
        hessian.block(n,0,1,n)=partgrad.transpose();
        hessian(n,n)=h_t;

    }

    static void target_gradient(const Eigen::Vector3d& endpoint, const Eigen::Vector3d& target,
                                Eigen::VectorXd& grad, Eigen::MatrixXd& hessian)
    {
        Eigen::Matrix3d I;
        I.setIdentity();
        grad=endpoint-target;
        hessian=I;
    }

    static void dynamic_gradient(const Data& p_part, const double& t_part, 
                                 Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
                                 double& g_t, double& h_t, 
                                 Eigen::VectorXd& partgrad) //one piece
    {
        double dynamic_energy=0;

        //dynamic energy
        Eigen::Matrix3d I;
        I.setIdentity();

        Eigen::MatrixXd B = Eigen::kroneckerProduct(M_dynamic,I);
        Eigen::MatrixXd x=M_dynamic*p_part;

        x.transposeInPlace();
        Eigen::VectorXd v1(Eigen::Map<Eigen::VectorXd>(x.data(), 3*(order_num+1)));

        grad = ks/std::pow(t_part,2*der_num-1)*v1;
        hessian = ks/std::pow(t_part,2*der_num-1)*B;
        
        for(int j=0;j<3;j++)
        {
            Eigen::VectorXd x1=p_part.col(j);
            dynamic_energy+=ks/std::pow(t_part,2*der_num-1)*0.5*x1.transpose()*M_dynamic*x1;
        }

        //dynamic_energy+=kt*t_part;
        
        g_t=-(2*der_num-1)*dynamic_energy/t_part;
        //g_t+=kt;
        g_t+=kt*1.1*std::pow(t_part,0.1);

        h_t=(2*der_num-1)*(2*der_num)*dynamic_energy/(t_part*t_part);  
        
        h_t+=kt*0.11*std::pow(t_part,-0.9);

        partgrad=-(2*der_num-1)*grad/t_part;

    }
    
};

PRJ_END

#endif