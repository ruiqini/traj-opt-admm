#ifndef OPTIMAL_PLANE_H
#define OPTIMAL_PLANE_H

#include "Utils/CCDUtils.h"

PRJ_BEGIN

class Optimal_plane
{
  public:
    typedef Eigen::MatrixXd Data;
   
    static void optimal_d(const Data& position, const Data & _position, 
                          const Eigen::Vector3d& c, double& d)
    {

      double energy;
      double grad_d, hessian_d;
      double dist;

      
      double step=1.0;
      while(true)
      {
        energy=0;
        grad_d=0;
        hessian_d=0;
        for(int j=0;j<position.rows();j++)
        {
            //d=P[j].dot(c_list[k])+d_list[k];
            dist=position.row(j).dot(c)+d-0.5*offset;
            if(dist<margin)
            { 
                energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 
                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);  

                double e1=-(2*(dist-margin)*log(dist/margin)+(dist-margin)*(dist-margin)/dist);

                double e2=-(2*log(dist/margin)+4*(dist-margin)/dist-(dist-margin)*(dist-margin)/(dist*dist)); 

                grad_d+=e1;
                hessian_d+=e2;
            }
        }
        for(int j=0;j<_position.rows();j++)
        {
            dist=-_position.row(j).dot(c)-d-0.5*offset;
            if(dist<margin)
            { 
                energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 
                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   

                double e1=-(2*(dist-margin)*log(dist/margin)+(dist-margin)*(dist-margin)/dist);

                double e2=-(2*log(dist/margin)+4*(dist-margin)/dist-(dist-margin)*(dist-margin)/(dist*dist)); 

                grad_d+=-e1;
                hessian_d+=e2;
            }
        }

        double direction_d=-grad_d/hessian_d;
        
        d=d+step*direction_d;
        //std::cout<<grad_d<<"\n";
        if(std::abs(grad_d)<1e-2)
          break;

      }
      //std::cout<<"\n";
    }


    static Eigen::Vector3d current_c(const Eigen::Vector3d& c, const Eigen::Vector3d& c0, const Eigen::Vector3d& c1,
                                     const double& theta, const double& phi)
    {
         Eigen::Vector3d temp_c=std::cos(theta)*c+std::sin(theta)*(std::cos(phi)*c0+std::sin(phi)*c1);
         //temp_c.normalize();
         return temp_c;
    }

    static double current_d(const Eigen::Vector3d& c, const Data & _position)
    {
        //double *c_data=c.data();
        
        double d=-c.dot(_position.row(0));
        

        d-=offset;
        return d;
    }
    
    static double barrier_energy(const Data& position, const Data & _position, 
                                 Eigen::Vector3d& c, double& d)
    {
      double energy=0;
      double dist;
      d=current_d(c,_position);
        for(int j=0;j<position.rows();j++)
        {
            //d=P[j].dot(c_list[k])+d_list[k];
            dist=position.row(j).dot(c)+d;
            if(dist<=0)
            {
                return INFINITY;
            }
            if(dist<margin)
            { 
                energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 

            }
        }
        
        return energy;

    }

    static void barrier_grad(const Data& position, const Data & _position, 
                                  Eigen::Vector3d& c, Eigen::Vector3d& c0, Eigen::Vector3d& c1, double& d,
                                  Eigen::Vector2d& grad, Eigen::Matrix2d& hessian)
    {
        double dist;
        Eigen::Matrix2d tmp_h;
        double p_c,p_c0,p_c1;
        //d=INFINITY;
        Eigen::RowVector3d tmp;
        
        tmp=-_position;
        //tmp/=id+1;
        
        //d-=offset;
        
        for(int j=0;j<position.rows();j++)
        {
            //d=P[j].dot(c_list[k])+d_list[k];
            dist=(position.row(j)+tmp).dot(c)-offset;
            if(dist<margin)
            { 
                //energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 
                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);  
                p_c=(position.row(j)+tmp).dot(c);
                p_c0=(position.row(j)+tmp).dot(c0);
                p_c1=(position.row(j)+tmp).dot(c1);

                double e1=-(2*(dist-margin)*log(dist/margin)+(dist-margin)*(dist-margin)/dist);

                double e2=-(2*log(dist/margin)+4*(dist-margin)/dist-(dist-margin)*(dist-margin)/(dist*dist)); 
                
                grad(0)+=e1*p_c0;
                grad(1)+=0;//e1*p_c1;

                tmp_h<<e2*p_c0*p_c0-e1*p_c,   e1*p_c1,  
                          e1*p_c1,            0        ;
                hessian+=tmp_h;
            }
        }
        
    }

    static void optimal_cd(const Data& position, const Data & _position, 
                             Eigen::Vector3d& c, double& d) //order_num+1, order_num+1
    {      
      Eigen::Vector3d c0, c1, temp_c;

      Eigen::Matrix2d hessian;
      Eigen::Vector2d grad, tmp_g;

      double theta, phi; // [-Pi/2,Pi/2] [-Pi/2,Pi/2]
      
      double step;
      while(true)
      //for(int kk=0;kk<200;kk++)
      {

        c0(0)=c(1); c0(1)=-c(0); c0(2)=0;
        c0.normalize();

        c1=c0.cross(c);
        c1.normalize();

        theta=0;
        phi=0;

        hessian.setZero();
        grad.setZero();

        barrier_grad( position, _position, 
                      c,  c0,  c1,  d,
                      grad,  hessian);
        
        if(grad.norm()<1e-2)
        {
          d=current_d(c, _position);
          //std::cout<<"BREAK:"<<c.norm()<<" "<<c.transpose()<<" "<<d<<"\n";

          break; 
        }
        

        Eigen::LLT<Eigen::Matrix2d> solver; 
        Eigen::Matrix2d I; I.setIdentity();
        hessian=hessian+1e-2*I; 
        solver.compute(hessian);
        //std::cout<<hessian<<"\n";
        if(solver.info() == Eigen::NumericalIssue)//50
        {
          
          Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(hessian);
          Eigen::MatrixXd eigenvalue=eigensolver.eigenvalues();
          if(eigenvalue(0)<0)
          {
            hessian=hessian-eigenvalue(0)*I+1e-8*I; 
          }
          solver.compute(hessian);

        }
       
   
        //std::cout<<hessian<<"\n\n";
        double e0,e1;
        double tmp_d, tmp_theta, tmp_phi;
        Eigen::Vector3d tmp_c;
        
        Eigen::Vector2d direction=-solver.solve(grad);
        //Eigen::Vector3d direction=-grad;
        //double w= -grad.dot(direction);
        double w= -grad.dot(direction);

        //std::cout<<w<<" direction:"<<direction.transpose()<<"\n";

        step=1.0;
        if(std::abs(direction(0))>0.5*M_PI || std::abs(direction(1))> 0.5*M_PI)
        {
          step=0.95*std::min(0.5*std::abs(M_PI/direction(0)),0.5*std::abs(M_PI/direction(1)));
        }
        //std::cout<<"step:"<<step<<"\n";
        tmp_c=current_c(c,  c0,  c1,
                          theta,  phi);
        tmp_d=d;
        tmp_theta=theta;
        tmp_phi=phi;
        e0=barrier_energy(position, _position, 
                                 tmp_c,  tmp_d);
        tmp_theta= theta+step*direction(0);    
        tmp_phi=phi+step*direction(1);
        tmp_c=current_c(c,  c0,  c1,
                       tmp_theta,  tmp_phi);  
        tmp_d=current_d(tmp_c, _position);                   
        e1=barrier_energy(position, _position, 
                                 tmp_c,  tmp_d);
        //std::cout<<c.norm()<<" "<<c.transpose()<<" "<<d<<"\n"; 

        while(e0-1e-4*w*step< e1)
        {
          step*=0.8;
          
          tmp_theta= theta+step*direction(0);    
          tmp_phi=phi+step*direction(1);
          tmp_c=current_c(c,  c0,  c1,
                          tmp_theta,  tmp_phi);   
          tmp_d=current_d(tmp_c, _position);                       
          e1=barrier_energy(position, _position, 
                                 tmp_c,  tmp_d);
        }
        /*
        std::cout<<position<<"\n\n";
        std::cout<<_position<<"\n\n";
        std::cout<<e0<<" "<<e1<<" "<<grad.norm()<<"\n";
        std::cout<<step<<" "<<tmp_theta<<" "<<tmp_phi<<"\n";
        std::cout<<tmp_c.transpose()<<"\n";
        */
        //std::cout<<kk<<" "<<grad.norm()<<" "<<step<<" "<<e1<<"\n";
        c=tmp_c;
        d=current_d(c, _position);
        /*
        if(std::isnan(d))
        {
          std::cout<<"err "<<c.transpose()<<" "<<d<<"\n";
          double t;
          std::cin>>t;
        }

        std::cout<<c.transpose()<<" "<<d<<"\n";
        */
        
        
        if(std::abs((e1-e0)/e0)<1e-1)
           break;
        //if(std::abs((e1-e0)/e0)<1e-4)
        //    break;
      }
      //std::cout<<"??\n";
    }
    
    /*
    static double barrier_energy(const Data& position, const Data & _position, 
                         Eigen::Vector3d& c, double& d)
    {
      double energy=0;
      double dist;
      
        for(int j=0;j<position.rows();j++)
        {
            //d=P[j].dot(c_list[k])+d_list[k];
            dist=position.row(j).dot(c)+d-0.5*offset;
            if(dist<=0)
            {
                return INFINITY;
            }
            if(dist<margin)
            { 
                energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 

            }
        }
        double mar=margin;
        for(int j=0;j<_position.rows();j++)
        {
            dist=-_position.row(j).dot(c)-d-0.5*offset;
            if(dist<=0)
            {
                return INFINITY;
            }
            if(dist<mar)
            { 
                energy+=-(dist-mar)*(dist-mar)*log(dist/mar); 
               
            }
        }
        return energy;

    }

    static void barrier_grad(const Data& position, const Data & _position, 
                                  Eigen::Vector3d& c, Eigen::Vector3d& c0, Eigen::Vector3d& c1, double& d,
                                  Eigen::Vector3d& grad, Eigen::Matrix3d& hessian)
    {
        double dist;
        Eigen::Matrix3d tmp_h;
        double p_c,p_c0,p_c1;

        //double mar=margin;
        for(int j=0;j<position.rows();j++)
        {
            //d=P[j].dot(c_list[k])+d_list[k];
            dist=position.row(j).dot(c)+d-0.5*offset;
            if(dist<margin)
            { 
                //energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 
                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);  
                p_c=position.row(j).dot(c);
                p_c0=position.row(j).dot(c0);
                p_c1=position.row(j).dot(c1);

                double e1=-(2*(dist-margin)*log(dist/margin)+(dist-margin)*(dist-margin)/dist);

                double e2=-(2*log(dist/margin)+4*(dist-margin)/dist-(dist-margin)*(dist-margin)/(dist*dist)); 
                
                grad(0)+=e1*p_c0;
                grad(1)+=0;//e1*p_c1;
                grad(2)+=e1;

                tmp_h<<e2*p_c0*p_c0-e1*p_c,   e1*p_c1,  e2*p_c0,
                          e1*p_c1,            0,         0,
                          e2*p_c0,            0,        e2;
                hessian+=tmp_h;
                
               
            }
        }

        double mar=margin;
        for(int j=0;j<_position.rows();j++)
        {
            dist=-_position.row(j).dot(c)-d-0.5*offset;
            if(dist<mar)
            { 
                //energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 
                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);  

                p_c=-_position.row(j).dot(c);
                p_c0=-_position.row(j).dot(c0);
                p_c1=-_position.row(j).dot(c1); 

                double e1=-(2*(dist-mar)*log(dist/mar)+(dist-mar)*(dist-mar)/dist);

                double e2=-(2*log(dist/mar)+4*(dist-mar)/dist-(dist-mar)*(dist-mar)/(dist*dist)); 
                
                grad(0)+=e1*p_c0;
                grad(1)+=0;//e1*p_c1;
                grad(2)+=-e1;
               
                tmp_h<<e2*p_c0*p_c0-e1*p_c,   e1*p_c1,  -e2*p_c0,
                          e1*p_c1,            0,         0,
                          -e2*p_c0,            0,        e2;
                hessian+=tmp_h;
                
             
            }
        }

    }

    static void optimal_cd(const Data& position, const Data & _position, 
                             Eigen::Vector3d& c, double& d) //order_num+1, order_num+1
    {      
      Eigen::Vector3d c0, c1, temp_c;

      Eigen::Matrix3d hessian;
      Eigen::Vector3d grad, tmp_g;

      double theta, phi; // [-Pi/2,Pi/2] [-Pi/2,Pi/2]
      
      double step;
      while(true)
      //for(int kk=0;kk<200;kk++)
      {

        c0(0)=c(1); c0(1)=-c(0); c0(0)=0;
        c0.normalize();

        c1=c0.cross(c);

        theta=0;
        phi=0;

        hessian.setZero();
        grad.setZero();

        barrier_grad( position, _position, 
                           c,  c0,  c1,  d,
                           grad,  hessian);

        //std::cout<<grad.norm()<<" ";
        if(grad.norm()<1e-2)
          break;
        Eigen::LLT<Eigen::Matrix3d> solver; 
        Eigen::Matrix3d I; I.setIdentity();

        solver.compute(hessian);
        //std::cout<<hessian<<"\n";
        if(solver.info() == Eigen::NumericalIssue)//50
        {
          
          Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(hessian);
          Eigen::MatrixXd eigenvalue=eigensolver.eigenvalues();
          if(eigenvalue(0)<0)
          {
            hessian=hessian-eigenvalue(0)*I+1e-8*I;
          }
          solver.compute(hessian);

        }
       

        //std::cout<<hessian<<"\n\n";
        double e0,e1;
        double tmp_d, tmp_theta, tmp_phi;
        Eigen::Vector3d tmp_c;
        
        
       
        
        Eigen::Vector3d direction=-solver.solve(grad);
        //Eigen::Vector3d direction=-grad;
        //double w= -grad.dot(direction);
        double w= -grad.dot(direction);

        //std::cout<<w<<" direction:"<<direction.transpose()<<"\n";

        step=1.0;
        if(std::abs(direction(0))>0.5*M_PI || std::abs(direction(1))> 0.5*M_PI)
        {
          step=0.95*std::min(0.5*std::abs(M_PI/direction(0)),0.5*std::abs(M_PI/direction(1)));
        }
        //std::cout<<"step:"<<step<<"\n";
        tmp_c=current_c(c,  c0,  c1,
                          theta,  phi);
        tmp_d=d;
        tmp_theta=theta;
        tmp_phi=phi;
        e0=barrier_energy(position, _position, 
                                 tmp_c,  tmp_d);
        tmp_theta= theta+step*direction(0);    
        tmp_phi=phi+step*direction(1);
        tmp_c=current_c(c,  c0,  c1,
                       tmp_theta,  tmp_phi);  
        tmp_d=d+step*direction(2);                 
        e1=barrier_energy(position, _position, 
                                 tmp_c,  tmp_d);
         

        while(e0-1e-4*w*step< e1)
        {
          step*=0.8;
          
          tmp_theta= theta+step*direction(0);    
          tmp_phi=phi+step*direction(1);
          tmp_c=current_c(c,  c0,  c1,
                          tmp_theta,  tmp_phi);   
          tmp_d=d+step*direction(2);                       
          e1=barrier_energy(position, _position, 
                                 tmp_c,  tmp_d);
        }
        //std::cout<<kk<<" "<<grad.norm()<<" "<<step<<" "<<e1<<"\n";
        
        c=tmp_c;
        d=tmp_d;
        //std::cout<<hessian<<"\n\n"<<std::flush;
        //std::cout<<grad<<"\n\n"<<std::flush;
        //std::cout<<c.transpose()<< " "<<d<<"\n\n"<<std::flush;
        

      }
      //std::cout<<"\n";
    }
    */
    static double self_barrier_energy(const Data& position, const Data & _position, 
                         Eigen::Vector3d& c, double& d)
    {
      double energy=0;
      double dist;
        for(int j=0;j<position.rows();j++)
        {
            //d=P[j].dot(c_list[k])+d_list[k];
            dist=position.row(j).dot(c)+d-0.5*offset;
            if(dist<=0)
            {
                return INFINITY;
            }
            if(dist<margin)
            { 
                energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 

            }
        }
        for(int j=0;j<_position.rows();j++)
        {
            dist=-_position.row(j).dot(c)-d-0.5*offset;
            if(dist<=0)
            {
                return INFINITY;
            }
            if(dist<margin)
            { 
                energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 
               
            }
        }
        return energy;

    }

    static void self_barrier_grad(const Data& position, const Data & _position, 
                                  Eigen::Vector3d& c, Eigen::Vector3d& c0, Eigen::Vector3d& c1, double& d,
                                  Eigen::Vector3d& grad, Eigen::Matrix3d& hessian)
    {
        double dist;
        Eigen::Matrix3d tmp_h;
        double p_c,p_c0,p_c1;
        for(int j=0;j<position.rows();j++)
        {
            //d=P[j].dot(c_list[k])+d_list[k];
            dist=position.row(j).dot(c)+d-0.5*offset;
            if(dist<margin)
            { 
                //energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 
                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);  
                p_c=position.row(j).dot(c);
                p_c0=position.row(j).dot(c0);
                p_c1=position.row(j).dot(c1);

                double e1=-(2*(dist-margin)*log(dist/margin)+(dist-margin)*(dist-margin)/dist);

                double e2=-(2*log(dist/margin)+4*(dist-margin)/dist-(dist-margin)*(dist-margin)/(dist*dist)); 
                
                grad(0)+=e1*p_c0;
                grad(1)+=0;//e1*p_c1;
                grad(2)+=e1;

                tmp_h<<e2*p_c0*p_c0-e1*p_c,   e1*p_c1,  e2*p_c0,
                          e1*p_c1,            0,         0,
                          e2*p_c0,            0,        e2;
                hessian+=tmp_h;
                
               
            }
        }
        for(int j=0;j<_position.rows();j++)
        {
            dist=-_position.row(j).dot(c)-d-0.5*offset;
            if(dist<margin)
            { 
                //energy+=-(dist-margin)*(dist-margin)*log(dist/margin); 
                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);  

                p_c=-_position.row(j).dot(c);
                p_c0=-_position.row(j).dot(c0);
                p_c1=-_position.row(j).dot(c1); 

                double e1=-(2*(dist-margin)*log(dist/margin)+(dist-margin)*(dist-margin)/dist);

                double e2=-(2*log(dist/margin)+4*(dist-margin)/dist-(dist-margin)*(dist-margin)/(dist*dist)); 
                
                grad(0)+=e1*p_c0;
                grad(1)+=0;//e1*p_c1;
                grad(2)+=-e1;
               
                tmp_h<<e2*p_c0*p_c0-e1*p_c,   e1*p_c1,  -e2*p_c0,
                          e1*p_c1,            0,         0,
                          -e2*p_c0,            0,        e2;
                hessian+=tmp_h;
                
             
            }
        }

    }

    static void self_optimal_cd(const Data& position, const Data & _position, 
                             Eigen::Vector3d& c, double& d) //order_num+1, order_num+1
    {      
      Eigen::Vector3d c0, c1, temp_c;

      Eigen::Matrix3d hessian;
      Eigen::Vector3d grad, tmp_g;

      double theta, phi; // [-Pi/2,Pi/2] [-Pi/2,Pi/2]
      
      double step;
      while(true)
      //for(int kk=0;kk<200;kk++)
      {

        c0(0)=c(1); c0(1)=-c(0); c0(2)=0;
        c0.normalize();

        c1=c0.cross(c);
        c1.normalize();

        theta=0;
        phi=0;

        hessian.setZero();
        grad.setZero();

        self_barrier_grad( position, _position, 
                           c,  c0,  c1,  d,
                           grad,  hessian);
        if(grad.norm()<1e-2)
          break;
        Eigen::LLT<Eigen::Matrix3d> solver; 
        Eigen::Matrix3d I; I.setIdentity();

        solver.compute(hessian);
        //std::cout<<hessian<<"\n";
        if(solver.info() == Eigen::NumericalIssue)//50
        {
          
          Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(hessian);
          Eigen::MatrixXd eigenvalue=eigensolver.eigenvalues();
          if(eigenvalue(0)<0)
          {
            hessian=hessian-eigenvalue(0)*I+1e-8*I;
          }
          solver.compute(hessian);

        }
       

        //std::cout<<hessian<<"\n\n";
        double e0,e1;
        double tmp_d, tmp_theta, tmp_phi;
        Eigen::Vector3d tmp_c;
        
        
        /*finite difference*/
        /*
        tmp_c=current_c(c,  c0,  c1,
                          theta,  phi);
        tmp_d=d;
        tmp_theta=theta;
        tmp_phi=phi;
        e0=barrier_energy(position, _position, 
                          tmp_c,  tmp_d);

        Eigen::Vector3d num_grad;
        
        tmp_theta= theta+1e-8;    
        tmp_phi=phi;
        tmp_c=current_c(c,  c0,  c1,
                        tmp_theta,  tmp_phi);  
        tmp_d=d;                 
        e1=barrier_energy(position, _position, 
                          tmp_c,  tmp_d);
        num_grad(0)=(e1-e0)/1e-8;

        tmp_theta = theta;    
        tmp_phi = phi-1e-8;
        tmp_c=current_c(c,  c0,  c1,
                        tmp_theta,  tmp_phi);  
        tmp_d=d;                 
        e1=barrier_energy(position, _position, 
                          tmp_c,  tmp_d);
        num_grad(1)=-(e1-e0)/1e-8;

        tmp_theta= theta;    
        tmp_phi=phi;
        tmp_c=current_c(c,  c0,  c1,
                        tmp_theta,  tmp_phi);  
        tmp_d=d+1e-8;                 
        e1=barrier_energy(position, _position, 
                          tmp_c,  tmp_d);
        num_grad(2)=(e1-e0)/1e-8;

        std::cout<<"g:"<<(grad).transpose()<<"\n";
        std::cout<<"ng:"<<(num_grad).transpose()<<"\n";
        std::cout<<"dg:"<<(num_grad-grad).transpose()<<"\n";
        */
        
        Eigen::Vector3d direction=-solver.solve(grad);
        //Eigen::Vector3d direction=-grad;
        //double w= -grad.dot(direction);
        double w= -grad.dot(direction);

        //std::cout<<w<<" direction:"<<direction.transpose()<<"\n";

        step=1.0;
        if(std::abs(direction(0))>0.5*M_PI || std::abs(direction(1))> 0.5*M_PI)
        {
          step=0.95*std::min(0.5*std::abs(M_PI/direction(0)),0.5*std::abs(M_PI/direction(1)));
        }
        //std::cout<<"step:"<<step<<"\n";
        tmp_c=current_c(c,  c0,  c1,
                          theta,  phi);
        tmp_d=d;
        tmp_theta=theta;
        tmp_phi=phi;
        e0=self_barrier_energy(position, _position, 
                                 tmp_c,  tmp_d);
        tmp_theta= theta+step*direction(0);    
        tmp_phi=phi+step*direction(1);
        tmp_c=current_c(c,  c0,  c1,
                       tmp_theta,  tmp_phi);  
        tmp_d=d+step*direction(2);                 
        e1=self_barrier_energy(position, _position, 
                                 tmp_c,  tmp_d);
         

        while(e0-1e-4*w*step< e1)
        {
          step*=0.8;
          
          tmp_theta= theta+step*direction(0);    
          tmp_phi=phi+step*direction(1);
          tmp_c=current_c(c,  c0,  c1,
                          tmp_theta,  tmp_phi);   
          tmp_d=d+step*direction(2);                       
          e1=self_barrier_energy(position, _position, 
                                 tmp_c,  tmp_d);
        }
        //std::cout<<kk<<" "<<grad.norm()<<" "<<step<<" "<<e1<<"\n";
        
        c=tmp_c;
        d=tmp_d;
        //std::cout<<hessian<<"\n\n"<<std::flush;
        //std::cout<<grad<<"\n\n"<<std::flush;
        //std::cout<<c.transpose()<< " "<<d<<"\n\n"<<std::flush;
        

      }
      //std::cout<<"\n";
    }

    
};

PRJ_END

#endif
