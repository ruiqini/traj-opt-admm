#ifndef SEPARATE_H
#define SEPARATE_H

#include "Utils/CCDUtils.h"
#include <stdio.h>
#include <Eigen/Sparse>
//#include "QPwarper.h"
/*
// computes the distance between two cubes in R^3 using double
// as input type and some internal EXACT floating point type
#include <CGAL/Polytope_distance_d.h>
#include <CGAL/Polytope_distance_d_traits_3.h>
#include <CGAL/Homogeneous.h>
#include <iostream>
#include <cassert>
#include <CGAL/Gmpzf.h>
#include <CGAL/MP_Float.h>
*/
extern "C" {
#include <openGJK/openGJK.h>
}


PRJ_BEGIN

class Separate
{
  public:
    typedef Eigen::MatrixXd Data;
    /*
    typedef CGAL::MP_Float ET;
    //typedef CGAL::Homogeneous<ET>                         K;
    typedef CGAL::Homogeneous<double>                         K;
    typedef K::Point_3                                        Point;
    // ... and the EXACT traits class based on the inexcat kernel
    typedef CGAL::Polytope_distance_d_traits_3<K, ET, double> Traits;
    //typedef CGAL::Polytope_distance_d<Traits>                 Polytope_distance;
    //typedef CGAL::Polytope_distance_d_traits_3<K> Traits;
    typedef CGAL::Polytope_distance_d<Traits>     Polytope_distance;
    */
    /*
    static void svm(const Data& position, const Data & _position, 
                          Eigen::Vector3d& c, double& d)
    {
        int row=_position.rows();
        Eigen::MatrixXd A(order_num+1+row, 4);
        
        A.block(0,0,order_num+1,3)=position;
        A.block(order_num+1,0,row,3)=-_position;
        A.col(3).setOnes();
        A.block(order_num+1,3,row,1)*=-1;

        Eigen::MatrixXd Q(4,4);
        Q.setIdentity();
        Q(3,3)=0;

        Eigen::VectorXd  lc(order_num+1+row);
        lc.setOnes();

        Eigen::VectorXd b(4);

        Eigen::VectorXd b_=b; b_.setZero();
        QPwarper::QP(Q, b_,
                     A, lc, b);
        c=b.head(3);
        d=b(3);
        
        //std::cout<<"---\n";
        //std::cout<<c<<"\n"<<d<<"\n\n";

        d+=1;

        double c_n=c.norm();
        c/=c_n;
        d/=c_n;

        

        d-=offset;

        std::cout<<"---\n";
        std::cout<<c<<"\n"<<d<<"\n\n";
    }
    */
   /*
    static bool cgal(const Data& position, const Data & _position, const double & distance,
                      Eigen::Vector3d& c, double& d)
    {
        // the cube [0,1]^3
        std::vector<Point> PList;
        PList.resize((order_num+1));
        for(int i=0;i<=order_num;i++)
        {
          PList[i]=Point(position(i,0),position(i,1),position(i,2));
        }
        Point* P=PList.data();
        // the cube [2,3]^3
        Point Q[3] = { Point(_position(0,0),_position(0,1),_position(0,2)), 
                      Point(_position(1,0),_position(1,1),_position(1,2)), 
                      Point(_position(2,0),_position(2,1),_position(2,2))};
        Polytope_distance pd(P, P+(order_num+1), Q, Q+3);
        assert (pd.is_valid());
        
        // get squared distance (2,2,2)-(1,1,1))^2 = 3
        double cgal_distance2=
        CGAL::to_double (pd.squared_distance_numerator()) /
        CGAL::to_double (pd.squared_distance_denominator()) ;
        //return cgal_distance2 <= d*d;
        if(cgal_distance2 > distance*distance)
           return false;

        // get points that realize the distance
        Eigen::Vector3d p, q;
        int iter=0;
        Polytope_distance::Coordinate_iterator  coord_it;
        //std::cout << "p:"; // homogeneous point from first cube, (1,1,1,1)
        for (coord_it = pd.realizing_point_p_coordinates_begin();
            coord_it != pd.realizing_point_p_coordinates_end();
            ++coord_it)
          {
            if(iter<3)
            {
              double coord=CGAL::to_double(*coord_it);
              p(iter)=coord;
            }
            iter++;
          }
        //std::cout << " " << *coord_it;
        //std::cout << std::endl;

        iter=0;
        //std::cout << "q:"; // homogeneous point from second cube, (2,2,2,1)
        for (coord_it = pd.realizing_point_q_coordinates_begin();
            coord_it != pd.realizing_point_q_coordinates_end();
            ++coord_it)
          {
            if(iter<3)
            {
              double coord=CGAL::to_double(*coord_it);
              q(iter)=coord;
            }
            iter++;

          }
          //std::cout << " " << *coord_it;
          //std::cout << std::endl;

          c=p-q;
          d=INFINITY;
          for(int i=0;i<3;i++)
          {
            double d0=-c.dot(_position.row(i));
            if(d>d0)
              d=d0;
          }

          double c_n=c.norm();
          c/=c_n;
          d/=c_n;

          d-=offset;
          //std::cout<<c<<"\n"<<d<<"\n\n";
          return true;
    }
    */
    static bool opengjk(const Data& position, const Data & _position, const double & distance,
                      Eigen::Vector3d& c, double& d)
    {
      struct simplex  s;
      // Number of vertices defining body 1 and body 2, respectively.          
      int    nvrtx1,
             nvrtx2;

      // Structures of body 1 and body 2, respectively.                        
      struct bd       bd1;
      struct bd       bd2;
      // Specify name of input files for body 1 and body 2, respectively.      

      // Pointers to vertices' coordinates of body 1 and body 2, respectively. 

      //double (**vrtx1) = NULL,
      //       (**vrtx2) = NULL;
      
      // For importing openGJK this is Step 2: adapt the data structure for the
      // two bodies that will be passed to the GJK procedure. 
      nvrtx1 = order_num+1;
      nvrtx2 = 3;

      double **vrtx1 = (double **)malloc(nvrtx1 * sizeof(double *));
      for (int i=0; i<nvrtx1; i++)
      {
        vrtx1[i] = (double *)malloc(3 * sizeof(double));

        for(int j=0; j<3; j++) 
        {
          vrtx1[i][j] = position(i,j);
        }
      }
        
        
      double **vrtx2 = (double **)malloc(nvrtx2 * sizeof(double *));
      for (int i=0; i<nvrtx2; i++)
      {
        vrtx2[i] = (double *)malloc(3 * sizeof(double));

        for(int j=0; j<3; j++) 
        {
          vrtx2[i][j] = _position(i,j);
        }
      }
        

      // Import coordinates of object 1. 
      bd1.coord = vrtx1;
      bd1.numpoints = nvrtx1;

      // Import coordinates of object 2. 
      bd2.coord = vrtx2;
      bd2.numpoints = nvrtx2;

      // Initialise simplex as empty 
      s.nvrtx = 0;

    #ifdef DEBUG
      // Verify input of body A. 
      for (int i = 0; i < bd1.numpoints; ++i) {
        printf ( "%.2f ", vrtx1[i][0]);
        printf ( "%.2f ", vrtx1[i][1]);
        printf ( "%.2f\n", bd1.coord[i][2]);
      }

      // Verify input of body B. 
      for (int i = 0; i < bd2.numpoints; ++i) {
        printf ( "%.2f ", bd2.coord[i][0]);
        printf ( "%.2f ", bd2.coord[i][1]);
        printf ( "%.2f\n", bd2.coord[i][2]);
      }
    #endif

      // For importing openGJK this is Step 3: invoke the GJK procedure. 
      // Compute squared distance using GJK algorithm. 
      double *c0;
      
      c0=gjk (bd1, bd2, &s);
      
      // Free memory 
      for (int i=0; i<bd1.numpoints; i++)
        free(bd1.coord[i]);
      free(bd1.coord);
      for (int i=0; i<bd2.numpoints; i++)
        free(bd2.coord[i]);
      free(bd2.coord);

      c(0)=c0[0];
      c(1)=c0[1];
      c(2)=c0[2];

      double c_n=c.norm();
      if(c_n>distance)
        return false;
      
      c/=c_n;

      d=INFINITY;
      for(int i=0;i<3;i++)
      {
        double d0=-c.dot(_position.row(i));
        if(d>d0)
          d=d0;
      }

      d-=offset;
      //std::cout<<c<<"\n"<<d<<"\n\n";
      return true;
      
      // Print distance between objects. 
      //printf ("Distance between bodies %i\n", intersection);

      
    }


    static bool selfgjk(const Data& position, const Data & _position, const double & distance,
                        Eigen::Vector3d& c, double& d0, double& d1)
    {
      struct simplex  s;
      // Number of vertices defining body 1 and body 2, respectively.          
      int    nvrtx1,
             nvrtx2;

      // Structures of body 1 and body 2, respectively.                        
      struct bd       bd1;
      struct bd       bd2;
      // Specify name of input files for body 1 and body 2, respectively.      

      // Pointers to vertices' coordinates of body 1 and body 2, respectively. 

      //double (**vrtx1) = NULL,
      //       (**vrtx2) = NULL;
      
      // For importing openGJK this is Step 2: adapt the data structure for the
      // two bodies that will be passed to the GJK procedure. 
      nvrtx1 = order_num+1;
      nvrtx2 = order_num+1;

      double **vrtx1 = (double **)malloc(nvrtx1 * sizeof(double *));
      for (int i=0; i<nvrtx1; i++)
      {
        vrtx1[i] = (double *)malloc(3 * sizeof(double));

        for(int j=0; j<3; j++) 
        {
          vrtx1[i][j] = position(i,j);
        }
      }
        
        
      double **vrtx2 = (double **)malloc(nvrtx2 * sizeof(double *));
      for (int i=0; i<nvrtx2; i++)
      {
        vrtx2[i] = (double *)malloc(3 * sizeof(double));

        for(int j=0; j<3; j++) 
        {
          vrtx2[i][j] = _position(i,j);
        }
      }
        

      // Import coordinates of object 1. 
      bd1.coord = vrtx1;
      bd1.numpoints = nvrtx1;

      // Import coordinates of object 2. 
      bd2.coord = vrtx2;
      bd2.numpoints = nvrtx2;

      // Initialise simplex as empty 
      s.nvrtx = 0;

    #ifdef DEBUG
      // Verify input of body A. 
      for (int i = 0; i < bd1.numpoints; ++i) {
        printf ( "%.2f ", vrtx1[i][0]);
        printf ( "%.2f ", vrtx1[i][1]);
        printf ( "%.2f\n", bd1.coord[i][2]);
      }

      // Verify input of body B. 
      for (int i = 0; i < bd2.numpoints; ++i) {
        printf ( "%.2f ", bd2.coord[i][0]);
        printf ( "%.2f ", bd2.coord[i][1]);
        printf ( "%.2f\n", bd2.coord[i][2]);
      }
    #endif

      // For importing openGJK this is Step 3: invoke the GJK procedure. 
      // Compute squared distance using GJK algorithm. 
      double *c0;
      
      c0=gjk (bd1, bd2, &s);
      
      // Free memory 
      for (int i=0; i<bd1.numpoints; i++)
        free(bd1.coord[i]);
      free(bd1.coord);
      for (int i=0; i<bd2.numpoints; i++)
        free(bd2.coord[i]);
      free(bd2.coord);

      c(0)=c0[0];
      c(1)=c0[1];
      c(2)=c0[2];

      double c_n=c.norm();
      if(c_n>distance)
        return false;
      
      c/=c_n;

      d0=INFINITY;
      for(int i=0;i<order_num+1;i++)
      {
        double d_=-c.dot(_position.row(i));
        if(d0>d_)
          d0=d_;
      }

      //d0-=offset;

      d1=-INFINITY;
      for(int i=0;i<order_num+1;i++)
      {
        double d_=-c.dot(position.row(i));
        if(d1<d_)
          d1=d_;
      }

      //d1+=offset;

      //std::cout<<c<<"\n"<<d<<"\n\n";
      return true;
      
      // Print distance between objects. 
      //printf ("Distance between bodies %i\n", intersection);

      
    }

    static bool trigjk(const Data& position, const Data & _position, const double & distance,
                      Eigen::Vector3d& c, double& d0, double& d1)
    {
      struct simplex  s;
      // Number of vertices defining body 1 and body 2, respectively.          
      int    nvrtx1,
             nvrtx2;

      // Structures of body 1 and body 2, respectively.                        
      struct bd       bd1;
      struct bd       bd2;
      // Specify name of input files for body 1 and body 2, respectively.      

      // Pointers to vertices' coordinates of body 1 and body 2, respectively. 

      //double (**vrtx1) = NULL,
      //       (**vrtx2) = NULL;
      
      // For importing openGJK this is Step 2: adapt the data structure for the
      // two bodies that will be passed to the GJK procedure. 
      nvrtx1 = 3;
      nvrtx2 = 3;

      double **vrtx1 = (double **)malloc(nvrtx1 * sizeof(double *));
      for (int i=0; i<nvrtx1; i++)
      {
        vrtx1[i] = (double *)malloc(3 * sizeof(double));

        for(int j=0; j<3; j++) 
        {
          vrtx1[i][j] = position(i,j);
        }
      }
        
        
      double **vrtx2 = (double **)malloc(nvrtx2 * sizeof(double *));
      for (int i=0; i<nvrtx2; i++)
      {
        vrtx2[i] = (double *)malloc(3 * sizeof(double));

        for(int j=0; j<3; j++) 
        {
          vrtx2[i][j] = _position(i,j);
        }
      }
        

      // Import coordinates of object 1. 
      bd1.coord = vrtx1;
      bd1.numpoints = nvrtx1;

      // Import coordinates of object 2. 
      bd2.coord = vrtx2;
      bd2.numpoints = nvrtx2;

      // Initialise simplex as empty 
      s.nvrtx = 0;

    #ifdef DEBUG
      // Verify input of body A. 
      for (int i = 0; i < bd1.numpoints; ++i) {
        printf ( "%.2f ", vrtx1[i][0]);
        printf ( "%.2f ", vrtx1[i][1]);
        printf ( "%.2f\n", bd1.coord[i][2]);
      }

      // Verify input of body B. 
      for (int i = 0; i < bd2.numpoints; ++i) {
        printf ( "%.2f ", bd2.coord[i][0]);
        printf ( "%.2f ", bd2.coord[i][1]);
        printf ( "%.2f\n", bd2.coord[i][2]);
      }
    #endif

      // For importing openGJK this is Step 3: invoke the GJK procedure. 
      // Compute squared distance using GJK algorithm. 
      double *c0;
      
      c0=gjk (bd1, bd2, &s);
      
      // Free memory 
      for (int i=0; i<bd1.numpoints; i++)
        free(bd1.coord[i]);
      free(bd1.coord);
      for (int i=0; i<bd2.numpoints; i++)
        free(bd2.coord[i]);
      free(bd2.coord);

      c(0)=c0[0];
      c(1)=c0[1];
      c(2)=c0[2];

      double c_n=c.norm();
      if(c_n>distance)
        return false;
      
      c/=c_n;

      d0=INFINITY;
      for(int i=0;i<order_num+1;i++)
      {
        double d_=-c.dot(_position.row(i));
        if(d0>d_)
          d0=d_;
      }


      d1=-INFINITY;
      for(int i=0;i<order_num+1;i++)
      {
        double d_=-c.dot(position.row(i));
        if(d1<d_)
          d1=d_;
      }

      //std::cout<<c<<"\n"<<d<<"\n\n";
      return true;
      
      // Print distance between objects. 
      //printf ("Distance between bodies %i\n", intersection);

      
    }
};

PRJ_END

#endif
