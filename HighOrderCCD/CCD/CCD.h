#ifndef CCD_H
#define CCD_H

#include "HighOrderCCD/Utils/CCDUtils.h"

extern "C" {
    #include <openGJK/openGJK.h>
}

PRJ_BEGIN

class CCD 
{
  public:
    typedef Eigen::MatrixXd Data;
    
    static bool GJKDCD(const Data& position, const Data& _position ,const double& d) //TT
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
      nvrtx1 = position.rows();
      nvrtx2 = _position.rows();

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

     
      double gjk_distance2=c0[0]*c0[0]+c0[1]*c0[1]+c0[2]*c0[2];
      // Print distance between objects. 
      //printf ("Distance between bodies %i\n", intersection);
      return gjk_distance2 <= d*d;
    }

    static bool GJKCCD(const Data& position, const Data& direction, const Data& _position, 
                       const double& d, const double& tMin, const double& tMax)
    {
      Data position0=position+tMin*direction;
      Data position1=position+tMax*direction;

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
      nvrtx1 = 2*position.rows();
      nvrtx2 = _position.rows();

      double **vrtx1 = (double **)malloc(nvrtx1 * sizeof(double *));
      for (int i=0; i<order_num+1; i++)
      {
        vrtx1[i] = (double *)malloc(3 * sizeof(double));
        
        for(int j=0; j<3; j++) 
        {
          vrtx1[i][j] = position0(i,j);
        }
      }

      for (int i=0; i<order_num+1; i++)
      {
        vrtx1[i+order_num+1] = (double *)malloc(3 * sizeof(double));
        
        for(int j=0; j<3; j++) 
        {
          vrtx1[i+order_num+1][j] = position1(i,j);
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

     
      double gjk_distance2=c0[0]*c0[0]+c0[1]*c0[1]+c0[2]*c0[2];
      // Print distance between objects. 
      //printf ("Distance between bodies %i\n", intersection);
      return gjk_distance2 <= d*d;
    }

    static bool SelfGJKCCD(const Data& position, const Data& direction, 
                           const Data& _position, const Data& _direction , const double& d, 
                           const double& tMin, const double& tMax, const double& _tMin, const double& _tMax)
    {
      Data position0=position+tMin*direction;
      Data position1=position+tMax*direction;

      Data _position0=_position+_tMin*_direction;
      Data _position1=_position+_tMax*_direction;

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
      nvrtx1 = 2*position.rows();
      nvrtx2 = 2*_position.rows();

      double **vrtx1 = (double **)malloc(nvrtx1 * sizeof(double *));
      for (int i=0; i<order_num+1; i++)
      {
        vrtx1[i] = (double *)malloc(3 * sizeof(double));
        
        for(int j=0; j<3; j++) 
        {
          vrtx1[i][j] = position0(i,j);
        }
      }

      for (int i=0; i<order_num+1; i++)
      {
        vrtx1[i+order_num+1] = (double *)malloc(3 * sizeof(double));
        
        for(int j=0; j<3; j++) 
        {
          vrtx1[i+order_num+1][j] = position1(i,j);
        }
      }
        
        
      double **vrtx2 = (double **)malloc(nvrtx2 * sizeof(double *));
      for (int i=0; i<order_num+1; i++)
      {
        vrtx2[i] = (double *)malloc(3 * sizeof(double));
        
        for(int j=0; j<3; j++) 
        {
          vrtx2[i][j] = _position0(i,j);
        }
      }

      for (int i=0; i<order_num+1; i++)
      {
        vrtx2[i+order_num+1] = (double *)malloc(3 * sizeof(double));
        
        for(int j=0; j<3; j++) 
        {
          vrtx2[i+order_num+1][j] = _position1(i,j);
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

     
      double gjk_distance2=c0[0]*c0[0]+c0[1]*c0[1]+c0[2]*c0[2];
      // Print distance between objects. 
      //printf ("Distance between bodies %i\n", intersection);
      return gjk_distance2 <= d*d;
    }

    static bool KDOPDCD(const Data& position, const Data& _position, const double& d)
    {
      //Data A((order_num+1),3);
      //A<<position;
      const double *position_data=position.data();
      const double *_position_data=_position.data();
      double *kdop_data=kdop_matrix.data();
      int dim = kdop_axis.size();
      double level;
      double x,y,z;
      for(int k=0;k<dim;k++)
      {
        double upperA=-INFINITY;
        double lowerA=INFINITY;
        x=kdop_data[3*k];//kdop_axis[k](0);
        y=kdop_data[3*k+1];//kdop_axis[k](1);
        z=kdop_data[3*k+2];//kdop_axis[k](2);
        for(int i=0;i<(order_num+1);i++)
        {
          level = x*position_data[i] +
                  y*position_data[i+(order_num+1)]+
                  z*position_data[i+2*(order_num+1)];//kdop_axis[k].dot(A.row(i));
          if(level<lowerA)
            lowerA=level;
          if(level>upperA)
            upperA=level;
        }

        double upperB=-INFINITY;
        double lowerB=INFINITY;
        for(int i=0;i<3;i++)
        {
          level = x*_position_data[i] +
                  y*_position_data[i+3]+
                  z*_position_data[i+2*3];//kdop_axis[k].dot(_position.row(i));
          if(level<lowerB)
            lowerB=level;
          if(level>upperB)
            upperB=level;
        }
        if(upperB<lowerA-d || upperA<lowerB-d)
          return false;
      }
      return true;
      
      //return KDOPConvexCollision<double,3>::hasCollision(A,_position,d);

    }

    
    static bool KDOPCCD(const Data& position, const Data& direction, const Data& _position,  
                       const double& d, const double& tMin, const double& tMax)
    {
      Data A(2*(order_num+1),3);
      A<<position+tMin*direction,position+tMax*direction;
      double *A_data=A.data();
      const double *_position_data=_position.data();
      double *kdop_data=kdop_matrix.data();
      int dim = kdop_axis.size();
      double level;
      double x,y,z;
      for(int k=0;k<dim;k++)
      {
        double upperA=-INFINITY;
        double lowerA=INFINITY;
        x=kdop_data[3*k];//kdop_axis[k](0);
        y=kdop_data[3*k+1];//kdop_axis[k](1);
        z=kdop_data[3*k+2];//kdop_axis[k](2);
        for(int i=0;i<2*(order_num+1);i++)
        {
          level = x*A_data[i] +
                  y*A_data[i+2*(order_num+1)]+
                  z*A_data[i+4*(order_num+1)];//kdop_axis[k].dot(A.row(i));
          if(level<lowerA)
            lowerA=level;
          if(level>upperA)
            upperA=level;
        }

        double upperB=-INFINITY;
        double lowerB=INFINITY;
        for(int i=0;i<3;i++)
        {
          level = x*_position_data[i] +
                  y*_position_data[i+3]+
                  z*_position_data[i+2*3];//kdop_axis[k].dot(_position.row(i));
          if(level<lowerB)
            lowerB=level;
          if(level>upperB)
            upperB=level;
        }
        if(upperB<lowerA-d || upperA<lowerB-d)
          return false;
      }
      return true;
      
      //return KDOPConvexCollision<double,3>::hasCollision(A,_position,d);

    }

    static bool SelfKDOPCCD(const Data& position, const Data& direction, 
                            const Data& _position, const Data& _direction , const double& d, 
                            const double& tMin, const double& tMax, const double& _tMin, const double& _tMax)
    {
      Data A(2*(order_num+1),3);
      A<<position+tMin*direction,position+tMax*direction;

      Data B(2*(order_num+1),3);
      B<<_position+_tMin*_direction,_position+_tMax*_direction;
      
      int dim = kdop_axis.size();
      double level;
      for(int k=0;k<dim;k++)
      {
        double upperA=-INFINITY;
        double lowerA=INFINITY;
        for(int i=0;i<2*(order_num+1);i++)
        {
          level = kdop_axis[k].dot(A.row(i));
          if(level<lowerA)
            lowerA=level;
          if(level>upperA)
            upperA=level;
        }

        double upperB=-INFINITY;
        double lowerB=INFINITY;
        for(int i=0;i<2*(order_num+1);i++)
        {
          level = kdop_axis[k].dot(B.row(i));
          if(level<lowerB)
            lowerB=level;
          if(level>upperB)
            upperB=level;
        }
        if(upperB<lowerA-d || upperA<lowerB-d)
          return false;
      }
      return true;
      
      //return KDOPConvexCollision<double,3>::hasCollision(A,_position,d);

    }

    static bool SelfKDOPDCD(const Data& position,const Data& _position, const double& d)
    {
      
      int dim = kdop_axis.size();
      double level;
      for(int k=0;k<dim;k++)
      {
        
        double upperA=-INFINITY;
        double lowerA=INFINITY;
        for(int i=0;i<(order_num+1);i++)
        {
          level = kdop_axis[k].dot(position.row(i));
          if(level<lowerA)
            lowerA=level;
          if(level>upperA)
            upperA=level;
        }

        double upperB=-INFINITY;
        double lowerB=INFINITY;
        for(int i=0;i<(order_num+1);i++)
        {
          level = kdop_axis[k].dot(_position.row(i));
          if(level<lowerB)
            lowerB=level;
          if(level>upperB)
            upperB=level;
        }

        
        if(upperB<lowerA-d || upperA<lowerB-d)
          return false;
      }
      
      return true;
      
      //return KDOPConvexCollision<double,3>::hasCollision(A,_position,d);

    }
  
    static bool KDOPCCD_point(Data position, Data direction, Data _position,  double d, double tMin, double tMax)
    {
      Data A(2*(order_num+1),3);
      A<<position+tMin*direction,position+tMax*direction;
      
      int dim = kdop_axis.size();
      double level;
      for(int k=0;k<dim;k++)
      {
        double upperA=-INFINITY;
        double lowerA=INFINITY;
        for(int i=0;i<2*(order_num+1);i++)
        {
          level = kdop_axis[k].dot(A.row(i));
          if(level<lowerA)
            lowerA=level;
          if(level>upperA)
            upperA=level;
        }

        double upperB=-INFINITY;
        double lowerB=INFINITY;
    
        level = kdop_axis[k].dot(_position.row(0));
        lowerB=level;
        upperB=level;
        
        if(upperB<lowerA-d || upperA<lowerB-d)
          return false;
      }
      return true;
      
      //return KDOPConvexCollision<double,3>::hasCollision(A,_position,d);

    }

    static bool KDOPDCD_point(Data position, Data _position,  double d)
    {
      Data A((order_num+1),3);
      A<<position;
      
      int dim = kdop_axis.size();
      double level;
      for(int k=0;k<dim;k++)
      {
        double upperA=-INFINITY;
        double lowerA=INFINITY;
        for(int i=0;i<(order_num+1);i++)
        {
          level = kdop_axis[k].dot(A.row(i));
          if(level<lowerA)
            lowerA=level;
          if(level>upperA)
            upperA=level;
        }

        double upperB=-INFINITY;
        double lowerB=INFINITY;
     
        level = kdop_axis[k].dot(_position.row(0));
        lowerB=level;
        upperB=level;

        if(upperB<lowerA-d || upperA<lowerB-d)
          return false;
      }
      return true;
      
      //return KDOPConvexCollision<double,3>::hasCollision(A,_position,d);

    }

};

PRJ_END

#endif