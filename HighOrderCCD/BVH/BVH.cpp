
#include "BVH.h"

PRJ_BEGIN

BVH::BVH()
{
	
}
BVH::~BVH()
{
	
}

void BVH::InitObstacle(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    unsigned int n_ob=F.rows();
    int dim = aabb_axis.size();
    
    ob_tree=aabb::Tree(dim,0.0, n_ob,true);
    
    std::cout << "\nInserting ob particles into AABB tree ...\n";
    for (unsigned int i=0;i<n_ob;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }

        for(int j=0;j<3;j++)
        {
          for(int k=0;k<dim;k++)
          {
            double level = aabb_axis[k].dot(V.row(F(i,j)));
            if(level<lowerBound[k])
              lowerBound[k]=level;
            if(level>upperBound[k])
              upperBound[k]=level;
          }

        }


        ob_tree.insertParticle(i, lowerBound, upperBound);
    }

}

void BVH::InitPointcloud(const Eigen::MatrixXd& V)
{
    unsigned int n_pc=V.rows();
    int dim = aabb_axis.size();
    /*
    for(int k=0;k<dim;k++)
    {
      aabb_axis[k].normalize();
    }
    */
    pc_tree=aabb::Tree(dim,0.0, n_pc,true);
    
    std::cout << "\nInserting pc particles into AABB tree ...\n";
    for (unsigned int i=0;i<n_pc;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }

       
        for(int k=0;k<dim;k++)
        {
          double level = aabb_axis[k].dot(V.row(i));
          if(level<lowerBound[k])
            lowerBound[k]=level;
          if(level>upperBound[k])
            upperBound[k]=level;
        }

        


        pc_tree.insertParticle(i, lowerBound, upperBound);
    }

}

void BVH::EdgeCollision(const Data& edge, std::vector<unsigned int>& collision_pair, double margin)
{
    int dim = aabb_axis.size();
    std::vector<double> lowerBound(dim);
    std::vector<double> upperBound(dim);

       
    std::vector<Eigen::RowVector3d> P(2);

        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }

        for(int j=0;j<2;j++)
        {

          P[j]=edge.row(j);

          for(int k=0;k<dim;k++)
          {
            double level = aabb_axis[k].dot(P[j]);
            if(level<lowerBound[k])
              lowerBound[k]=level;
            if(level>upperBound[k])
              upperBound[k]=level;
          }
           
        }
        
        //tr_tree.insertParticle(i, lowerBound, upperBound);
        aabb::AABB convex=aabb::AABB(lowerBound, upperBound);
        collision_pair=ob_tree.query(convex, margin);

}

void BVH::DCDCollision(const Data& spline, std::vector<std::vector<unsigned int>>& collision_pairs, double margin)
{
    unsigned int n_tr=subdivide_tree.size();
    int dim = aabb_axis.size();
    for (unsigned int i=0;i<n_tr;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        int sp_id=std::get<0>(subdivide_tree[i]);
        Eigen::MatrixXd basis=std::get<2>(subdivide_tree[i]);
        std::vector<Eigen::RowVector3d> P(order_num+1);
        
        Eigen::MatrixXd bz;
        bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);

        
        
        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }

        for(int j=0;j<=order_num;j++)
        {

            P[j].setZero();
            for(int j0=0;j0<=order_num;j0++)
            {
              P[j]+=basis(j,j0)*bz.row(j0);
            } 

            for(int k=0;k<dim;k++)
            {
              double level = aabb_axis[k].dot(P[j]);
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;
            }
           
        }
        
        //tr_tree.insertParticle(i, lowerBound, upperBound);
        aabb::AABB convex=aabb::AABB(lowerBound, upperBound);
        collision_pairs.push_back(ob_tree.query(convex, margin));
    }
    
}

void BVH::CCDCollision(const Data& spline, const Data& direction, std::vector<std::vector<unsigned int>>& collision_pairs, double margin)
{
    unsigned int n_tr=subdivide_tree.size();
    int dim = aabb_axis.size();
    for (unsigned int i=0;i<n_tr;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        int sp_id=std::get<0>(subdivide_tree[i]);
        Eigen::MatrixXd basis=std::get<2>(subdivide_tree[i]);
        
        std::vector<Eigen::RowVector3d> P(order_num+1),D(order_num+1);
        
        Eigen::MatrixXd bz, bz_d;
        bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
        bz_d=direction.block<order_num+1,3>(sp_id*(order_num-2),0);
       
        
        
        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }
        for(int j=0;j<=order_num;j++)
        {   
            P[j].setZero();
            D[j].setZero();
            for(int j0=0;j0<=order_num;j0++)
            {
              P[j]+=basis(j,j0)*bz.row(j0);
              D[j]+=basis(j,j0)*bz_d.row(j0);
            } 
            for(int k=0;k<dim;k++)
            {
              double level = aabb_axis[k].dot(P[j]);
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;

              level = aabb_axis[k].dot(P[j]+D[j]);
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;
            }
        }
        
        //tr_tree.insertParticle(i, lowerBound, upperBound);
        aabb::AABB convex=aabb::AABB(lowerBound, upperBound);
        collision_pairs.push_back(ob_tree.query(convex, margin));
    }
}

void BVH::SelfDCDCollision(const std::vector<Data>& P, std::vector<id_pair>& collision_pair, double margin)
{
    aabb::Tree self_tree;
    unsigned int n_tr=P.size();
    int dim = aabb_axis.size();
    self_tree=aabb::Tree(dim, 0.0, n_tr,true);
 
    //std::cout << n_tr<<"\nInserting tr particles into AABB tree ...\n";
    for (unsigned int i=0;i<n_tr;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }

        for(int j=0;j<=order_num;j++)
        {
            for(int k=0;k<dim;k++)
            {
              double level = aabb_axis[k].dot(P[i].row(j));
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;
            }
           
        }
        self_tree.insertParticle(i, lowerBound, upperBound);
    }
    collision_pair=self_tree.query(margin);

}

void BVH::SelfCCDCollision(const std::vector<Data>& P, const std::vector<Data>& D, std::vector<id_pair>& collision_pair, double margin)
{
    aabb::Tree self_tree;
    unsigned int n_tr=P.size();
    int dim = aabb_axis.size();
    self_tree=aabb::Tree(dim, 0.0, n_tr,true);
 
    //std::cout << n_tr<<"\nInserting tr particles into AABB tree ...\n";
    for (unsigned int i=0;i<n_tr;i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        for(int k=0;k<dim;k++)
        {
          upperBound[k]=-INFINITY;
          lowerBound[k]=INFINITY;
        }

        for(int j=0;j<=order_num;j++)
        {
            for(int k=0;k<dim;k++)
            {
              double level = aabb_axis[k].dot(P[i].row(j));
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;

              level = aabb_axis[k].dot((P[i]+D[i]).row(j));
              if(level<lowerBound[k])
                lowerBound[k]=level;
              if(level>upperBound[k])
                upperBound[k]=level;
            }

        }
        self_tree.insertParticle(i, lowerBound, upperBound);
    }
    collision_pair=self_tree.query(margin);

}

PRJ_END

