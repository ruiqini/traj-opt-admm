#include <igl/opengl/glfw/Viewer.h>

//#include <igl/read_triangle_mesh.h>
//#include <igl/write_triangle_mesh.h>

#include "HighOrderCCD/Utils/CCDUtils.h"
#include "HighOrderCCD/OMPL/OMPL.h"

#include "HighOrderCCD/Optimization/Optimization3D_admm.h"
#include "HighOrderCCD/Optimization/Optimization3D_am.h"

#include "lib/nlohmann/json.hpp"    // https://github.com/nlohmann/json/tree/develop/single_include/nlohmann/json.hpp

USE_PRJ_NAMESPACE

using json = nlohmann::json;

typedef Eigen::MatrixXd Data;

Eigen::Vector3d getPosFromBezier(const Eigen::MatrixXd & polyCoeff, double t_now, int seg_now , int _traj_order) 
{
    Eigen::Vector3d ret = Eigen::VectorXd::Zero(3);
    Eigen::VectorXd ctrl_now = polyCoeff.row(seg_now);
    int ctrl_num1D = polyCoeff.cols() / 3;

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < ctrl_num1D; j++)
            ret(i) += combination[_traj_order][j] * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) * pow((1 - t_now), (_traj_order - j) ); 

    return ret;  
}

void log_data(std::string meshfile, Eigen::MatrixXd spline, double piece_time)
{
              Eigen::VectorXd my_time(piece_num);
              Eigen::MatrixXd coeff(piece_num, 3*(order_num+1));
              for(int i=0;i<piece_num;i++)
              {
                my_time(i)=time_weight[i]*piece_time;
                Eigen::MatrixXd bz;
                bz=spline.block<order_num+1,3>(i*(order_num-2),0);
                
                bz=convert_list[i]*bz;
                coeff.block<1,order_num+1>(i,0)=bz.col(0).transpose();
                coeff.block<1,order_num+1>(i,order_num+1)=bz.col(1).transpose();
                coeff.block<1,order_num+1>(i,2*(order_num+1))=bz.col(2).transpose();
              }
              /*
              std::ofstream file(meshfile+"_ccd_test.txt");
              if (file.is_open())
              {
                file << "coeff:\n" << coeff << '\n';
                file << "time:" << '\n' <<  my_time << '\n';
              }
              file.close();
              */
              std::vector<Eigen::Vector3d> ccd_traj;
              
             for (double t = 0.0; t < piece_num; t += 0.05 / piece_time){
                      int i = floor(t);
                      double cur_t=t-i;
                      Eigen::Vector3d state = getPosFromBezier( coeff, cur_t, i ,order_num);
                      Eigen::Vector3d cur;
                      cur(0) =   state(0);
                      cur(1) =   state(1);
                      cur(2) =   state(2);
                      ccd_traj.push_back(cur);
                  }
              std::cout<<"ccd time:"<<my_time.sum()<<std::endl;
              double  len_ccd=0;
              for(int i=0;i<(int)ccd_traj.size()-1;i++)
              {
                  len_ccd+=(ccd_traj[i+1]-ccd_traj[i]).norm();;
              }
              std::cout<<"ccd len:"<<len_ccd<<std::endl;
              
              std::ofstream curve_file;
              curve_file.open (meshfile + "_curve_file.txt");
              for(int i=0;i<ccd_traj.size();i++)
              {
                curve_file<<ccd_traj[i].transpose()<<"\n";
              }

              Eigen::MatrixXd v_ccd(2*ccd_traj.size()-1,3);
              Eigen::MatrixXi f_ccd(ccd_traj.size()-1,3);
              for(int i=0;i<(int)ccd_traj.size();i++)
              {
                  v_ccd.row(i)=ccd_traj[i].transpose();
              }
              for(int i=0;i<(int)ccd_traj.size()-1;i++)
              {
                  v_ccd.row(i+ccd_traj.size())=0.5*(ccd_traj[i].transpose()+ccd_traj[i+1].transpose());
              }
              for(int i=0;i<(int)ccd_traj.size()-1;i++)
              {
                  f_ccd(i,0)=i; f_ccd(i,1)=i+1; f_ccd(i,2)=i+ccd_traj.size();
              }
              
              

              //igl::write_triangle_mesh("ccd_traj_"+meshfile+".obj",v_ccd,f_ccd);
}

void way_point_init(const std::string& mesh_file, std::vector<Eigen::Vector3d>& way_points)
{
    std::string line;
    std::ifstream myfile ("init/" +mesh_file + "_init_file.txt");
    
    //V*=0.5;
    Eigen::Vector3d p0;
    std::vector<Eigen::Vector3d> read_points;
    if (myfile.is_open())
    {
      int i=0;
      while ( getline (myfile,line) )
      {
        std::istringstream iss(line);
        
        iss>>p0(0)>>p0(1)>>p0(2);
        read_points.push_back(p0);
        //result_file<<2*p0.transpose()<<"\n";
        i++;
      }
      myfile.close();
      //piece_num=i-1;
    }

    way_points=read_points;
    /*
    std::ofstream ff;
    ff.open("init/b2.obj_init_file.txt");
    for(int i=0;i<way_points.size();i++)
    {
      ff<<2*way_points[i].transpose()<<"\n";
    }
    */
}

bool edge_collision(const Eigen::MatrixXd& V,BVH& bvh, const std::vector<std::vector<Eigen::MatrixXd>>& edges, 
                    const Eigen::MatrixXd& edge)
{
           std::vector<unsigned int> collision_pair;
           bvh.EdgeCollision(edge, collision_pair,offset+0.5*margin);

            
           int collision_size=collision_pair.size();
           bool is_collided=false;
          for(int i=0;i<collision_size;i++)
          {
              //int ob_id=*it;
              int ob_id=collision_pair[i];

              Eigen::RowVector3d _position=V.row(ob_id);

              is_collided= CCD::GJKDCD(edge,_position, offset+0.5*margin);  //cgal
              if(is_collided)
              {
                return true;
              }     
              
          }

          for(int i=0;i<(int)edges.size();i++)
          {
            for(int j=0;j<(int)edges[i].size();j++)
            {
              is_collided= CCD::GJKDCD(edge,edges[i][j], offset+0.5*margin);  //cgal
              if(is_collided)
              {
                  return true;
              }
            }
          }

          
          return false;
}

void simplify_path(const Eigen::MatrixXd& V,BVH& bvh, const std::vector<std::vector<Eigen::MatrixXd>>& edges,
                   const std::vector<Eigen::Vector3d> & path, std::vector<Eigen::Vector3d>& tmp_path)
{
     std::vector<int> id_vec; id_vec.resize(path.size());
     std::vector<bool> rm_vec; rm_vec.resize(path.size());
     for(int i=0;i<(int)id_vec.size();i++)
     {
        id_vec[i]=i;
        rm_vec[i]=false;
     }
     int prev=0;
     int next=2;
     bool is_collided;
     for(int i=1;i<(int)id_vec.size()-1;i++)
     {
        Eigen::MatrixXd edge;
        edge.resize(2,3);
        edge.row(0)=path[prev].transpose();
        edge.row(1)=path[next].transpose();
        is_collided=edge_collision( V, bvh, edges, 
                                    edge);
        if(is_collided)
        {
          prev=i;
          next+=1;
        }
        else
        {
          next+=1;
          rm_vec[i]=true;
        }
     }

     tmp_path.clear();
     for(int i=0;i<(int)id_vec.size();i++)
     {
       if(rm_vec[i]==false)
         tmp_path.push_back(path[i]);
     }

}

void ompl_init(const Eigen::MatrixXd& V,BVH& bvh, std::vector<Eigen::Vector3d>& way_points)
{
      Eigen::VectorXd minV=V.colwise().minCoeff().transpose();
      Eigen::VectorXd maxV=V.colwise().maxCoeff().transpose();
      std::cout<<minV<<"\n";
      std::cout<<maxV<<"\n";
      Eigen::VectorXd lowerBound, upperBound;
      lowerBound=1.2*minV;
      upperBound=1.2*maxV;

     
      
      std::vector<Eigen::Vector3d> path;
      std::vector<Eigen::Vector3d> tmp_path;

      Eigen::Vector3d start,end;
      /*
      start(0)= 4;
      start(1)= 8;
      start(2)= 3.3;

      end(0)= -6;
      end(1)= -5.5;
      end(2)= 4.3;
      */

      start(0)= 2.7;
      start(1)= 0;
      start(2)= 0;

      end(0)= -2.7;
      end(1)= 0;
      end(2)= 0;

      std::vector<std::vector<Eigen::MatrixXd>> edges;
      edges.clear();
      OMPL ompl(lowerBound, upperBound,
                V,edges, bvh);
      
      if(ompl.planRRT(start, end, V, edges, bvh))
      {
        ompl.getPath(path);
      }

      simplify_path( V, bvh,  edges,
                      path,  tmp_path);
      std::cout<<"size:"<<path.size()<<" "<<tmp_path.size()<<"\n";
      way_points=tmp_path;
    
    
  
}

void init_variable(const std::vector<Eigen::Vector3d>& way_points,
                   const std::vector<Eigen::RowVector3d>& vertex_list, 
                   Data& spline, Data& p_slack, Data& p_lambda, 
                   double piece_time, Eigen::VectorXd& t_slack, Eigen::VectorXd& t_lambda)
{

    trajectory_num = (order_num+1)+(piece_num-1)*(order_num+1-3);
    spline.resize(trajectory_num,3);

    spline.row(0)=way_points[0].transpose();
    for(int i=0;i<piece_num;i++)
    { 
      
      Eigen::Vector3d head=0.9*way_points[i]+0.1*way_points[i+1];
      Eigen::Vector3d tail=0.9*way_points[i+1]+0.1*way_points[i];
      spline.row(i*(order_num-2)+1)=way_points[i].transpose();
      for(int j=1;j<order_num-2;j++)
      {
        spline.row(j+i*(order_num-2)+1)=double(order_num-3-j)/(order_num-4)*head.transpose()
                                         +(double)(j-1)/(order_num-4)*tail.transpose();
      }
      spline.row((i+1)*(order_num-2)+1)=way_points[i+1].transpose();
    }
    spline.row(trajectory_num-1)=way_points[piece_num].transpose();

  spline.row(1)=spline.row(0);
  spline.row(trajectory_num-2)=spline.row(trajectory_num-1);
  
  
  p_lambda.resize((order_num+1)*piece_num,3); p_lambda.setZero();
  p_slack.resize((order_num+1)*piece_num,3);
  for(int sp_id=0; sp_id<piece_num; sp_id++)
  {
    p_slack.block<order_num+1,3>(sp_id*(order_num+1),0)=convert_list[sp_id]*spline.block<order_num+1,3>(sp_id*(order_num-2),0);
  }

  
  t_lambda.resize(piece_num); t_lambda.setZero();
  t_slack.resize(piece_num);
  for(int sp_id=0; sp_id<piece_num; sp_id++)
  {
    t_slack(sp_id)=piece_time;
  }

  M_dynamic=Dynamic3D<order_num, der_num>::dynamic_matrix();//Dynamic<order_num,1>::dynamic_matrix()+
  
  subdivide_tree.resize(piece_num*res);
  A_list.resize(piece_num*res);
  A_vel_list.resize(piece_num*res);
  A_acc_list.resize(piece_num*res);

  Eigen::MatrixXd basis, tmp_basis;
  
  Eigen::Matrix3d I; I.setIdentity();
  for(int k=0;k<res;k++)
  {
    double a=k/double(res),b=(k+1)/double(res);

    Blossom<order_num>::coefficient(basis, a, b);

    for(int i=0;i<piece_num;i++)
    { 
      std::pair<double,double> range(a,b);
      //basis*=;
      subdivide_tree[i*res+k]=std::make_tuple(i,range,basis*convert_list[i]);
      tmp_basis=basis*convert_list[i];
      
      A_list[i*res+k].resize(order_num+1);
      A_vel_list[i*res+k].resize(order_num);
      A_acc_list[i*res+k].resize(order_num-1);
      
      for(int j=0;j<=order_num;j++)
      {
        Eigen::MatrixXd A=Eigen::kroneckerProduct(tmp_basis.row(j),I);
        A.transposeInPlace();
        A_list[i*res+k][j]=A;
        if(j<order_num)
        {
          A=Eigen::kroneckerProduct(tmp_basis.row(j+1),I)-
            Eigen::kroneckerProduct(tmp_basis.row(j),I);
          A_vel_list[i*res+k][j]=A;
        }
        if(j<order_num-1)
        {
          A=Eigen::kroneckerProduct(tmp_basis.row(j+2),I)-
            2 * Eigen::kroneckerProduct(tmp_basis.row(j+1),I)+
            Eigen::kroneckerProduct(tmp_basis.row(j),I);
          A_acc_list[i*res+k][j]=A;
        }
        
      }
    }
  }  

  is_seperate.resize(piece_num*res);
  seperate_c.resize(piece_num*res);
  seperate_d.resize(piece_num*res);
  for(int i=0;i<piece_num*res;i++)
  {
    is_seperate[i].resize(vertex_list.size());
    seperate_c[i].resize(vertex_list.size());
    seperate_d[i].resize(vertex_list.size());
  }

}

std::vector<Eigen::Vector3d> next_point(Eigen::Vector3d cur_point, Eigen::Vector3d goal_point, 
                                        std::vector<std::pair<Eigen::Vector3d,Eigen::Vector3d>> exist_paths)
{
    std::vector<Eigen::Vector3d> point_list;
    for(int i=0;i<3;i++)
    {
        if(cur_point(i)!=goal_point(i))
        {
            Eigen::Vector3d point=cur_point;
            if(goal_point(i)-cur_point(i)>0)
              point(i)+=1;
            else
              point(i)+=-1;
            bool valid=true;
            for(int j=0;j<exist_paths.size();j++)
            {
                if(point==exist_paths[j].second)
                   valid=false;
                if(point==exist_paths[j].first && cur_point==exist_paths[j].second)
                   valid=false; 
            }
            if(valid==true)
              point_list.push_back(point);
        }
    }
    return point_list;
}

std::vector<std::vector<Eigen::Vector3d>> next_state(std::vector<Eigen::Vector3d> cur, std::vector<Eigen::Vector3d> goal)
{
    Eigen::Vector3d cur_point; Eigen::Vector3d goal_point; 
    std::vector<std::pair<Eigen::Vector3d,Eigen::Vector3d>> exist_paths;
  
    exist_paths.clear();
    std::pair<Eigen::Vector3d,Eigen::Vector3d> tmp;
    tmp.first=Eigen::Vector3d(-10,-10,-10);
    tmp.second=Eigen::Vector3d(-10,-10,-10);
    exist_paths.resize(4,tmp);
  
    cur_point=cur[0];
    goal_point=goal[0];
    std::vector<Eigen::Vector3d> j0_points = next_point(cur_point, goal_point, 
                             exist_paths);//

    std::vector<std::vector<Eigen::Vector3d>> all_path;
    std::vector<Eigen::Vector3d> path;
    path.resize(4);
    
    for(int j0=0;j0<j0_points.size();j0++)
    {
        path[0]=j0_points[j0];
        std::cout<<"j0:"<<j0_points[j0].transpose()<<"\n";
        std::pair<Eigen::Vector3d,Eigen::Vector3d> pair_path;
        pair_path.first=cur_point;
        pair_path.second=j0_points[j0];
        exist_paths[0]=(pair_path);

        cur_point=cur[1];
        goal_point=goal[1];
        std::vector<Eigen::Vector3d> j1_points = next_point(cur_point, goal_point, 
                                 exist_paths);//
        for(int j1=0;j1<j1_points.size();j1++)
        {
            path[1]=j1_points[j1];
            std::cout<<"j1:"<<j1_points[j1].transpose()<<"\n";
            std::pair<Eigen::Vector3d,Eigen::Vector3d> pair_path;
            pair_path.first=cur_point;
            pair_path.second=j1_points[j1];
            exist_paths[1]=(pair_path);

            cur_point=cur[2];
            goal_point=goal[2];
            std::vector<Eigen::Vector3d> j2_points = next_point(cur_point, goal_point, 
                                                                exist_paths);//
            for(int j2=0;j2<j2_points.size();j2++)
            {
                path[2]=j2_points[j2];
                std::cout<<"j2:"<<j2_points[j2].transpose()<<"\n";
                std::pair<Eigen::Vector3d,Eigen::Vector3d> pair_path;
                pair_path.first=cur_point;
                pair_path.second=j2_points[j2];
                exist_paths[2]=(pair_path);
                
                cur_point=cur[3];
                goal_point=goal[3];
                std::vector<Eigen::Vector3d> j3_points = next_point(cur_point, goal_point, 
                                        exist_paths);//
                for(int j3=0;j3<j3_points.size();j3++)
                {
                    path[3]=j3_points[j3];
                    std::cout<<"j3:"<<j3_points[j3].transpose()<<"\n";
                    all_path.push_back(path);
                }
            }
        }


    }
    for(int i=0;i<all_path.size();i++)
    {
        for(int j=0;j<4;j++)
        {
            std::cout<<all_path[i][j].transpose()<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"end\n";
    return all_path;
}


int main(int argc, char *argv[])
{
  std::vector<Eigen::Vector3d> start, goal, cur;
  start.resize(4);
  goal.resize(4);
  
  start[0]=Eigen::Vector3d(-1,-1,-1);
  goal[0] =Eigen::Vector3d(1,1,1);
  
  start[1]=Eigen::Vector3d(-1,1,-1);
  goal[1] =Eigen::Vector3d(1,-1,1);
  
  start[2]=Eigen::Vector3d(-1,-1,1);
  goal[2] =Eigen::Vector3d(1,1,-1);

  start[3]=Eigen::Vector3d(-1,1,1);
  goal[3] =Eigen::Vector3d(1,-1,-1);

  cur=start;
  
  //std::vector<std::vector<Eigen::Vector3d>> all_paths;
  //std::vector<Eigen::Vector3d> path;
  Eigen::Vector3d cur_point; Eigen::Vector3d goal_point; 
  std::vector<std::pair<Eigen::Vector3d,Eigen::Vector3d>> exist_paths;
  std::vector<Eigen::Vector3d> next_points;


  std::vector<std::vector<Eigen::Vector3d>> all_path0;
  
  std::vector<std::vector<std::vector<Eigen::Vector3d>>> allpath;
  std::vector<std::vector<Eigen::Vector3d>> path;
  path.resize(6);
  all_path0 = next_state(cur, goal);

  for(int i0=0;i0<all_path0.size();i0++)
  {
      
        cur=all_path0[i0];
        path[0]=cur;
        std::vector<std::vector<Eigen::Vector3d>> all_path1= next_state(cur, goal);
        for(int i1=0;i1<all_path1.size();i1++)
        {
            cur=all_path1[i1];
            path[1]=cur;
            std::vector<std::vector<Eigen::Vector3d>> all_path2= next_state(cur, goal);
            for(int i2=0;i2<all_path2.size();i2++)
            {
                cur=all_path2[i2];
                path[2]=cur;
                std::vector<std::vector<Eigen::Vector3d>> all_path3= next_state(cur, goal);
                for(int i3=0;i3<all_path3.size();i3++)
                {
                    cur=all_path3[i3];
                    path[3]=cur;
                    std::vector<std::vector<Eigen::Vector3d>> all_path4= next_state(cur, goal);
                    for(int i4=0;i4<all_path4.size();i4++)
                    {
                        cur=all_path4[i4];
                        path[4]=cur;
                        std::vector<std::vector<Eigen::Vector3d>> all_path5= next_state(cur, goal);
                        for(int i5=0;i5<all_path5.size();i5++)
                        {
                            cur=all_path5[i5];
                            path[5]=cur;

                            allpath.push_back(path);
                            //std::vector<std::vector<Eigen::Vector3d>> all_path5= next_state(cur, goal);
                        }
                    }
                }
            }
        }
        
  }
    std::cout<<allpath.size()<<"\n";
    iter=0;
    for(int k=0;k<60000;k=k+100)
    {
        std::ofstream myfile("multi_init/init_file"+std::to_string(iter)+".txt");
        for(int j=0;j<4;j++)
        {
            myfile<<start[j].transpose()<<" ";
        }
        myfile<<"\n";
        for(int i=0;i<6;i++)
        {
            for(int j=0;j<4;j++)
            {
                myfile<<allpath[k][i][j].transpose()<<" ";
            }
            myfile<<"\n";
        }
        iter++;
    }
    
    

  return 0;
}

