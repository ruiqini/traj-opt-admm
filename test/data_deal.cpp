
#include <igl/opengl/glfw/Viewer.h>

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/png/readPNG.h>

#include "HighOrderCCD/Utils/CCDUtils.h"

USE_PRJ_NAMESPACE

int main(int argc, char *argv[])
{
  

  

  typedef Eigen::MatrixXd Data;

  
  igl::opengl::glfw::Viewer viewer;
 
  viewer.core().background_color<< 1.0f, 1.0f, 1.0f, 1.0f;
  viewer.core().is_animating = true;
  viewer.core().camera_zoom = 5.0f;
  viewer.data().line_width = 1.0f;
  viewer.data().point_size = 3.0f;

  Eigen::MatrixXd V, BV;
  Eigen::MatrixXi F, BF;
  //std::ifstream myfile;
    std::string line;
    std::vector<std::vector<Eigen::Vector3d> > point_lists;
    int uav_num=4;
    point_lists.resize(uav_num);
    //std::string name("/home/n/ADMM/draw_curve/multi_traj/12_multi_data/"+std::to_string(0)+".txt");
     // std::cout<<name<<" ";
    
    for(int id=0;id<uav_num;id++)
    {
      
      //std::string name("/home/n/ADMM/draw_curve/multi_traj/12_new/"+std::to_string(id)+".txt");
      //std::string name("/home/n/ADMM/point_cloud/NLC_ADMM-UAV/build/multi/empty.obj_curve_file"+std::to_string(id)+".txt");

      //std::string name("/home/n/ADMM/draw_curve/multi_traj/4_new/"+std::to_string(id)+".txt");
      std::string name("/home/n/ADMM/point_cloud/NLC_ADMM-UAV/build/multi/multi_uav.obj_curve_file"+std::to_string(id)+".txt");
      //std::cout<<name<<" ";
      std::ifstream myfile=std::ifstream(name);
      point_lists[id].clear();
      while ( getline (myfile,line) )
      {
        std::istringstream iss(line);
        
        Eigen::Vector3d p0;
        iss>>p0(0)>>p0(1)>>p0(2);
        //p0(2)-=2;
        //p0(2)-=1;
        point_lists[id].push_back(p0);
      }
      myfile.close(); 

    }
      
    viewer.data().line_width = 5.0f;
      viewer.data().point_size = 30.0f;
      
      /*
      Eigen::RowVector3d C0(1.0,0.1,0.1);
      Eigen::RowVector3d C1(0.1,1.0,0.1);
      Eigen::RowVector3d C2(0.1,0.1,1.0);
      Eigen::RowVector3d C3(0.1,1.0,1.0);
      Eigen::RowVector3d C4(1.0,0.1,1.0);
      Eigen::RowVector3d C5(1.0,1.0,0.1);
      */

      Eigen::RowVector3d C0(1.0,0.1,0.1);
      Eigen::RowVector3d C1(0.1,1.0,0.1);
      Eigen::RowVector3d C2(0.1,0.1,1.0);
      Eigen::RowVector3d C3(0.2,0.95,0.95);
      Eigen::RowVector3d C4(0.95,0.2,0.95);
      Eigen::RowVector3d C5(0.95,0.95,0.2);

      //igl::read_triangle_mesh(mesh_file,V,F);//32770 cylinder
      Mesh::readOBJ("/home/n/ADMM/point_cloud/NLC_ADMM-UAV/build/multi_uav.obj", V);
      
      Eigen::MatrixXd C=V;
      double x_up=V.col(0).maxCoeff();
      double x_down=V.col(0).minCoeff();

      double y_up=V.col(1).maxCoeff();
      double y_down=V.col(1).minCoeff();

      double z_up=V.col(2).maxCoeff();
      double z_down=V.col(2).minCoeff();

      std::cout<<viewer.core().camera_eye<<"\n";
      viewer.core().camera_eye = Eigen::Vector3f(0,-10,15);
      //viewer.core().align_camera_center(V);
      
      for(int i=0;i<V.rows();i++)
      {
       
        double x=V(i,0);
        double y=V(i,1);
        double z=V(i,2);

        Eigen::RowVector3d x_tmp=((x-x_down)*C0+(x_up-x)*C3)/(x_up-x_down);
        Eigen::RowVector3d y_tmp=((y-y_down)*C1+(y_up-y)*C4)/(y_up-y_down);
        Eigen::RowVector3d z_tmp=((z-z_down)*C2+(z_up-z)*C5)/(z_up-z_down);

        C.row(i)=(x_tmp+y_tmp+z_tmp)/3.0;
      }
      //if(V.rows()>1)
      //viewer.data().set_points(V,C);
    
    Eigen::RowVector3d CC;
    igl::read_triangle_mesh("/home/n/untitled.obj",BV,BF);
    viewer.data().set_mesh(BV,BF);
    /*
    for(int i=0;i<uav_num;i++)
    {
      if(i<4)
          CC=(4.0-i)/4.0*C0+i/4.0*C5;
      else if(i>8)
          CC=(4.0-(i-8.0))/4.0*C3+(i-8.0)/4.0*C2;
      else
          CC=(4.0-(i-4.0))/4.0*C5+(i-4.0)/4.0*C3;
      //if(point_lists[i].size()>0)
      viewer.data().add_points(point_lists[i][0].transpose(),CC);
      for(int j=0;j<point_lists[i].size()-1;j++)
      {
        viewer.data().add_edges(point_lists[i][j].transpose(), point_lists[i][j+1].transpose(), CC);
      }
    }
    */
    
    
   for(int i=0;i<uav_num;i++)
    {
      if(i==0)
        CC=C0;
      else if(i==1)
        CC=C5;
      else if(i==2)
        CC=C3;
      else if(i==3)
        CC=C2;
      //if(point_lists[i].size()>0)
      viewer.data().add_points(point_lists[i][0].transpose(),CC);
      for(int j=0;j<point_lists[i].size()-1;j++)
      {
        viewer.data().add_edges(point_lists[i][j].transpose(), point_lists[i][j+1].transpose(), CC);
      }
    }

     Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R,G,B,A;
       igl::png::readPNG("/home/n/skin2.png",R,G,B,A);

     viewer.data().set_texture(R,G,B,A);
     viewer.data().use_matcap = true;

    /*
    for(int i=0;i<3;i++)
    {
      for(int j=0;j<=3;j++)
      {
        double x0=-3;
        double x1=-1;
        double y_=-3;
        double z=-2;
        Eigen::RowVector3d p0,p1;
        x0+=i*2;
        x1+=i*2;
        y_+=j*2;
        p0(0)=x0; p0(1)=y_; p0(2)=z;
        p1(0)=x1; p1(1)=y_; p1(2)=z;
        viewer.data().add_edges(p0, p1, 
                                Eigen::RowVector3d(0.5,0.5,0.5));
        
        double y0=-3;
        double y1=-1;
        double x_=-3;

        y0+=i*2;
        y1+=i*2;
        x_+=j*2;
        p0(0)=x_; p0(1)=y0; p0(2)=z;
        p1(0)=x_; p1(1)=y1; p1(2)=z;
        viewer.data().add_edges(p0, p1, 
                                Eigen::RowVector3d(0.5,0.5,0.5));
      }
    }
    */
    
    viewer.launch();
   



  return 0;
}

