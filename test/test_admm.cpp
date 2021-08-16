
#include <igl/opengl/glfw/Viewer.h>

#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>

#include "HighOrderCCD/Utils/CCDUtils.h"
#include "HighOrderCCD/OMPL/OMPL.h"

#include "HighOrderCCD/Optimization/Optimization3D_admm.h"
#include "HighOrderCCD/Optimization/Optimization3D_am.h"

#include "lib/nlohmann/json.hpp"    // https://github.com/nlohmann/json/tree/develop/single_include/nlohmann/json.hpp

USE_PRJ_NAMESPACE

using json = nlohmann::json;

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
              std::ofstream file(meshfile+"_ccd_test.txt");
              if (file.is_open())
              {
                file << "coeff:\n" << coeff << '\n';
                file << "time:" << '\n' <<  my_time << '\n';
              }
              file.close();

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

              igl::write_triangle_mesh("ccd_traj_"+meshfile+".obj",v_ccd,f_ccd);
}

int main(int argc, char *argv[])
{
  if (argc < 2)
	{
		std::cerr << "Syntax: " << argv[0] << " <mesh file>" << std::endl;
		return -1;
  }

  #if 0
  
  #else

  typedef Eigen::MatrixXd Data;

  
  igl::opengl::glfw::Viewer viewer;
 
  viewer.core().background_color<< 1.0f, 1.0f, 1.0f, 1.0f;
  viewer.core().is_animating = true;
  viewer.core().camera_zoom = 2.0f;
  viewer.data().line_width = 1.0f;
  viewer.data().point_size = 3.0f;

  std::ifstream fin("Config File/3D.json");   
  json j = json::parse(fin);
	fin.close();


  lambda=j["lambda"].get<double>();

  epsilon=j["epsilon"].get<double>();
  margin=j["margin"].get<double>();//d hat
  automove = j["auto"].get<int>();
  step_choose = j["step_choose"].get<int>();

  int init_spline=j["init"].get<int>();
  
  offset = j["offset"].get<double>();

  res = j["res"].get<int>();

  int if_exit=j["exit"].get<int>();

  int if_init_ob=j["init_ob"].get<int>();

  double scale=j["scale"].get<double>();

  double stop=j["stop"].get<double>();

  mu=j["mu"].get<double>();
  //lambda=1.0/margin2;

  double whole_time=0;
  
  adaptive_change=false;
  optimize_time=false;
  
  
  
  int dim = kdop_axis.size();
  
  kdop_matrix.resize(3, dim);
  for(int k=0;k<dim;k++)
  {
    kdop_axis[k].normalize();
    kdop_matrix.col(k) = kdop_axis[k];
  }

  aabb_matrix.resize(3, 3);
  for(int k=0;k<3;k++)
  {
    aabb_matrix.col(k) = aabb_axis[k];
  }
  //lambda/=epsilon;

  Eigen::MatrixXd V, BV;
  Eigen::MatrixXi F, BF;
  
  const std::string mesh_file = argv[1];
  igl::read_triangle_mesh(mesh_file,V,F);//32770 cylinder

  std::cout<<"before bvh init\n";
  BVH bvh;
  clock_t time1 = clock();
  
  

  std::vector<Eigen::Matrix3d> face_list;
  face_list.resize(F.rows());
  for(int i=0;i<(int)face_list.size();i++)
  {
    int f0=F(i,0); int f1=F(i,1); int f2=F(i,2);
              
    Eigen::Matrix3d _position;
    _position<<V.row(f0),V.row(f1),V.row(f2);
    face_list[i]=_position;
  }

  clock_t time2 = clock();

  std::cout<<"time_bvh:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
  
  vel_limit=j["vel_limit"].get<double>();;//2
  acc_limit=j["acc_limit"].get<double>();;//2
    
  result_file.open ("result/" +mesh_file + "_result_file_admm.txt");
 
  double init_time=-INFINITY;
  Data spline;



   




  Eigen::MatrixXd  P0,P1;
  P0.resize(6,3);
  P0<<-0.4504048407, -5.624279037e-05 ,-0.0008498807838,
   -0.4504048407 ,-5.624279037e-05 ,-0.0008498807838,
   -0.4504048407 ,-5.624279037e-05 ,-0.0008498807838,
   -0.4504048407 ,-5.624279037e-05, -0.0008498807838,
   -0.4504048407, -5.624279037e-05, -0.0008498807838,
   -0.4504048407, -5.624279037e-05, -0.0008498807838;
    P1.resize(6,3);
    P1<<-0.4484194089,   0.186800794 , 0.2546158723,
-0.4484194089  , 0.186800794 , 0.2546158723,
-0.4484194089  , 0.186800794 , 0.2546158723,
-0.4484194089 ,  0.186800794 , 0.2546158723,
-0.4484194089 ,  0.186800794 , 0.2546158723,
-0.4484194089 ,  0.186800794 , 0.2546158723;

  Eigen::Vector3d c;
              double d;
              c<<-0.02741941545,  -0.6010291643,  -0.7990760005;
              d= 0.1375445874;
  Eigen::MatrixXd  D0,D1;
  D0=P0; D0.setZero();
  D1=P1; D1.setZero();
  bool is_collided= CCD::SelfGJKCCD(P0,D0,P1,D1, offset,0,0, 0, 0); 
  std::cout<<is_collided<<"\n"<<std::flush;
  std::cout<<Optimal_plane::self_barrier_energy(P0,  P1, 
                                  c,  d)<<"\n"<<std::flush;
  Optimal_plane::self_optimal_cd(P0,  P1, 
                                  c,  d);
  std::cout<<c.transpose()<< " "<<d<<"\n\n"<<std::flush;
  #endif
  return 0;
}

