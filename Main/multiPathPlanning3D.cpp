#include <igl/opengl/glfw/Viewer.h>

#include <igl/read_triangle_mesh.h>

#include "HighOrderCCD/Utils/CCDUtils.h"

#include "HighOrderCCD/Optimization/Optimization3D_multi.h"

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

  //int init_spline=j["init"].get<int>();
  
  offset = j["offset"].get<double>();

  res = j["res"].get<int>();

  int if_exit=j["exit"].get<int>();

  int if_init_ob=j["init_ob"].get<int>();  

  double stop=j["stop"].get<double>();

  mu=j["mu"].get<double>();
  //lambda=1.0/margin2;

  double whole_time=0;
  
  adaptive_change=false;
  optimize_time=false;

  int dim = kdop_axis.size();
  for(int k=0;k<dim;k++)
  {
    kdop_axis[k].normalize();
  }
  //lambda/=epsilon;

  Eigen::MatrixXd V, BV;
  Eigen::MatrixXi F, BF;
  
  const std::string mesh_file = argv[1];
  igl::read_triangle_mesh(mesh_file,V,F);//32770 cylinder
  
  vel_limit=2.0;
  acc_limit=2.0;
    
  result_file.open ("result/" +mesh_file + "_result_file_multi.txt");
  
  
  uav_num=4;//2
  piece_num=3;
  trajectory_num = (order_num+1)+(piece_num-1)*(order_num+1-3);

  time_weight.resize(piece_num);
  whole_weight=0;
  for(int k=0;k<piece_num;k++)
  {
    time_weight[k]=1;
    //time_weight[i]=(way_points[i+1]-way_points[i]).norm()/vel_limit;
    whole_weight+=time_weight[k];
  }
  std::vector<std::vector<Eigen::RowVector3d>> way_points_list;
  way_points_list.resize(uav_num);
  
  way_points_list[0].push_back(Eigen::Vector3d(2, 2 ,1));
  way_points_list[0].push_back(Eigen::Vector3d(1, 0, -0.2));
  way_points_list[0].push_back(Eigen::Vector3d(-2, 0, -0.2));
  way_points_list[0].push_back(Eigen::Vector3d(-2, -2, -1));

  way_points_list[1].push_back(Eigen::Vector3d(2.5, 2 ,-1));
  way_points_list[1].push_back(Eigen::Vector3d(2, 0, 0.2));
  way_points_list[1].push_back(Eigen::Vector3d(-2, 0, 0.2));
  way_points_list[1].push_back(Eigen::Vector3d(-2.5, -2, 1));

  way_points_list[2].push_back(Eigen::Vector3d(-2, 2 ,1));
  way_points_list[2].push_back(Eigen::Vector3d(-1, 0, -0.4));
  way_points_list[2].push_back(Eigen::Vector3d(2, 0, -0.4));
  way_points_list[2].push_back(Eigen::Vector3d(2, -2, -1));

  way_points_list[3].push_back(Eigen::Vector3d(-2.5, 2 ,-1));
  way_points_list[3].push_back(Eigen::Vector3d(-2, 0, 0.4));
  way_points_list[3].push_back(Eigen::Vector3d(2, 0, 0.4));
  way_points_list[3].push_back(Eigen::Vector3d(2.5, -2, 1));

  max_step=1.0;

  combination = Combination<40>::value();

  //Eigen::MatrixXd V,O;
  
  //std::vector<Eigen::MatrixXd> p0,p1,d0;
  //V*=0.000000001;
  int num=1000000, turns=0;

  gnorm=1;
  
  ks=1e-3;
  kt=1;

  Conversion<order_num>::convert_matrix();
  
  std::cout<<convert_list[0]<<std::endl;

  std::vector<Data> spline_list; std::vector<double> piece_time_list;
  std::vector<Data> p_slack_list; std::vector<Eigen::VectorXd> t_slack_list;
  std::vector<Data> p_lambda_list; std::vector<Eigen::VectorXd> t_lambda_list;
  
  spline_list.resize(uav_num); piece_time_list.resize(uav_num);
  p_slack_list.resize(uav_num); t_slack_list.resize(uav_num);
  p_lambda_list.resize(uav_num); t_lambda_list.resize(uav_num);

  for(int i=0;i<uav_num;i++)
  {
      double piece_time=20;

      Data spline;
      spline.resize(trajectory_num,3);
      
      std::vector<Eigen::RowVector3d> way_points=way_points_list[i];

      spline.row(0)=way_points[0];
      for(int k=0;k<piece_num;k++)
      {
        for(int j=0;j<=order_num-2;j++)
        {
          spline.row(j+k*(order_num-2)+1)=double(order_num-2-j)/(order_num-2)*way_points[k]+(double)j/(order_num-2)*way_points[k+1];
        }
      }
      spline.row(trajectory_num-1)=way_points[piece_num];

      spline.row(1)=spline.row(0);
      spline.row(trajectory_num-2)=spline.row(trajectory_num-1);
      
      Data p_slack, p_lambda;
      p_lambda.resize((order_num+1)*piece_num,3); p_lambda.setZero();
      p_slack.resize((order_num+1)*piece_num,3);
      for(int sp_id=0; sp_id<piece_num; sp_id++)
      {
        p_slack.block<order_num+1,3>(sp_id*(order_num+1),0)=convert_list[sp_id]*spline.block<order_num+1,3>(sp_id*(order_num-2),0);
      }

      Eigen::VectorXd t_slack, t_lambda;
      t_lambda.resize(piece_num); t_lambda.setZero();
      t_slack.resize(piece_num);
      for(int sp_id=0; sp_id<piece_num; sp_id++)
      {
        t_slack(sp_id)=piece_time;
      }

      

      spline_list[i]=spline; piece_time_list[i]=piece_time;
      p_slack_list[i]=p_slack; t_slack_list[i]=t_slack;
      p_lambda_list[i]=p_lambda; t_lambda_list[i]=t_lambda;
  }


  M_dynamic=Dynamic<order_num, der_num>::dynamic_matrix();//Dynamic<order_num,1>::dynamic_matrix()+
  
  subdivide_tree.resize(piece_num*res);
  Eigen::MatrixXd basis;
  
  for(int k=0;k<res;k++)
  {
    double a=k/double(res),b=(k+1)/double(res);

    Blossom<order_num>::coefficient(basis, a, b);

    for(int i=0;i<piece_num;i++)
    { 
      std::pair<double,double> range(a,b);
      subdivide_tree[i*res+k]=std::make_tuple(i,range,basis*convert_list[i]);
    }
  }

  std::cout<<"before bvh init\n";
  BVH bvh;
  clock_t time1 = clock();
  if(if_init_ob)
  {
    bvh.InitObstacle(V,F);
  }
    

  clock_t time2 = clock();

  std::cout<<"time_obstacle:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
  //std::cout<<F_<<std::endl;
  
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IR(640,800);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IG(640,800);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IB(640,800);
  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IA(640,800);

  //bool draw_init=true;

  bool is_write=true;

  const auto &key_down = [&](igl::opengl::glfw::Viewer &,unsigned char key,int mod)->bool
  {
    switch(key)
    {
      case 'R':
      {
        iter=0;
        turns=0;
        break;
      }
      case 'C':
      {
        break;
      }
      case ' ':
      {
        //viewer.core().draw_buffer(viewer.data(),false,IR,IG,IB,IA);
        //igl::png::writePNG(IR,IG,IB,IA, "../png/" + std::to_string(turns) + ".png");
        turns=iter+1;
        break;
      }
      case '.':
      {
        automove=false;
        break;
      }
      case ',':
      {
        automove=true;
        break;
      }
      default:
        return false;
    }

    return true;
  };
  viewer.data().line_width = 1.0f;
        
  viewer.data().set_mesh(V,F);

  for(unsigned int k=0;k<subdivide_tree.size();k++)
  {
    int sp_id=std::get<0>(subdivide_tree[k]);
    Eigen::MatrixXd basis=std::get<2>(subdivide_tree[k]);
    Eigen::MatrixXd bz;
    for(int i=0;i<uav_num;i++)
    {
      bz=spline_list[i].block<order_num+1,3>(sp_id*(order_num-2),0);
    
      std::vector<Eigen::RowVector3d> P(order_num+1);
      for(int j=0;j<=order_num;j++)
      {
        P[j].setZero();
        for(int j0=0;j0<=order_num;j0++)
        {
          P[j]+=basis(j,j0)*bz.row(j0);
        }
      }
      if(k%2==1)
      {
        for(int j=0;j<=order_num;j++)
        {
          viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], Eigen::RowVector3d(0.2,0.8,0.8));
        }
      }
      else
      {
        for(int j=0;j<=order_num;j++)
        {
          viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], Eigen::RowVector3d(0.8,0.2,0.8));
        }
      }

    }
    
  }
   
  const auto &pre_draw = [&](igl::opengl::glfw::Viewer & )->bool
  {  
     if(iter<num) 
     {
       if(iter<turns||automove)
       {
          viewer.data().clear_edges();
                    
          //viewer.data().line_width = 5.0f;
          for(unsigned int k=0;k<subdivide_tree.size();k++)
          {
            int sp_id=std::get<0>(subdivide_tree[k]);
            Eigen::MatrixXd basis=std::get<2>(subdivide_tree[k]);
            Eigen::MatrixXd bz;
            for(int i=0;i<uav_num;i++)
            {
              bz=spline_list[i].block<order_num+1,3>(sp_id*(order_num-2),0);
            
              std::vector<Eigen::RowVector3d> P(order_num+1);
              for(int j=0;j<=order_num;j++)
              {
                P[j].setZero();
                for(int j0=0;j0<=order_num;j0++)
                {
                  P[j]+=basis(j,j0)*bz.row(j0);
                }
              }
              if(k%2==1)
              {
                for(int j=0;j<=order_num;j++)
                {
                  viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], Eigen::RowVector3d(0.2,0.8,0.8));
                }
              }
              else
              {
                for(int j=0;j<=order_num;j++)
                {
                  viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], Eigen::RowVector3d(0.8,0.2,0.8));
                }
              }

            }
            
          }
        }
        
        if(gnorm<stop && is_write)
        {
          if(optimize_time)
          {
            if(adaptive_change)
            {
              is_write = false;
              result_file<<iter<<std::endl;
              result_file<<whole_time<<std::endl;
              result_file<<gnorm<<std::endl;
              result_file<<V.rows()<<" "<<F.rows()<<std::endl;
              //result_file<<spline<<std::endl;
              //result_file<<piece_time<<std::endl;

              if(if_exit)
                exit(0);
              else
                automove=false;

            }

            adaptive_change=true;

          }
          else
          {
            result_file<<iter<<std::endl;
            result_file<<whole_time<<std::endl;
            result_file<<gnorm<<std::endl;
            result_file<<V.rows()<<" "<<F.rows()<<std::endl;
            //result_file<<spline<<std::endl;

          
            optimize_time=true;
          }
          
          
        }
        
        if(iter<turns||automove)
        {

          clock_t time0 = clock();
          std::cout<<"iter: "<<iter<<std::endl;
              
         
            
          Optimization3D_multi::optimization(spline_list, piece_time_list, 
                                              p_slack_list, t_slack_list, 
                                              p_lambda_list, t_lambda_list,
                                              V, F, bvh);
            
          
                    
          clock_t time1 = clock();
          whole_time+=(time1-time0)/(CLOCKS_PER_SEC/1000);
          std::cout<<"time:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
          //std::cout<<p0.size()<<std::endl;
          iter++;
        }
     }  
     
      return false;
  };
  
  
  viewer.callback_pre_draw = pre_draw;
  viewer.callback_key_down = key_down;
  
  
  
  
  viewer.launch();

  #endif
  return 0;
}

