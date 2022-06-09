
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

int main(int argc, char *argv[])
{
  if (argc < 2)
	{
		std::cerr << "Syntax: " << argv[0] << " <mesh file>" << std::endl;
		return -1;
  }

  const std::string mesh_file = argv[1];

  std::ifstream fin("Config File/3D.json");   
  json j = json::parse(fin);
	fin.close();

  lambda=j["lambda"].get<double>();

  epsilon=j["epsilon"].get<double>();
  margin=j["margin"].get<double>();//d hat
  automove = j["auto"].get<int>();

  is_optimal_plane=j["optimal_plane"].get<int>();

  int init_spline=j["init"].get<int>();
  
  offset = j["offset"].get<double>();

  res = j["res"].get<int>();

  int if_exit=j["exit"].get<int>();

  int if_init_ob=j["init_ob"].get<int>();

  double stop=j["stop"].get<double>();

  mu=j["mu"].get<double>();
  
  int gui=j["gui"].get<int>();
    
  vel_limit=j["vel_limit"].get<double>();
  acc_limit=j["acc_limit"].get<double>();
  //lambda=1.0/margin2;
  
  result_file.open ("result/" +mesh_file + "_result_file_admm.txt");

  
  
  
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



  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  
  //igl::read_triangle_mesh(mesh_file,V,F);//32770 cylinder
  Mesh::readOBJ(mesh_file, V);

  std::cout<<"before bvh init\n";
  BVH bvh;
  clock_t time1 = clock();
  
  if(if_init_ob)
  {
    bvh.InitPointcloud(V);
  }

  std::vector<Eigen::RowVector3d> vertex_list;
  vertex_list.resize(V.rows());
  for(int i=0;i<(int)vertex_list.size();i++)
  {              
    vertex_list[i]=V.row(i);
  }

  clock_t time2 = clock();

  std::cout<<"time_bvh:"<<(time2-time1)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
     
 
  std::vector<Eigen::Vector3d> way_points;

  uav_num=1;
  
  if(init_spline==1)
  {
    way_point_init( mesh_file, way_points);
  }
  else if(init_spline==2)
  {
    ompl_init( V, bvh,  way_points);
  }
  piece_num=way_points.size()-1;
  time_weight.resize(piece_num); 
  whole_weight=0;
  for(int i=0;i<piece_num;i++)
  {
    time_weight[i]=1;
    //time_weight[i]=(way_points[i+1]-way_points[i]).norm()/vel_limit;
    
    whole_weight+=time_weight[i];
  }

  combination = Combination<40>::value();

  int num=1000000, turns=0;

  gnorm=1;
  
  ks=1e-8;//1e-3
  kt=1;

  

  double piece_time=20;//init_time;//20

  Conversion<order_num>::convert_matrix();
  
  std::cout<<convert_list[0]<<std::endl;
  //std::cout<<M_head<<std::endl;
  //std::cout<<M_tail<<std::endl;
  //M_convert/=double(Factorial<order_num>::value());
  Data spline, p_slack, p_lambda; 
  Eigen::VectorXd t_slack, t_lambda;
  init_variable(way_points, vertex_list, 
                spline,p_slack,p_lambda, 
                piece_time,t_slack,t_lambda);
  
  double whole_time=0;

  energy_file.open (mesh_file + "_energy_file_admm.txt");

  double energy=  lambda*Energy_admm::bound_energy(spline,piece_time);//lambda*Energy_admm::plane_barrier_energy(spline,c_lists, d_lists) +
    
    for(int sp_id=0;sp_id<piece_num;sp_id++)
    {
        int init=sp_id*(order_num-2);

        Data c_spline = convert_list[sp_id]*spline.block<order_num+1,3>(init,0);

        energy+=Energy_admm::dynamic_energy(c_spline, piece_time);
    }
   double e = ks*Energy::dynamic_energy(spline,piece_time)+ kt*whole_weight*piece_time;
   std::cout<<energy<<" "<<e<<"\n";
  energy_file<<energy;
  energy_file<<" "<<whole_time<<"\n";
  
  if(!gui)
  {
    automove=true;
    while(iter<num) 
    {
      
      if(iter>1 && gnorm<stop)
      {
        
            result_file<<"iter: "<<iter<<std::endl;
            result_file<<"running time: "<<whole_time<<std::endl;
            //result_file<<gnorm<<std::endl;
            result_file<<"point cloud size: "<<V.rows()<<std::endl;
            result_file<<"spline\n"<<spline<<std::endl;
            //result_file<<piece_time<<std::endl;

            log_data(mesh_file, spline, piece_time);

            if(if_exit)
              exit(0);
            else
              automove=false;
            break;
            //log_data(mesh_file, spline, piece_time);
      }
      
      if(iter<turns||automove)
      {

        clock_t time0 = clock();
        std::cout<<"iter: "<<iter<<std::endl;
            
        
        Optimization3D_admm::optimization(spline, piece_time, 
                                          p_slack, t_slack, 
                                          p_lambda, t_lambda,
                                          vertex_list, bvh);
    
        clock_t time1 = clock();
        whole_time+=(time1-time0)/(CLOCKS_PER_SEC/1000);

        energy_file<<" "<<whole_time<<"\n";

        std::cout<<"time:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
        //std::cout<<p0.size()<<std::endl;
        iter++;
      }
    }

  }
  else
  {
    
      igl::opengl::glfw::Viewer viewer;
 
      viewer.core().background_color<< 1.0f, 1.0f, 1.0f, 1.0f;
      viewer.core().is_animating = true;
      viewer.core().camera_zoom = 2.0f;
      viewer.data().line_width = 1.0f;
      viewer.data().point_size = 3.0f;

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
      
      viewer.data().line_width = 5.0f;
      viewer.data().point_size = 1.5f;
      
      Eigen::MatrixXd C=V;
      
      Eigen::RowVector3d C0(1.0,0.1,0.1);
      Eigen::RowVector3d C1(0.1,1.0,0.1);
      Eigen::RowVector3d C2(0.1,0.1,1.0);
      Eigen::RowVector3d C3(0.1,1.0,1.0);
      Eigen::RowVector3d C4(1.0,0.1,1.0);
      Eigen::RowVector3d C5(1.0,1.0,0.1);
      
     

      double x_up=V.col(0).maxCoeff();
      double x_down=V.col(0).minCoeff();

      double y_up=V.col(1).maxCoeff();
      double y_down=V.col(1).minCoeff();

      double z_up=V.col(2).maxCoeff();
      double z_down=V.col(2).minCoeff();
      
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
      viewer.data().set_points(V,C);

      double tree_size=(subdivide_tree.size()-1)/5.0;
      for(unsigned int k=0;k<subdivide_tree.size();k++)
      {
        int sp_id=std::get<0>(subdivide_tree[k]);
        Eigen::MatrixXd basis=std::get<2>(subdivide_tree[k]);
        Eigen::MatrixXd bz;
        Eigen::RowVector3d CC;
        if(k<2*tree_size)
          CC=(2*tree_size-k)/(2*tree_size)*C0+k/(2*tree_size)*C5;
        else if(k>3*tree_size)
          CC=(2*tree_size-(k-3*tree_size))/(2*tree_size)*C3+(k-3*tree_size)/(2*tree_size)*C2;
        else
          CC=(tree_size-(k-2*tree_size))/tree_size*C5+(k-2*tree_size)/tree_size*C3;
        
        bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);
      
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
            viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], CC);//Eigen::RowVector3d(0.2,0.8,0.8));
          }
        }
        else
        {
          for(int j=0;j<=order_num;j++)
          {
            viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], CC);//Eigen::RowVector3d(0.8,0.2,0.8));
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
              //int edge_iter=0;
              for(unsigned int k=0;k<subdivide_tree.size();k++)
              {
                int sp_id=std::get<0>(subdivide_tree[k]);
                Eigen::MatrixXd basis=std::get<2>(subdivide_tree[k]);
                Eigen::MatrixXd bz;
                bz=spline.block<order_num+1,3>(sp_id*(order_num-2),0);

                Eigen::RowVector3d CC;
                if(k<2*tree_size)
                  CC=(2*tree_size-k)/(2*tree_size)*C0+k/(2*tree_size)*C5;
                else if(k>3*tree_size)
                  CC=(2*tree_size-(k-3*tree_size))/(2*tree_size)*C3+(k-3*tree_size)/(2*tree_size)*C2;
                else
                  CC=(tree_size-(k-2*tree_size))/tree_size*C5+(k-2*tree_size)/tree_size*C3;
                
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
                    viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], CC);//Eigen::RowVector3d(0.2,0.8,0.8));
                  }
                }
                else
                {
                  for(int j=0;j<=order_num;j++)
                  {
                    viewer.data().add_edges(P[j], P[(j+1)%(order_num+1)], CC);//Eigen::RowVector3d(0.8,0.2,0.8));
                  }
                }
              }
            }
            
            if(iter>1 && gnorm<stop && is_write)
            {
              
                  is_write = false;
                  result_file<<iter<<std::endl;
                  result_file<<whole_time<<std::endl;
                  result_file<<gnorm<<std::endl;
                  result_file<<V.rows()<<" "<<F.rows()<<std::endl;
                  result_file<<spline<<std::endl;
                  result_file<<piece_time<<std::endl;

                  log_data(mesh_file, spline, piece_time);


                  if(if_exit)
                    exit(0);
                  else
                    automove=false;
                  
                  //log_data(mesh_file, spline, piece_time);
            }
            
            if(iter<turns||automove)
            {

              clock_t time0 = clock();
              std::cout<<"iter: "<<iter<<std::endl;
                  
              
              Optimization3D_admm::optimization(spline, piece_time, 
                                                p_slack, t_slack, 
                                                p_lambda, t_lambda,
                                                vertex_list, bvh);
        
              clock_t time1 = clock();
              whole_time+=(time1-time0)/(CLOCKS_PER_SEC/1000);
              std::cout<<"time:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl<<std::endl;
              //std::cout<<p0.size()<<std::endl;

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
              /*
              std::ofstream curve_file;
              curve_file.open (mesh_file + "_curve_file"+std::to_string(iter)+".txt");
              for(int i=0;i<ccd_traj.size();i++)
              {
                curve_file<<ccd_traj[i].transpose()<<"\n";
              }
              */
              iter++;
            }
        }  
        
          return false;
      };
    
      viewer.callback_pre_draw = pre_draw;
      viewer.callback_key_down = key_down;
    
      viewer.launch();
     
  }
  
  
  return 0;
}

