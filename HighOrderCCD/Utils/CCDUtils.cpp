#include "CCDUtils.h"

PRJ_BEGIN

std::vector<std::vector<long int>> combination;
//std::vector<double> element_len;
//int element_num;
int trajectory_num, piece_num, uav_num;
int res;
int iter;
double epsilon , 
       distance, 
       margin, 
       max_step, lambda, wolfe, offset, 
       gnorm, gtime, tnorm,
       mu;

std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd> > subdivide_tree;
std::vector<std::vector< Eigen::MatrixXd>> A_list, A_vel_list, A_acc_list;
bool automove, step_choose, adaptive_change, optimize_time;
std::ofstream  result_file, init_file;
Eigen::MatrixXd M_dynamic;

std::vector<bool> reach_target;

std::vector<double> time_weight;
double whole_weight;
std::vector<Eigen::MatrixXd> convert_list;

std::vector<std::vector<std::vector<bool>>> is_self_seperate;
std::vector<std::vector<std::vector<Eigen::Vector3d>>> self_seperate_c;
std::vector<std::vector<std::vector<double>>> self_seperate_d;

bool is_optimal_plane;

double piece_time, kt, ks;
double vel_limit, acc_limit;

Eigen::MatrixXd aabb_matrix;
Eigen::MatrixXd kdop_matrix;

std::vector<Eigen::Vector3d> aabb_axis={
                                          Eigen::Vector3d( 1, 0, 0),
                                          Eigen::Vector3d( 0, 1, 0),
                                          Eigen::Vector3d( 0, 0, 1),//x y z
                                   };
std::vector<Eigen::Vector3d> axis={
                                          Eigen::Vector3d( 1, 0, 0),
                                          Eigen::Vector3d( 0, 1, 0),
                                          Eigen::Vector3d( 0, 0, 1),//x y z
                                   };
std::vector<Eigen::Vector3d> kdop_axis={
                                          Eigen::Vector3d( 1, 0, 0),
                                          Eigen::Vector3d( 0, 1, 0),
                                          Eigen::Vector3d( 0, 0, 1),//x y z

                                          Eigen::Vector3d( 1, 1, 1),
                                          Eigen::Vector3d( 1,-1, 1),
                                          Eigen::Vector3d( 1, 1,-1),
                                          Eigen::Vector3d( 1,-1,-1),//

                                          Eigen::Vector3d( 0, 1, 1),
                                          Eigen::Vector3d( 0, 1,-1),
                                          Eigen::Vector3d( 1, 0, 1),
                                          Eigen::Vector3d( 1, 0,-1),
                                          Eigen::Vector3d( 1, 1, 0),
                                          Eigen::Vector3d( 1,-1, 0),//
                                          
                                          
                                          Eigen::Vector3d( 0, 2, 1),
                                          Eigen::Vector3d( 0, 2,-1),
                                          Eigen::Vector3d( 0, 1, 2),
                                          Eigen::Vector3d( 0, 1,-2),

                                          Eigen::Vector3d( 2, 0, 1),
                                          Eigen::Vector3d( 2, 0,-1),
                                          Eigen::Vector3d( 1, 0, 2),
                                          Eigen::Vector3d( 1, 0,-2),

                                          Eigen::Vector3d( 2, 1, 0),
                                          Eigen::Vector3d( 2,-1, 0),
                                          Eigen::Vector3d( 1, 2, 0),
                                          Eigen::Vector3d( 1,-2, 0),
                                          
                                          Eigen::Vector3d( 1, 2, 1),
                                          Eigen::Vector3d( 1, 2,-1),
                                          Eigen::Vector3d( 1,-2, 1),
                                          Eigen::Vector3d(-1, 2, 1),

                                          Eigen::Vector3d( 1, 1, 2),
                                          Eigen::Vector3d( 1, 1,-2),
                                          Eigen::Vector3d( 1,-1, 2),
                                          Eigen::Vector3d(-1, 1, 2),

                                          Eigen::Vector3d( 2, 1, 1),
                                          Eigen::Vector3d( 2, 1,-1),
                                          Eigen::Vector3d( 2,-1, 1),
                                          Eigen::Vector3d(-2, 1, 1),

                                          Eigen::Vector3d( 2, 2, 1),
                                          Eigen::Vector3d( 2, 2,-1),
                                          Eigen::Vector3d( 2,-2, 1),
                                          Eigen::Vector3d(-2, 2, 1),

                                          Eigen::Vector3d( 2, 1, 2),
                                          Eigen::Vector3d( 2, 1,-2),
                                          Eigen::Vector3d( 2,-1, 2),
                                          Eigen::Vector3d(-2, 1, 2),

                                          Eigen::Vector3d( 1, 2, 2),
                                          Eigen::Vector3d( 1, 2,-2),
                                          Eigen::Vector3d( 1,-2, 2),
                                          Eigen::Vector3d(-1, 2, 2),

                                          };

PRJ_END