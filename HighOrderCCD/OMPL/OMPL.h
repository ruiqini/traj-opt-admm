#ifndef OMPL_H
#define OMPL_H

#include "HighOrderCCD/Utils/CCDUtils.h"

#include "HighOrderCCD/CCD/CCD.h"
#include "HighOrderCCD/BVH/BVH.h"

#include <vector>
#include <memory>
#include <ctime>
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/spaces/RealVectorBounds.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/objectives/StateCostIntegralObjective.h>

PRJ_BEGIN

namespace ob=ompl::base;
class OMPL 
{
public:

  

  OMPL(Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound,
        Eigen::MatrixXd _V, std::vector<std::vector<Eigen::MatrixXd>> _edges,
        BVH& _bvh);
  bool checkMotion(const ob::State *s1, const ob::State *s2);
  int nrBroad() const;
  int nrNarrow() const;
  void getPath(std::vector<Eigen::Vector3d> &path);//const;const 
  bool planRRT( Eigen::Vector3d start,  Eigen::Vector3d goal,
                      Eigen::MatrixXd _V, std::vector<std::vector<Eigen::MatrixXd>> _edges,
                      BVH& _bvh,
                       int time=1200);//std::function<bool(const Mesh&)> fn
  ob::OptimizationObjectivePtr getClearanceObjective(const ob::SpaceInformationPtr& si);

  //void returnToStart(const Mesh& start);
protected:
  //std::shared_ptr<ob::RealVectorStateSpace> _state;
  ob::StateSpacePtr _state;
  //std::shared_ptr<ob::SpaceInformation> _stateInfo;
  ob::SpaceInformationPtr _stateInfo;
  //options
  bool _AECC,_env,_self;
  //result
  std::vector<Eigen::Vector3d> _path;
  int _nrBroad;
  int _nrNarrow;
};

PRJ_END

#endif
