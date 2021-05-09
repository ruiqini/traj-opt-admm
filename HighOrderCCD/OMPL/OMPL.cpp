#include "OMPL.h"
#include <ompl/geometric/planners/rrt/RRTConnect.h>
#include <ompl/geometric/planners/rrt/RRTstar.h>
#include <ompl/geometric/SimpleSetup.h>
#include <ompl/base/MotionValidator.h>


USE_PRJ_NAMESPACE

namespace og=ompl::geometric;


class myMotionValidator : public ob::MotionValidator
{
public:
   
    myMotionValidator(ob::SpaceInformation* si, 
                      Eigen::MatrixXd _V,
                      Eigen::MatrixXi _F,
                       BVH& _bvh) : ob::MotionValidator(si)
    {
      V=_V;
      F=_F;
      bvh=&_bvh;
    }
 
    myMotionValidator(const ob::SpaceInformationPtr &si, 
                      Eigen::MatrixXd _V,
                      Eigen::MatrixXi _F,
                       BVH& _bvh) : ob::MotionValidator(si)
    {
      V=_V;
      F=_F;
      bvh=&_bvh;
    }

    ~myMotionValidator() override = default;
    bool checkMotion (const ob::State *s1, const ob::State *s2) const override
    {
      bool valid=true;

      const ob::RealVectorStateSpace::StateType* sv=dynamic_cast<const ob::RealVectorStateSpace::StateType*>(s1);
      const ob::RealVectorStateSpace::StateType* gv=dynamic_cast<const ob::RealVectorStateSpace::StateType*>(s2);

      double* sDouble=sv->values;
      double* gDouble=gv->values;
      //Cold ret=Eigen::Map<const Eigen::Matrix<double,-1,1>>(sv->values,_body.nrDOF()).cast<scalarD>();
      //int num=3;
      Eigen::MatrixXd edge(2,3);
      edge(0,0)=sDouble[0];
      edge(0,1)=sDouble[1];
      edge(0,2)=sDouble[2];

      edge(1,0)=gDouble[0];
      edge(1,1)=gDouble[1];
      edge(1,2)=gDouble[2];
      std::cout<<edge<<"\n\n";
      std::vector<unsigned int> collision_pair;
        //bvh.CheckCollision(collision_pair,margin);
      bvh->EdgeCollision(edge, collision_pair,offset);

      int collision_size=collision_pair.size();
      std::cout<<collision_size<<"\n\n";
      if(collision_size==0)
         return valid;
      else
      {
        for(int i=0;i<collision_size;i++)
        {
            //int ob_id=*it;
            int ob_id=collision_pair[i];

            int f0=F(ob_id,0); int f1=F(ob_id,1); int f2=F(ob_id,2);
            
            Eigen::Matrix3d _position;
            _position<<V.row(f0),V.row(f1),V.row(f2);

            bool is_collided= CCD::GJKDCD(edge,_position, offset);  //cgal
            if(is_collided)
            {
                return !valid;
            }     
            
        }

      }

      return valid;
      //std::cout<<"check3"<<std::endl;
      //INFOV("%s",valid?"true":"false")
      
       
    }
    //This function assumes s1 is valid. Check if the path between two states is valid. Also compute the last state that was valid and the time of that state. The time is used to parametrize the motion from s1 to s2, s1 being at t = 0 and s2 being at t = 1
    bool checkMotion (const ob::State *s1, const ob::State *s2, std::pair< ob::State *, double > &lastValid) const override
    {
      bool valid=true;
      
      return valid;
    }
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    BVH* bvh;
    
};
//OMPL
OMPL::OMPL(Eigen::VectorXd lowerBound, Eigen::VectorXd upperBound, Eigen::MatrixXd _V,
           Eigen::MatrixXi _F,
           BVH& _bvh)
{
  int dim=lowerBound.size();
  _state=ob::StateSpacePtr(new ob::RealVectorStateSpace(dim));
  ompl::base::RealVectorBounds bounds(dim);
  //lower
  for(int i=0; i<dim; i++)
    bounds.setLow(i,lowerBound(i));
  //upper
  for(int i=0; i<dim; i++)
    bounds.setHigh(i,upperBound(i));
  //setup
  _state->as<ob::RealVectorStateSpace>()->setBounds(bounds);
  
  _stateInfo=ob::SpaceInformationPtr(new ob::SpaceInformation(_state));

  myMotionValidator *motion_validator=new myMotionValidator(_stateInfo,_V,_F,_bvh);
  ob::MotionValidatorPtr mv(motion_validator);
  
  _stateInfo->setMotionValidator(mv);//std::make_shared<myMotionValidator>(_stateInfo));

  //_stateInfo->setMotionValidator(std::make_shared<myMotionValidator>(_stateInfo));
  _stateInfo->setup();
  /*
  _stateInfo->setStateValidityChecker([&](const ob::State* s) {
    return isValid(s);
  });
  */
}
/*
bool OMPL::isValid(const ob::State* s)
{
  bool valid=true;
  const ob::RealVectorStateSpace::StateType* sv=dynamic_cast<const ob::RealVectorStateSpace::StateType*>(s);
  double* stateDouble=sv->values;
  //Cold ret=Eigen::Map<const Eigen::Matrix<double,-1,1>>(sv->values,_body.nrDOF()).cast<scalarD>();
  
  //INFOV("%s",valid?"true":"false")
  return valid;
}
*/

void OMPL::getPath(std::vector<Eigen::VectorXd> &path) //const const 
{
  path=_path;
}

bool OMPL::planRRT( Eigen::VectorXd start,  Eigen::VectorXd goal, 
                      Eigen::MatrixXd _V,
                      Eigen::MatrixXi _F,
                      BVH& _bvh, int time)//,std::function<bool(const Mesh&)> fn
{
  ob::ProblemDefinitionPtr prob(new ob::ProblemDefinition(_stateInfo));
  //set start
  //std::shared_ptr<ob::Goal> goal(new GoalFn(_stateInfo));
  std::shared_ptr<ob::RealVectorStateSpace::StateType> goalState(new ob::RealVectorStateSpace::StateType);
  std::shared_ptr<ob::RealVectorStateSpace::StateType> startState(new ob::RealVectorStateSpace::StateType);
  
  startState->values=start.data();
  prob->addStartState(startState.get());
  //set goal
    std::cout<<"goal"<<std::endl;

  
  goalState->values=NULL;
  //clock_t stamp=clock();
    
  goalState->values=goal.data();

  prob->setGoalState(goalState.get());
    
    
  myMotionValidator mv(_stateInfo,_V,_F,_bvh);
  std::cout<<"test:"<<mv.checkMotion(startState.get(),goalState.get())<<std::endl;
    /*
    if(goalState->values)
      prob->setGoalState(goalState.get());
    else 
      return false;
    */
  
  std::cout<<"planner"<<std::endl;
  //planner
  std::shared_ptr<ob::Planner> planner;

  //prob->setOptimizationObjective(getClearanceObjective(_stateInfo));

  
  planner.reset(new og::RRTConnect(_stateInfo));

  planner->setProblemDefinition(prob);
  planner->setup();
  _stateInfo->printSettings(std::cout);
  prob->print(std::cout);


  //solve
  ob::PlannerStatus solved=planner->ob::Planner::solve(time);

  if(solved==ob::PlannerStatus::EXACT_SOLUTION) {
    std::cout << "Found solution:" << std::endl;

    og::PathGeometric path( dynamic_cast< const og::PathGeometric& >( *prob->getSolutionPath()));

    const std::vector< ob::State* > &states = path.getStates();
    ob::State *state;
    int dim=start.size();
    for( size_t i = 0 ; i < states.size( ) ; ++i )
    {
        state = states[ i ]->as< ob::State >( );

        Eigen::VectorXd ss(dim);
        ss(0) = state->as<ob::RealVectorStateSpace::StateType>()->values[0];
        ss(1) = state->as<ob::RealVectorStateSpace::StateType>()->values[1];
        ss(2) = state->as<ob::RealVectorStateSpace::StateType>()->values[2];
       
        _path.push_back(ss);
    }
    
    return true;
  } else {
    std::cout << "No solution found, status=" << solved << "!" <<std::endl;
    return false;
  }
}
/*
void OMPL::returnToStart(const Cold& start)
{
  sizeType id=(sizeType)_path.size()-2;
  ASSERT_MSG(id>=0,"Cannot return to start!")
  while((_path[id]-start).cwiseAbs().maxCoeff()>1E-6f) {
    ASSERT_MSG(id>=0,"Cannot return to start!")
    _path.push_back(_path[id--]);
  }
  _path.push_back(_path[id]);
}
*/

