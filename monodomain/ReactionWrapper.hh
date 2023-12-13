
#ifndef REACTIONWRAPPER_HH
#define REACTIONWRAPPER_HH

#include "mfem.hpp"
#include "ReactionManager.hh"
#include "ThreadServer.hh"
#include <memory>
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

class ReactionWrapper
{
 public:
   ReactionWrapper(const double new_dt, const std::vector<std::string> new_objectNames, const ThreadTeam& group, const std::vector<int>& cellTypes);
   void Initialize();
   void Calc(const double _dt);
   Vector& getVmReadwrite();
   Vector& getIionReadwrite();
   const Vector& getVmReadonly() const;
   const Vector& getIionReadonly() const;
   std::vector<std::string> Filednames;
   std::vector<std::string> FieldUnits;
   std::vector<int> handle;
   std::vector<double> value;

   //void GetCheckpointInfo(std::vector<std::string>& fieldNames,
                 //         std::vector<std::string>& fieldUnits) const;
   void GetValue(int iCell) ;
   vector<int> allCellTypes() ;

   //private:
   
   lazy_array<double> Vm;
   lazy_array<double> iStim;
   lazy_array<double> dVm;
   
   int nCells;
   double dt;
   std::string objectName;
   ThreadTeam threadGroup;

   ReactionManager reaction;

   Vector Vm_vector;
   Vector Iion_vector;
   //Vector Output;
};

#endif