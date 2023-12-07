#include "mfem.hpp"
#include "ReactionWrapper.hh"
#include "reactionFactory.hh"
#include "ThreadServer.hh"
#include "object.h"
#include "object_cc.hh"
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mfem;

ReactionWrapper::ReactionWrapper(const double new_dt, const std::vector<std::string> new_objectNames, const ThreadTeam& group, const std::vector<int>& cellTypes)
: dt(new_dt), threadGroup(group)
{
      
   nCells = cellTypes.size();

   Vm.resize(nCells);
   dVm.resize(nCells);
   iStim.resize(nCells);

   for (auto objectName : new_objectNames) {
      reaction.addReaction(objectName);
   }
   reaction.create(dt, cellTypes, threadGroup);

   reaction.getCheckpointInfo(Filednames,FieldUnits);
   for (unsigned ii=0; ii<Filednames.size(); ++ii)
   {
      handle.push_back(reaction.getVarHandle(Filednames[ii]));
   }
   value.resize(handle.size());

   rw_array_ptr<double> Vm_ptr = Vm.readwrite(CPU);
   Vm_vector.SetDataAndSize(Vm_ptr.raw(), nCells);
   rw_array_ptr<double> Iion_ptr = dVm.readwrite(CPU);
   Iion_vector.SetDataAndSize(Iion_ptr.raw(), nCells);

}

void ReactionWrapper::Initialize()
{
   reaction.initializeMembraneState(Vm);
}

void ReactionWrapper::Calc()
{
   reaction.calc(dt, Vm, iStim, dVm);
}

Vector& ReactionWrapper::getVmReadwrite()
{
   rw_array_ptr<double> Vm_ptr = Vm.useOn(CPU);
   return Vm_vector;
}

Vector& ReactionWrapper::getIionReadwrite()
{
   rw_array_ptr<double> Iion_ptr = dVm.useOn(CPU);
   return Iion_vector;
}

const Vector& ReactionWrapper::getVmReadonly() const
{
   ro_array_ptr<double> Vm_ptr = Vm.useOn(CPU);
   return Vm_vector;
}

const Vector& ReactionWrapper::getIionReadonly() const
{
   ro_array_ptr<double> Iion_ptr = dVm.useOn(CPU);
   return Iion_vector;
}


void ReactionWrapper::GetValue(int iCell) 
{
   for (unsigned ii=0; ii<handle.size(); ++ii){
      value[ii] = reaction.getValue(iCell, handle[ii]);
   };
}

vector<int> ReactionWrapper::allCellTypes() 
{
   return reaction.allCellTypes();
}

