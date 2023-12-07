#include "mfem.hpp"
#include "Stimulus.hh"
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

double StimulusCollection::Eval(ElementTransformation& T, const IntegrationPoint &ip)
{
   double x[3];
   Vector transip(x, 3);
 
   T.Transform(ip, transip);

   double result = 0;
   for (std::size_t istim=0; istim<stim_.size(); istim++)
   {
      result += dt_*stim_[istim].eval(time_, T.ElementNo, transip);
   }
   return result;
}

void StimulusCollection::add(Stimulus newStim)
{
   stim_.push_back(newStim);
}