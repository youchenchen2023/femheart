                                                                                                                                                                        #include "mfem.hpp"
#include "object.h"
#include "object_cc.hh"
#include "ddcMalloc.h"
#include "pio.h"
#include "pioFixedRecordHelper.h"
#include "units.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <cassert>
#include <memory>
#include <set>
#include <ctime>
#include <chrono>
#include <cstdio>
#include <dirent.h>
#include <regex.h>
#include <unistd.h>
#include <sys/stat.h>
#include "MatrixElementPiecewiseCoefficient.hpp"
#include "ReactionWrapper.hh"
#include "Stimulus.hh"


//#define StartTimer(x)
//#define EndTimer()

using namespace mfem;
using namespace std;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

//Stolen from SingleCell
class Timeline
{
 public:
   Timeline(double dt, double duration)
   {
      maxTimesteps_ = round(duration/dt);
      dt_ = duration/maxTimesteps_;
   }
   int maxTimesteps() const { return maxTimesteps_; };
   double dt() const { return dt_; }
   double maxTime() const { return dt_*maxTimesteps_; }
   double realTimeFromTimestep(int timestep) const
   {
      return timestep*dt_;
   }
   int timestepFromRealTime(double realTime) const
   {
      return round(realTime/dt_);
   }

 private:
   double dt_;
   int maxTimesteps_;
};

void  readGF(OBJECT* obj, const std::string dataname,
		mfem::Mesh* mesh, std::shared_ptr<mfem::GridFunction>& sp_gf) {
   std::string filename;
   objectGet(obj,dataname,filename," ");

   std::ifstream infile(filename);
   sp_gf = std::make_shared<mfem::GridFunction>(mesh, infile);
}

int main(int argc, char *argv[])
{
   MPI_Init(NULL, NULL);   // without specific command line
   int num_ranks, my_rank;
   MPI_Comm_size(COMM_LOCAL,&num_ranks); // number of processes in COMM_LOCAL
   MPI_Comm_rank(COMM_LOCAL,&my_rank); 

   units_internal(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);
   units_external(1e-3, 1e-9, 1e-3, 1e-3, 1, 1e-9, 1);

   std::string dataname = "femheart.data";
   std::cout<<"hhh"<<std::endl;
   std::cout<<"hhh"<<argc<<std::endl;
   std::cout <<"dataname1"<<dataname<<std::endl;

   OptionsParser args(argc, argv);
   args.AddOption(&dataname, "-m", "--mesh", "Data file to use.");
   args.ParseCheck();

   if (my_rank == 0)
   {  
      std::cout<< dataname << std::endl;
      std::cout << "Initializing with hhhh" << num_ranks << " MPI ranks." << std::endl;
   }
   
   int order = 1;

   std::vector<std::string> objectFilenames;
   //if (argc == 1)
   objectFilenames.push_back(dataname);

   //for (int iargCursor=1; iargCursor<argc; iargCursor++)
     // objectFilenames.push_back(argv[iargCursor]);

   if (my_rank == 0) {
      for (int ii=0; ii<objectFilenames.size(); ii++)
	 object_compilefile(objectFilenames[ii].c_str());
   }
   object_Bcast(0,MPI_COMM_WORLD);

   OBJECT* obj = object_find("femheart", "HEART");
   assert(obj != NULL);

   bool paraview;
   objectGet(obj, "paraview", paraview, "");

   bool glvis;
   objectGet(obj, "glvis", glvis, "");

   //StartTimer("Read the mesh");
   // Read shared global mesh
   std::string meshname;
   objectGet(obj, "mesh", meshname, "");
   mfem::Mesh *mesh = new mfem::Mesh(meshname.c_str(), 1, 1);  // generate_edges=1, refine=1 (fix_orientation=true by default)
   //mfem::Mesh *mesh = ecg_readMeshptr(obj, "mesh");
   //EndTimer();
   //Load fiber from file
   std::shared_ptr<GridFunction> flat_fiber;
   readGF(obj, "fiber", mesh, flat_fiber);

   std::shared_ptr<GridFunction> flat_sheet;
   readGF(obj, "sheet", mesh, flat_sheet);
 
   std::shared_ptr<GridFunction> flat_trans;
   readGF(obj, "trans", mesh, flat_trans);

   std::shared_ptr<GridFunction> flat_anatomy;
   readGF(obj, "anatomy", mesh, flat_anatomy);

   int dim = mesh->Dimension();

   //Fill in the MatrixElementPiecewiseCoefficients
   //std::vector<int> heartRegions;
   //objectGet(obj,"heart_regions", heartRegions);

   std::vector<double> sigma_m;
   objectGet(obj,"sigma_m",sigma_m);
   //assert(heartRegions.size()*3 == sigma_m.size());

   double dt;
   objectGet(obj,"dt",dt,"0.01 ms");
   double Bm;
   objectGet(obj,"Bm",Bm,"140"); // 1/mm
   double Cm;
   objectGet(obj,"Cm",Cm,"0.01"); // 1 uF/cm^2 = 0.01 uF/mm^2
 
   std::string reactionName;
   objectGet(obj, "reaction", reactionName, "BetterTT06");

   std::string outputDir;
   objectGet(obj, "outdir", outputDir, ".");
   
   double endTime;
   objectGet(obj, "end_time", endTime, "0 ms");

   double outputRate;
   objectGet(obj, "output_rate", outputRate, "1 ms");

   //double checkpointRate;
   //objectGet(obj, "checkpoint_rate", checkpointRate, "100 ms");

   double initVm;
   objectGet(obj, "init_vm", initVm, "-83");

   //bool useNodalIion = true;
   //objectGet(obj, "nodal_ion", useNodalIion, "1");

   StimulusCollection stims(dt);
   {
      std::vector<std::string> stimulusNames;
      objectGet(obj, "stimulus", stimulusNames);
      for (auto name : stimulusNames)
      {
         OBJECT* stimobj = object_find(name.c_str(), "STIMULUS");
         assert(stimobj != NULL);
         int numTimes;
         objectGet(stimobj, "n", numTimes, "1");
         double bcl;
         objectGet(stimobj, "bcl", bcl, "0 ms");
         assert(numTimes == 1 || bcl != 0);
         double startTime;
         objectGet(stimobj, "start", startTime, "0 ms");
         double duration;
         objectGet(stimobj, "duration", duration, "1 ms");
         double strength;
         objectGet(stimobj, "strength", strength, "0"); //uA/uF
         assert(strength >= 0);
         std::string location;
         objectGet(stimobj, "where", location, "");
         assert(!location.empty());
         OBJECT* locobj = object_find(location.c_str(), "REGION");
         assert(locobj != NULL);
         std::string regionType;
         objectGet(locobj, "type", regionType, "");
         assert(!regionType.empty());
         shared_ptr<StimulusLocation> stimLoc;
         if (regionType == "ball")
         {
            std::vector<double> center;
            objectGet(locobj, "center", center);
            assert(center.size() == 3);
            double radius;
            objectGet(locobj, "radius", radius, "-1");
            assert(radius >= 0);
            stimLoc = std::make_shared<CenterBallStimulus>(center[0],center[1],center[2],radius);
         }
         else if (regionType == "box")
         {
            std::vector<double> lower;
            objectGet(locobj, "lower", lower);
            assert(lower.size() == 3);
            vector<double> upper;
            objectGet(locobj, "upper", upper);
            assert(upper.size() == 3);
            stimLoc = std::make_shared<BoxStimulus>
               (lower[0], upper[0],
                lower[1], upper[1],
                lower[2], upper[2]);
         }
         shared_ptr<StimulusWaveform> stimWave(new SquareWaveform());
         stims.add(Stimulus(numTimes, startTime, duration, bcl, strength, stimLoc, stimWave));
      }
   }
   
   Timeline timeline(dt, endTime);   

   
   //StartTimer("Setting Attributes");
   //mesh->SetAttributes();
   //EndTimer();

   
  // StartTimer("Partition Mesh");
   // If I read correctly, pmeshpart will now point to an integer array
   //  containing a partition ID (rank!) for every element ID.
   int *pmeshpart = mesh->GeneratePartitioning(num_ranks);
  // EndTimer();


   ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   {
      int num_nodes = pmesh.GetNV();
      std::cout << "MPI Rank " << my_rank << " has "  << num_nodes << " nodes." << std::endl;
   }

   pmesh->Save("solution/mesh");

   std::shared_ptr<ParGridFunction> fiber;
   fiber = std::make_shared<mfem::ParGridFunction>(pmesh, flat_fiber.get(), pmeshpart);

   std::shared_ptr<ParGridFunction> sheet;
   sheet = std::make_shared<mfem::ParGridFunction>(pmesh, flat_sheet.get(), pmeshpart);

   std::shared_ptr<ParGridFunction> trans;
   trans = std::make_shared<mfem::ParGridFunction>(pmesh, flat_trans.get(), pmeshpart);

   std::shared_ptr<ParGridFunction> anatomy;
   anatomy = std::make_shared<mfem::ParGridFunction>(pmesh, flat_anatomy.get(), pmeshpart);
   
   // Build a new FEC...
   FiniteElementCollection *fec;
   if (my_rank == 0) { std::cout << "Creating new FEC..." << std::endl; }
   fec = new H1_FECollection(order, dim);
   // ...and corresponding FES
   ParFiniteElementSpace *pfespace = new ParFiniteElementSpace(pmesh, fec);
   std::cout << "[" << my_rank << "] Number of finite element unknowns: "
	     << pfespace->GetTrueVSize() << std::endl;

   // 5. Determine the list of true (i.e. conforming) essential boundary DOFs
   Array<int> ess_tdof_list;   // Essential true degrees of freedom
   // "true" takes into account shared vertices.
   {
      Array<int> ess_bdr(pmesh->bdr_attributes.Max());
      ess_bdr = 0;
      pfespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to pfespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   ParGridFunction gf_Vm(pfespace);
   gf_Vm = initVm;

   ParaViewDataCollection pd("V_m", pmesh);
   delete pmesh;
   pd.SetPrefixPath("ParaView");
   pd.RegisterField("solution", &gf_Vm);
   if(paraview){
      pd.SetLevelsOfDetail(order);
      pd.SetDataFormat(VTKFormat::BINARY);
      pd.SetHighOrderOutput(true);
      pd.SetCycle(0);
      pd.SetTime(0.0);
      pd.Save();
   }
   {
      char cwd[255];
      getcwd(cwd, sizeof(cwd));
      std::string folderPath = std::string(cwd) + "/" + "solution";
      mkdir(folderPath.c_str(), 0777);
   }
   if(glvis){
      gf_Vm.Save("solution/sol.0.gf");
   }
   
   // Load conductivity data
   MatrixElementPiecewiseCoefficient sigma_m_coeffs(fiber, sheet, trans);
   Vector sigma_m_vec(3);
   for (int jj=0; jj<3; jj++)
   {
      double value = sigma_m[jj]*dt/Bm/Cm;
      sigma_m_vec[jj] = value;
   }
   sigma_m_coeffs.heartConductivities_[1] = sigma_m_vec;
 

   // 8. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   // NOTICE THE FLIP IN SIGNS FOR SIGMA!  This is on purpose, Diffusion does -div(sigma*grad)

   //StartTimer("Forming bilinear system (RHS)");

   ConstantCoefficient one(1.0);
   //StartTimer("Forming bilinear system (LHS)");
   
   // Brought out of loop to avoid unnecessary duplication
   ParBilinearForm *a = new ParBilinearForm(pfespace);   // defines a.
   a->AddDomainIntegrator(new DiffusionIntegrator(sigma_m_coeffs)); //(Q grad u, grad v) Q means coeff
   a->AddDomainIntegrator(new MassIntegrator(one));//(Q u, v) Q means coeff
   a->Assemble();// Assembles the form i.e. sums over all domain/bdr integrators.
   HypreParMatrix LHS_mat;
   a->FormSystemMatrix(ess_tdof_list,LHS_mat);
   //EndTimer();

   ParBilinearForm *b = new ParBilinearForm(pfespace);
   b->AddDomainIntegrator(new MassIntegrator(one));
   b->Update(pfespace);  // Update the @a FiniteElementSpace and delete all data associated with the old one.
   b->Assemble();
   HypreParMatrix RHS_mat;
   b->FormSystemMatrix(ess_tdof_list, RHS_mat);
   //EndTimer();

   ParBilinearForm *Iion_blf = new ParBilinearForm(pfespace);
   ConstantCoefficient dt_coeff(dt);
   Iion_blf->AddDomainIntegrator(new MassIntegrator(dt_coeff));
   Iion_blf->Update(pfespace);
   Iion_blf->Assemble();
   HypreParMatrix Iion_mat;
   Iion_blf->FormSystemMatrix(ess_tdof_list,Iion_mat);


   //Set up the ionic models
   ParLinearForm *c = new ParLinearForm(pfespace);
   //positive dt here because the reaction models use dVm = -Iion
   c->AddDomainIntegrator(new DomainLFIntegrator(stims));

   //Set up the solve
   HyprePCG pcg(LHS_mat);
   pcg.SetTol(1e-12);
   pcg.SetMaxIter(2000);
   pcg.SetPrintLevel(2);
   HypreSolver *M_test = new HypreBoomerAMG(LHS_mat);
   pcg.SetPreconditioner(*M_test);
   
   ThreadServer& threadServer = ThreadServer::getInstance();
   ThreadTeam defaultGroup = threadServer.getThreadTeam(vector<unsigned>());

   std::vector<std::string> reactionNames;
   reactionNames.push_back(reactionName);

   std::vector<int> cellTypes(anatomy->Size());
   for (int i = 0; i <  anatomy->Size(); ++i) {
    cellTypes[i] = static_cast<int>((*anatomy)(i));
   }

   ReactionWrapper reactionWrapper(dt,reactionNames,defaultGroup,cellTypes); 

   {
      std::vector<int> types = reactionWrapper.allCellTypes();
      std::ofstream outputFile("type.txt");
      auto it1 = types.begin();
      while (it1 != types.end()){
            outputFile << *it1<<std::endl;
            ++ it1;      
         }
      outputFile.close();
   }
   {
      std::ofstream outputFile("handle.txt");
      auto it1 = reactionWrapper.handle.begin();
      while (it1 != reactionWrapper.handle.end()){
            outputFile << *it1<<std::endl;
            ++ it1;      
         }
      outputFile.close();
   }
   {
      std::ofstream outputFile("varNames.txt");
      auto it1 = reactionWrapper.Filednames.begin();
      while (it1 != reactionWrapper.Filednames.end()){
         outputFile << *it1 <<std::endl;
         ++ it1;
      }
      outputFile.close();
   }
   {
      std::ofstream outputFile("varUnits.txt");
      auto it1 = reactionWrapper.FieldUnits.begin();
      while (it1 != reactionWrapper.FieldUnits.end()){
         outputFile << *it1 <<std::endl;
         ++ it1;
      }
      outputFile.close();
   }

   reactionWrapper.Initialize();
   cellTypes.clear();
   reactionNames.clear();


   int pfespace_size = pfespace->GetTrueVSize();
   Vector actual_Vm(pfespace_size); 
   Vector actual_b(pfespace_size);
   Vector actual_old(pfespace_size);
   Vector actual_Iion(pfespace_size);


   actual_Vm = reactionWrapper.getVmReadonly();
   
   int itime = 0;
   int i_name = 1;
   clock_t time_start = clock();
   std::vector<vector<double>> data(reactionWrapper.handle.size()+1);
   while (itime != timeline.maxTimesteps())
   {  

      if (my_rank == 0)
      {  
         double time = (double)(clock()-time_start)/CLOCKS_PER_SEC;
         std::cout <<"times =" << time << "seconds." << std::endl;
         std::cout << "time = " << timeline.realTimeFromTimestep(itime) << std::endl;
      }

      reactionWrapper.getVmReadwrite() = actual_Vm; //should be a memcpy
      reactionWrapper.Calc(dt);
      
      //add stimulii
      stims.updateTime(timeline.realTimeFromTimestep(itime));
      
      //compute the Iion and stimulus contribution
      c->Update();
      c->Assemble();
      a->FormLinearSystem(ess_tdof_list, gf_Vm, *c, LHS_mat, actual_Vm, actual_b, 1);
      //compute the RHS matrix contribution
      RHS_mat.Mult(actual_Vm, actual_old);
      actual_b += actual_old;

      Iion_mat.Mult(reactionWrapper.getIionReadonly(), actual_old);
      actual_b += actual_old;

      //solve the matrix
      pcg.Mult(actual_b, actual_Vm);

      a->RecoverFEMSolution(actual_Vm, *c, gf_Vm);

      itime++;
      //output if appropriate
     if ((itime % timeline.timestepFromRealTime(outputRate)) == 0)
      {
         if(paraview){
            pd.SetCycle(itime);
            pd.SetTime(timeline.realTimeFromTimestep(itime));
            pd.Save();
         }
         if(glvis){

               std::string filename = "solution/sol." + std::to_string(i_name) + ".gf";
               
               gf_Vm.Save(filename.c_str());
               
         }
         reactionWrapper.GetValue(3);
         for (int ii = 0 ; ii < reactionWrapper.value.size(); ++ii)
            {
               data[ii].push_back(reactionWrapper.value[ii]);
            }
         data[reactionWrapper.value.size()].push_back(timeline.realTimeFromTimestep(itime));
         i_name++;
      }
   }


   std::string filename = "data.raw";
   std::ofstream file(filename, std::ios::binary);
   for(const auto& row : data)
   {
      file.write(reinterpret_cast<const char*>(row.data()),row.size() * sizeof(double));
   }
   file.close();


   // 14. Free the used memory.
   delete M_test;
   delete a;
   delete b;
   delete c;
   delete pfespace;
   if (order > 0) { delete fec; }
   MPI_Finalize();
   
   return 0;
}
