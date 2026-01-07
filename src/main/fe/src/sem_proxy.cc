//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "sem_proxy.h"
#include "ToPPMConverter.h"

#include <cartesian_struct_builder.h>
#include <cartesian_unstruct_builder.h>
#include <sem_solver_acoustic.h>
#include <source_and_receiver_utils.h>

#include <cxxopts.hpp>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <variant>

using namespace SourceAndReceiverUtils;

SEMproxy::SEMproxy(const SemProxyOptions& opt)
{
  const int order = opt.order;
  nb_elements_[0] = opt.ex;
  nb_elements_[1] = opt.ey;
  nb_elements_[2] = opt.ez;
  nb_nodes_[0] = opt.ex * order + 1;
  nb_nodes_[1] = opt.ey * order + 1;
  nb_nodes_[2] = opt.ez * order + 1;

  const float spongex = opt.boundaries_size;
  const float spongey = opt.boundaries_size;
  const float spongez = opt.boundaries_size;
  const std::array<float, 3> sponge_size = {spongex, spongey, spongez};
  src_coord_[0] = opt.srcx;
  src_coord_[1] = opt.srcy;
  src_coord_[2] = opt.srcz;

  domain_size_[0] = opt.lx;
  domain_size_[1] = opt.ly;
  domain_size_[2] = opt.lz;

  rcv_coord_[0] = opt.rcvx;
  rcv_coord_[1] = opt.rcvy;
  rcv_coord_[2] = opt.rcvz;

  bool isModelOnNodes = opt.isModelOnNodes;
  isElastic_ = opt.isElastic;
  cout << boolalpha;
  bool isElastic = isElastic_;

  const SolverFactory::methodType methodType = getMethod(opt.method);
  const SolverFactory::implemType implemType = getImplem(opt.implem);
  const SolverFactory::meshType meshType = getMesh(opt.mesh);
  const SolverFactory::modelLocationType modelLocation =
      isModelOnNodes ? SolverFactory::modelLocationType::OnNodes
                     : SolverFactory::modelLocationType::OnElements;
  const SolverFactory::physicType physicType = SolverFactory::physicType::Acoustic;

  float lx = domain_size_[0];
  float ly = domain_size_[1];
  float lz = domain_size_[2];
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  if (meshType == SolverFactory::Struct)
  {
    switch (order)
    {
      case 1: {
        model::CartesianStructBuilder<float, int, 1> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 2: {
        model::CartesianStructBuilder<float, int, 2> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      case 3: {
        model::CartesianStructBuilder<float, int, 3> builder(
            ex, lx, ey, ly, ez, lz, isModelOnNodes);
        m_mesh = builder.getModel();
        break;
      }
      default:
        throw std::runtime_error(
            "Order other than 1 2 3 is not supported (semproxy)");
    }
  }
  else if (meshType == SolverFactory::Unstruct)
  {
    model::CartesianParams<float, int> param(order, ex, ey, ez, lx, ly, lz,
                                             isModelOnNodes);
    model::CartesianUnstructBuilder<float, int> builder(param);
    m_mesh = builder.getModel();
  }
  else
  {
    throw std::runtime_error("Incorrect mesh type (SEMproxy ctor.)");
  }

  // Get unique data filder name

  namespace fs = std::filesystem;
  fs::path rootDir = "data";
  fs::create_directories(rootDir);


  time_t rawtime;
  struct tm * timeinfo;
  char buffer[100];
  time (&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,sizeof(buffer),"%d-%m-%Y-%H:%M:%S",timeinfo);
  date_ = buffer;


  std::string datedFolder = "data_" + date_;
  fs::path runDir = rootDir / datedFolder;
  fs::create_directories(runDir);

  fs::create_directories(runDir / "slices");
  fs::create_directories(runDir / "slice_images");
  fs::create_directories(runDir / "snapshots");
  fs::create_directories(runDir / "sismos");

  // save parameters
  is_snapshots_ = opt.isSnapshot;
  is_slices_ = opt.isSlice;
  snap_time_interval_ = opt.snapTimeInterval;
  is_in_situ = opt.isInSitu;
  data_folder_ = "data/data_" + date_ + "/";
  is_ppm_slices_ = opt.isPPM;
  plane_ = opt.plane;
  slice_position_[0] = opt.slice_posx;
  slice_position_[1] = opt.slice_posy;
  slice_position_[2] = opt.slice_posz;

  // time parameters
  if (opt.autodt)
  {
    float cfl_factor = (order == 2) ? 0.5 : 0.7;
    dt_ = find_cfl_dt(cfl_factor);
  }
  else
  {
    dt_ = opt.dt;
  }
  timemax_ = opt.timemax;
  num_sample_ = timemax_ / dt_;



  std::ifstream selectPointFile(opt.receiverfilename);

  if(selectPointFile.is_open()){
    std::string line;
    int x, y, z;
    int nodeIndex;
    while(std::getline(selectPointFile, line)){
      std::istringstream iss(line);
      if(!(iss >> x >> y >> z)){
        throw std::runtime_error("Error reading receiver coordinates from file.");
      }
      bool validPoint = false;
      for (nodeIndex = 0; nodeIndex < m_mesh->getNumberOfNodes(); nodeIndex++) {
        if(x == m_mesh->nodeCoord(nodeIndex, 0) && y == m_mesh->nodeCoord(nodeIndex, 1) && z == m_mesh->nodeCoord(nodeIndex, 2)){
          validPoint = true;
          break;
        }
      }
      if(!validPoint){
        throw std::runtime_error("Receiver coordinate does not match any mesh node.");
      }

      selectPoint.push_back(nodeIndex);
    }
    selectPointFile.close();


  }

  compression = opt.compression;

  m_solver = SolverFactory::createSolver(methodType, implemType, meshType,
                                         modelLocation, physicType, order);
  m_solver->computeFEInit(*m_mesh, sponge_size, opt.surface_sponge,
                          opt.taper_delta);

  initFiniteElem();

  std::cout << "Number of node is " << m_mesh->getNumberOfNodes() << std::endl;
  std::cout << "Number of element is " << m_mesh->getNumberOfElements()
            << std::endl;
  std::cout << "Launching the Method " << opt.method << ", the implementation "
            << opt.implem << " and the mesh is " << opt.mesh << std::endl;
  std::cout << "Model is on " << (isModelOnNodes ? "nodes" : "elements")
            << std::endl;
  std::cout << "Physics type is " << (isElastic ? "elastic" : "acoustic")
            << std::endl;
  std::cout << "Order of approximation will be " << order << std::endl;
  std::cout << "Time step is " << dt_ << "s" << std::endl;
  std::cout << "Simulated time is " << timemax_ << "s" << std::endl;
}

FILE* open_file(string filename){
  FILE *file = fopen(filename.c_str(), "a+");
  if (!file) {
    fprintf(stderr, "Couldn't open file %s\n", filename.c_str());
    exit(EXIT_FAILURE);
  }
  return file;
}

void SEMproxy::saveSismo(int timestep)
{
  string filename = data_folder_ + "sismos/sismo.csv";
  FILE *file = open_file(filename);

  if(timestep == 0){
    fprintf(file, "index timestep x y z pressure\n");
  }
  for (int i = 0; i < selectPoint.size(); i++) {
    fprintf(file, "%d %d %f %f %f %f\n",i, timestep,  m_mesh->nodeCoord(selectPoint[i], 0),  m_mesh->nodeCoord(selectPoint[i], 1), m_mesh->nodeCoord(selectPoint[i], 2), pnGlobal(selectPoint[i], i1));
  }
  fclose(file);
}

void SEMproxy::compresseRLESismo(int timestep)
{
  string filename = data_folder_ + "sismos/sismo.rle";
  FILE *file = open_file(filename);

  if(timestep == 0){
    fprintf(file, "pressure");
  }
  
  if(timestep < num_sample_){
    for (int i = 0; i < selectPoint.size(); i++) {
      if(accumulator == -1){
        prevPressure = pnGlobal(selectPoint[i], i1);
        accumulator = 1;
      }else{
        if(prevPressure == pnGlobal(selectPoint[i], i1)){
          accumulator++;
        }
        else{
          if(accumulator == 1){
            fprintf(file, " %f",prevPressure);
          }
          else{
            fprintf(file," %dx%f",accumulator,prevPressure);
          }
          accumulator = 1;
          prevPressure = pnGlobal(selectPoint[i], i1);
        }
      }
      
    }
  }
  if(timestep == num_sample_ - 1){
    //pressure
    if(accumulator == 1){
      fprintf(file, " %f",prevPressure);
    }
    else{
      fprintf(file," %dx%f",accumulator,prevPressure);
    }
    fprintf(file, "\n");

    // z
    fprintf(file,"z %dx(",num_sample_);
    for (int i = 0; i < selectPoint.size(); i++) {
      if(i != 0){
        fprintf(file, ",");
      }
      fprintf(file, "%f",m_mesh->nodeCoord(selectPoint[i], 2));
    }
    fprintf(file, ")\n");

    // y
    fprintf(file,"y %dx(",num_sample_);
    for (int i = 0; i < selectPoint.size(); i++) {
      if(i != 0){
        fprintf(file, ",");
      }
      fprintf(file, "%f",m_mesh->nodeCoord(selectPoint[i], 1));
    }
    fprintf(file, ")\n");

    // x
    fprintf(file,"x %dx(",num_sample_);
    for (int i = 0; i < selectPoint.size(); i++) {
      if(i != 0){
        fprintf(file, ",");
      }
      fprintf(file, "%f",m_mesh->nodeCoord(selectPoint[i], 0));
    }
    fprintf(file, ")\n");

    // timestep
    fprintf(file,"timestep");
    for (int i = 0; i < num_sample_; i++) {
      fprintf(file," %ldx%d",selectPoint.size(),timestep);
    }
    fprintf(file,"\n");

    // index
    fprintf(file,"index %dx(",num_sample_);
    for (int i = 0; i < selectPoint.size(); i++) {
      if(i != 0){
        fprintf(file, ",");
      }
      fprintf(file, "%d",i);
    }
    fprintf(file, ")\n");
  }
  fclose(file);
}

void SEMproxy::saveMeasure(float kerneltime_ms, float outputtime_ms, float traitementtime_ms = 0.0f) {
  string filename = data_folder_ + "measure.csv";
  FILE *file = open_file(filename);

  long int sizefile_snapshots = 0;
  if(is_snapshots_){
    for(int i = 0; i < num_sample_/snap_time_interval_; i++){
      std::string snapshotfile = data_folder_ + "snapshots/snapshot_" + to_string(i) + ".csv";
      FILE *snapshot = open_file(snapshotfile);
      fseek(snapshot, 0, SEEK_END);
      sizefile_snapshots += ftell(snapshot);
      fclose(snapshot);
    }
  }

  long int sizefile_slices = 0;
  if(is_slices_){
    for(int i = 0; i < num_sample_/snap_time_interval_; i++){
      std::string slicefile = data_folder_ + "slices/slice_" + to_string(i) + ".csv";
      FILE *slice = open_file(slicefile);
      fseek(slice, 0, SEEK_END);
      sizefile_slices += ftell(slice);
      fclose(slice);
    }
  }

  long int sizefile_sismos = 0;
  if(selectPoint.size() > 0){
    std::string sismofile = data_folder_ + "sismos/sismo.csv";
    FILE *sismo = open_file(sismofile);
    fseek(sismo, 0, SEEK_END);
    sizefile_sismos += ftell(sismo);
    fclose(sismo);
  }
  
  fseek(file, 0, SEEK_END);
  if (ftell(file) == 0) {
    fprintf(file, "kernel_time output_time traitement_time size_file_snapshots size_file_slices size_file_sismos\n");
  }

  fprintf(file, "%f %f %f %ld %ld %ld\n", kerneltime_ms, outputtime_ms, traitementtime_ms, sizefile_snapshots, sizefile_slices, sizefile_sismos);

  fclose(file);
}

void SEMproxy::saveAnalyse(int index, float min, float max, float mean, float median, float std ,const char* analysisType ) {
  string filename = data_folder_ + analysisType + "_" +  "analysis.csv";
  FILE *file = open_file(filename);

  fseek(file, 0, SEEK_END);
  if (ftell(file) == 0) {
    fprintf(file, "index min max mean median std\n");
  }

  fprintf(file, "%d %f %f %f %f %f\n",index, min, max, mean, median, std );

  fclose(file);
}

void SEMproxy::saveSliceAsPPM(int timestep) {
  const int slice_num = timestep / snap_time_interval_;
  std::string filename = data_folder_ + "slice_images/slice_"
                     + to_string(slice_num) + ".ppm";
  float node_size_x = floor(domain_size_[0] / (nb_nodes_[0] - 1));
  float node_size_y = floor(domain_size_[1] / (nb_nodes_[1] - 1));
  float node_size_z = floor(domain_size_[2] / (nb_nodes_[2] - 1));
  
  std::array<float, 3UL> plane_src;
  float srcx = (floor(slice_position_[0] / node_size_x) + 1) * node_size_x;
  float srcy = (floor(slice_position_[1] / node_size_y)  + 1) * node_size_y;
  float srcz = (floor(slice_position_[2] / node_size_z)  + 1) * node_size_z;
  std::array<float, 3UL> src_position;
  src_position[0] = srcx;
  src_position[1] = srcy;
  src_position[2] = srcz;
  ToPPMConverter::convert(filename, m_mesh, pnGlobal, i1, plane_, src_position, nb_nodes_);
}

void SEMproxy::saveSlice(int timestep) {
  const int slice_num = timestep / snap_time_interval_;
  std::string filename = data_folder_ + "slices/slice_"
                     + to_string(slice_num) + ".csv";
  FILE *file = open_file(filename);

  fprintf(file, "plane timestep i j pressure\n");

  float node_size_x = floor(domain_size_[0] / (nb_nodes_[0] - 1));
  float node_size_y = floor(domain_size_[1] / (nb_nodes_[1] - 1));
  float node_size_z = floor(domain_size_[2] / (nb_nodes_[2] - 1));

  float srcx = (floor(src_coord_[0] / node_size_x) + 1) * node_size_x;
  float srcy = (floor(src_coord_[1] / node_size_y)  + 1) * node_size_y;
  float srcz = (floor(src_coord_[2] / node_size_z)  + 1) * node_size_z;

  for (int nodeIndex = 0; nodeIndex < m_mesh->getNumberOfNodes(); nodeIndex++) {
    float x = m_mesh->nodeCoord(nodeIndex, 0);
    float y = m_mesh->nodeCoord(nodeIndex, 1);
    float z = m_mesh->nodeCoord(nodeIndex, 2);

    if (z == srcz) {
      fprintf(file, "xy %d %f %f %f\n", timestep, x, y, pnGlobal(nodeIndex, i1));
    }
    if (y == srcy) {
      fprintf(file, "xz %d %f %f %f\n", timestep, x, z, pnGlobal(nodeIndex, i1));
    }
    if (x == srcx) {
      fprintf(file, "yz %d %f %f %f\n", timestep, y, z, pnGlobal(nodeIndex, i1));
    }
  }

  fclose(file);
}

void SEMproxy::saveSnapshot(int timestep) {
  const int snapshot_num = timestep / snap_time_interval_;
  std::string filename = data_folder_ + "snapshots/snapshot_"
                    + to_string(snapshot_num) + ".csv";
  FILE *file = open_file(filename);

  fprintf(file, "snap timestep x y z pressure\n");
  const int order = m_mesh->getOrder();
  for (int nodeIndex = 0; nodeIndex < m_mesh->getNumberOfNodes(); nodeIndex++) {
    float x = m_mesh->nodeCoord(nodeIndex, 0);
    float y = m_mesh->nodeCoord(nodeIndex, 1);
    float z = m_mesh->nodeCoord(nodeIndex, 2);

    float pressure = pnGlobal(nodeIndex, i1);
    fprintf(file, "%d %d %f %f %f %f\n", snapshot_num, timestep, x, y, z, pressure);
  }
  fclose(file);
}

void SEMproxy::run()
{
  float maxPressurePerReceive[selectPoint.size()] = {-999};
  float minPressurePerReceive[selectPoint.size()] = {999};
  float meanPressurePerReceive[selectPoint.size()] = {0};

  float sd_sommePerReceive[selectPoint.size()] = {0};
  
  float sismos[selectPoint.size()][num_sample_];
  time_point<system_clock> startComputeTime, startOutputTime,startTraitementTime, totalTraitementTime, totalComputeTime,
      totalOutputTimeOneStep,totalOutputTime;

  SEMsolverDataAcoustic solverData(i1, i2, myRHSTerm, pnGlobal, rhsElement,
                                   rhsWeights);

  for (int indexTimeSample = 0; indexTimeSample < num_sample_;
       indexTimeSample++)
  {
    startComputeTime = system_clock::now();
    m_solver->computeOneStep(dt_, indexTimeSample, solverData);
    totalComputeTime += system_clock::now() - startComputeTime;
    
    if(is_in_situ){
      if (indexTimeSample % snap_time_interval_ == 0 && is_snapshots_){

          startTraitementTime = system_clock::now();

          float max = -999;
          float min = 999;
          float mean = 0;

          float list[m_mesh->getNumberOfNodes()];
          float median = 0;
          for (int nodeIndex = 0; nodeIndex < m_mesh->getNumberOfNodes(); nodeIndex++) {
            float pressure = pnGlobal(nodeIndex, i1);
            if(pressure > max){
              max = pressure;
            }
            if(pressure < min){
              min = pressure;
            }
              mean += pressure;
            
            list[nodeIndex] = pressure;
          }

          std::sort(list, list + m_mesh->getNumberOfNodes());

          median = list[m_mesh->getNumberOfNodes()/2];

          float sd_somme = 0;
          for (int nodeIndex = 0; nodeIndex < m_mesh->getNumberOfNodes(); nodeIndex++) {
            float pressure = pnGlobal(nodeIndex, i1);

            sd_somme += (pressure - mean) * (pressure - mean);
          }


          totalTraitementTime += system_clock::now() - startTraitementTime;
          
          startOutputTime = system_clock::now();

          saveAnalyse(indexTimeSample / snap_time_interval_,min,max,mean/m_mesh->getNumberOfNodes(),median,sqrt(sd_somme/m_mesh->getNumberOfNodes()),(std::string("snapshot")).c_str());

          totalOutputTime += system_clock::now() - startOutputTime;
      }

      if(selectPoint.size() > 0){

        startTraitementTime = system_clock::now();

        for (int i = 0; i < selectPoint.size(); i++) {

          float pressure = pnGlobal(selectPoint[i], i1);

          sismos[i][indexTimeSample] = pressure;

          if(pressure > maxPressurePerReceive[i]){
            maxPressurePerReceive[i] = pressure;
          }

          if(pressure < minPressurePerReceive[i]){
            minPressurePerReceive[i] = pressure;
          }

          meanPressurePerReceive[i] += pressure;
        }

        totalTraitementTime += system_clock::now() - startTraitementTime;
      }
    }

    startOutputTime = system_clock::now();

    if (indexTimeSample % snap_time_interval_ == 0)
    {
      m_solver->outputSolutionValues(indexTimeSample, i1, rhsElement[0],
                                     pnGlobal, "pnGlobal");
      if (is_snapshots_)
        saveSnapshot(indexTimeSample);    
      if (is_slices_)
        saveSlice(indexTimeSample);
      if (is_ppm_slices_)
        saveSliceAsPPM(indexTimeSample);
    }

    if(selectPoint.size() > 0){
      if(compression){
        compresseRLESismo(indexTimeSample);
      }
      else{
        saveSismo(indexTimeSample);
      }
    }

    // Save pressure at receiver
    const int order = m_mesh->getOrder();

    float varnp1 = 0.0;
    for (int i = 0; i < order + 1; i++)
    {
      for (int j = 0; j < order + 1; j++)
      {
        for (int k = 0; k < order + 1; k++)
        {
          int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
          int globalNodeOnElement =
              i + j * (order + 1) + k * (order + 1) * (order + 1);
          varnp1 +=
              pnGlobal(nodeIdx, i2) * rhsWeightsRcv(0, globalNodeOnElement);
        }
      }
    }

    pnAtReceiver(0, indexTimeSample) = varnp1;

    swap(i1, i2);

    auto tmp = solverData.m_i1;
    solverData.m_i1 = solverData.m_i2;
    solverData.m_i2 = tmp;
    
    totalOutputTime += system_clock::now() - startOutputTime;

  }

  if(is_in_situ){
      if(selectPoint.size() > 0){

        startTraitementTime = system_clock::now();

        float medianes [selectPoint.size()];


        float sd_somme[selectPoint.size()] = {0};
        for(int j = 0; j < selectPoint.size();j++){
          for(int i = 0; i < num_sample_;i++){
            sd_somme[j] += (sismos[j][i] - (meanPressurePerReceive[j]/num_sample_)) * (sismos[j][i] - (meanPressurePerReceive[j]/num_sample_));
          }
          std::sort(sismos[j], sismos[j] + num_sample_);
          medianes[j] = sismos[j][num_sample_/2];
        }

        totalTraitementTime += system_clock::now() - startTraitementTime;

        startOutputTime = system_clock::now();

        for(int i = 0; i < selectPoint.size(); i++){
          saveAnalyse(i,minPressurePerReceive[i],minPressurePerReceive[i],meanPressurePerReceive[i]/num_sample_,medianes[i],sqrt(sd_somme[i]/num_sample_),(std::string("receiver")).c_str());
        }

        totalOutputTime += system_clock::now() - startOutputTime;
    }
  }

  float kerneltime_ms = 
      time_point_cast<microseconds>(totalComputeTime).time_since_epoch().count();

  float outputtime_ms =
      time_point_cast<microseconds>(totalOutputTime).time_since_epoch().count();

  if(!is_in_situ){
    saveMeasure(kerneltime_ms,outputtime_ms);
  }
  else{
    float traitementtime_ms =
        time_point_cast<microseconds>(totalTraitementTime).time_since_epoch().count();
    saveMeasure(kerneltime_ms,outputtime_ms,traitementtime_ms);
  }


  cout << "------------------------------------------------ " << endl;
  cout << "\n---- Elapsed Kernel Time : " << kerneltime_ms / 1E6 << " seconds."
       << endl;
  cout << "---- Elapsed Output Time : " << outputtime_ms / 1E6 << " seconds."
       << endl;
  cout << "------------------------------------------------ " << endl;
}

// Initialize arrays
void SEMproxy::init_arrays()
{
  cout << "Allocate host memory for source and pressure values ..." << endl;

  rhsElement = allocateVector<vectorInt>(myNumberOfRHS, "rhsElement");
  rhsWeights = allocateArray2D<arrayReal>(
      myNumberOfRHS, m_mesh->getNumberOfPointsPerElement(), "RHSWeight");
  myRHSTerm = allocateArray2D<arrayReal>(myNumberOfRHS, num_sample_, "RHSTerm");
  pnGlobal =
      allocateArray2D<arrayReal>(m_mesh->getNumberOfNodes(), 2, "pnGlobal");
  pnAtReceiver = allocateArray2D<arrayReal>(1, num_sample_, "pnAtReceiver");
  // Receiver
  rhsElementRcv = allocateVector<vectorInt>(1, "rhsElementRcv");
  rhsWeightsRcv = allocateArray2D<arrayReal>(
      1, m_mesh->getNumberOfPointsPerElement(), "RHSWeightRcv");
}

// Initialize sources
void SEMproxy::init_source()
{
  arrayReal myRHSLocation = allocateArray2D<arrayReal>(1, 3, "RHSLocation");
  // std::cout << "All source are currently are coded on element 50." <<
  // std::endl;
  std::cout << "All source are currently are coded on middle element."
            << std::endl;
  int ex = nb_elements_[0];
  int ey = nb_elements_[1];
  int ez = nb_elements_[2];

  int lx = domain_size_[0];
  int ly = domain_size_[1];
  int lz = domain_size_[2];

  // Get source element index

  int source_index = floor((src_coord_[0] * ex) / lx) +
                     floor((src_coord_[1] * ey) / ly) * ex +
                     floor((src_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElement[i] = source_index;
  }

  // Get coordinates of the corners of the sourc element
  float cornerCoords[8][3];
  int I = 0;
  int nodes_corner[2] = {0, m_mesh->getOrder()};
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElement[0], i, j, k);
        cornerCoords[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoords[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoords[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  // initialize source term
  vector<float> sourceTerm =
      myUtils.computeSourceTerm(num_sample_, dt_, f0, sourceOrder);
  for (int j = 0; j < num_sample_; j++)
  {
    myRHSTerm(0, j) = sourceTerm[j];
    if (j % 100 == 0)
      cout << "Sample " << j << "\t: sourceTerm = " << sourceTerm[j] << endl;
  }

  // get element number of source term
  myElementSource = rhsElement[0];
  cout << "Element number for the source location: " << myElementSource << endl
       << endl;

  int order = m_mesh->getOrder();

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoords, src_coord_,
                                                   rhsWeights);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }

  // Receiver computation
  int receiver_index = floor((rcv_coord_[0] * ex) / lx) +
                       floor((rcv_coord_[1] * ey) / ly) * ex +
                       floor((rcv_coord_[2] * ez) / lz) * ey * ex;

  for (int i = 0; i < 1; i++)
  {
    rhsElementRcv[i] = receiver_index;
  }

  // Get coordinates of the corners of the receiver element
  float cornerCoordsRcv[8][3];
  I = 0;
  for (int k : nodes_corner)
  {
    for (int j : nodes_corner)
    {
      for (int i : nodes_corner)
      {
        int nodeIdx = m_mesh->globalNodeIndex(rhsElementRcv[0], i, j, k);
        cornerCoordsRcv[I][0] = m_mesh->nodeCoord(nodeIdx, 0);
        cornerCoordsRcv[I][2] = m_mesh->nodeCoord(nodeIdx, 2);
        cornerCoordsRcv[I][1] = m_mesh->nodeCoord(nodeIdx, 1);
        I++;
      }
    }
  }

  switch (order)
  {
    case 1:
      SourceAndReceiverUtils::ComputeRHSWeights<1>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    case 2:
      SourceAndReceiverUtils::ComputeRHSWeights<2>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    case 3:
      SourceAndReceiverUtils::ComputeRHSWeights<3>(cornerCoordsRcv, rcv_coord_,
                                                   rhsWeightsRcv);
      break;
    default:
      throw std::runtime_error("Unsupported order: " + std::to_string(order));
  }
}

SolverFactory::implemType SEMproxy::getImplem(string implemArg)
{
  if (implemArg == "makutu") return SolverFactory::MAKUTU;
  if (implemArg == "shiva") return SolverFactory::SHIVA;

  throw std::invalid_argument(
      "Implentation type does not follow any valid type.");
}

SolverFactory::meshType SEMproxy::getMesh(string meshArg)
{
  if (meshArg == "cartesian") return SolverFactory::Struct;
  if (meshArg == "ucartesian") return SolverFactory::Unstruct;

  std::cout << "Mesh type found is " << meshArg << std::endl;
  throw std::invalid_argument("Mesh type does not follow any valid type.");
}

SolverFactory::methodType SEMproxy::getMethod(string methodArg)
{
  if (methodArg == "sem") return SolverFactory::SEM;
  if (methodArg == "dg") return SolverFactory::DG;

  throw std::invalid_argument("Method type does not follow any valid type.");
}

float SEMproxy::find_cfl_dt(float cfl_factor)
{
  float sqrtDim3 = 1.73;  // to change for 2d
  float min_spacing = m_mesh->getMinSpacing();
  float v_max = m_mesh->getMaxSpeed();

  float dt = cfl_factor * min_spacing / (sqrtDim3 * v_max);

  return dt;
}
