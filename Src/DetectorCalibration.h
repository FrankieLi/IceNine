//============================================================================== 
// Copyright (c) 2014, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by S. F. Li (li31@llnl.gov)
// LLNL-CODE-657639
// All rights reserved.
//
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//    * Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//    * Neither the name of the Lawrence Livermore National Lab nor the
//      names of its contributors may be used to endorse or promote products
//      derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL LAB BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//============================================================================== 

//------------------------------------------------------------------------------------
//  Author:  S. F. Li (Frankie)
//  e-mail:  li31@llnl.gov; sfli@cmu.edu 
//------------------------------------------------------------------------------------
////////////////////////////////////////////////////////////
//
//
//   DetectorCalibration.h
//   Author:   Frankie Li  (li31@llnl.gov)
//
//   Purpose:  New, refactor detector calibration class.
//             
//
//
////////////////////////////////////////////////////////////
#ifndef DETECTOR_CALIBRATION_H_
#define DETECTOR_CALIBRATION_H_


#include "Reconstructor.h"
#include "DiscreteAdaptive.h"
#include <ctime>
#include <boost/shared_ptr.hpp>
#include "GrainReconstruction.h"
#include "Simulation.h"
#include "nlopt.hpp"
#include <ostream>


namespace Calibration
{

  using namespace Reconstruction;
  //---------------------------------------------------
  //  DetectorCost
  //---------------------------------------------------
  class SingleDetectorCost
  {
  public:
    typedef SVoxel SamplePointT;
    typedef boost::shared_ptr< GrainCostFunction >            GrainCostPtr;
    typedef boost::shared_ptr< std::vector< SamplePointT > >  SamplePointPtr;
    
    
    //private:
  public:
    GrainCostPtr   pCostCalculator;
    SamplePointPtr pSamplePointList;
    
    CDetector     *pDetector;  // local pointer referring to the detector to be optimized     
    
    //----------------------------------------
    //
    //----------------------------------------
    void SetDetectorParameters( const std::vector<double> & Params )
    {
      CXDMDetectorFactory::SDetParameters DetParams 
	= CXDMDetectorFactory::GetImageParameters( *pDetector );;
      
      GeneralLib::SQuaternion q;
      GeneralLib::SVector3 RF;

      RF.m_fX = Params[0];
      RF.m_fY = Params[1];
      RF.m_fZ = Params[2];

      q.CreateFromRodrigues( RF );
      DetParams.oOrientation = q.GetRotationMatrix3x3();
      
      DetParams.oPosition.m_fX = Params[3];
      DetParams.fBeamCenterJ   = Params[4];
      DetParams.fBeamCenterK   = Params[5];
      
      CXDMDetectorFactory::ModifyImageParameters( *pDetector,
						  DetParams.oPosition,
						  DetParams.fBeamCenterJ, DetParams.fBeamCenterK,
						  DetParams.fPixelWidth, DetParams.fPixelHeight,
						  DetParams.oOrientation );
    }


    SingleDetectorCost();
 
    std::vector<double> StepSize;
  public:

    SingleDetectorCost( GrainCostPtr _pCostCalculator,
			SamplePointPtr _pSamplePointList,
			int DetectorToOptimize ) :
      pCostCalculator( _pCostCalculator ),
      pSamplePointList( _pSamplePointList )
    {
      RUNTIME_ASSERT( DetectorToOptimize 
		      < pCostCalculator->GetSetup()->ExperimentalSetup().GetDetectorList().size(),
		      " SingleDetectorCost: trying to optimize detector that does not exist\n" );
      pDetector = pCostCalculator->GetSetup()->ExperimentalSetup().GetDetector( DetectorToOptimize ); 
    }


    void SetStepSize( const std::vector<double> & StepSize_ )
    {
      StepSize = StepSize_;
    }

    
    //-------------------------------------------
    //-------------------------------------------
    Float operator() ( const std::vector<double> &v )
    {
      SetDetectorParameters( v );
      return pCostCalculator->GetAverageVoxelOverlapRatio( *pSamplePointList );
    }

    //-------------------------------------------
    //-------------------------------------------
    Float operator() ( const std::vector<double> &v, std::vector<double> & grad )
    {
      if( ! grad.empty() )
	grad = GetGrad( v, StepSize );
      
      
      SetDetectorParameters( v );
      return pCostCalculator->GetAverageVoxelOverlapRatio( *pSamplePointList );
    }
    
    //----------------------------------------
    //
    //----------------------------------------
    std::vector<double> GetGrad( const std::vector<double> & v,
				 const std::vector<double> & steps )
    {
      std::vector<double> g( v.size() );
      
      for( int i = 0; i < v.size(); i ++)
      {
	std::vector<double> u = v;
	u[i] = v[i] + steps[i];
	double forward = operator()( u );

	u[i] = v[i] - steps[i];
	double backward = operator()( u );
	
	g[i] = ( forward - backward ) / ( double(2) * steps[i] );
      }
      
      return g;
    }

    
    //----------------------------------------
    //
    //----------------------------------------
    std::vector<double> GetDetectorParameters( )
    { 
      CXDMDetectorFactory::SDetParameters DetParams 
	= CXDMDetectorFactory::GetImageParameters( *pDetector );
      
      GeneralLib::SQuaternion q;
      q.Set( DetParams.oOrientation );
      GeneralLib::SVector3 RF = q.GetRodriguesVector();
      
      std::vector<double> v(6);

      v[0] = RF.m_fX;
      v[1] = RF.m_fY;
      v[2] = RF.m_fZ;
      
      v[3] = DetParams.oPosition.m_fX;
      v[4] = DetParams.fBeamCenterJ;
      v[5] = DetParams.fBeamCenterK;
      
      return v;
    }
    
  };
  

  //---------------------------------------------------------------
  //---------------------------------------------------------------
  static double fn_wrapper(const std::vector<double> &x, std::vector<double> &grad, void *data)
  {
    return (*reinterpret_cast<SingleDetectorCost*>(data))(x, grad); 
  }
  
  //-------------------------------------------
  //  class DetectorCalibration
  //        Given a region of interest, auto-magically figure out
  //        the geometric parameters of the detector based on
  //        some previously tested heuristics. 
  //-------------------------------------------
  class DetectorCalibration
  {

  private:
    typedef SVoxel SamplePointT;

    
  public:
    
    DetectorCalibration() {}
    
    
    //-------------------------------
    //  CalibrateDetector
    //    - calibrate a single detector
    //-------------------------------
    void CalibrateDetector( )
    {
    }
    
  };


  struct Null_Deleter { void operator()( void const*) const{} };

  //-------------------------------------------
  //  
  //
  //-------------------------------------------
  class TestCostFunction
  {
  private:
    
    typedef SVoxel SamplePointT;
//    GrainCostFunction CostEval;
    

    std::vector<SamplePointT> SamplePoints;
  
    CConfigFile ConfigFile;

  public:
    
    TestCostFunction( const CConfigFile & _ConfigFile ):ConfigFile( _ConfigFile ){}
    

    //--------------------------------------------
    //
    //--------------------------------------------
    void TestCostVariation( const std::string & OutputFilename )
    {
    }


    
    //--------------------------------------------
    //
    //--------------------------------------------
    Float DrawCost()
    {
      typedef MicFile<SamplePointT>                                 Mic;      

      ReconstructionSetup     oSetup;
      CSimulation             oSimulator;      
      oSetup.InitializeWithDataFiles( ConfigFile );
      boost::shared_ptr<Mic> pMic= boost::dynamic_pointer_cast<Mic>( oSetup.ReconstructionRegion());
      oSimulator.Initialize( oSetup.ExperimentalSetup() );
      
      std::vector<SVoxel> VoxelList = pMic->GetVoxels();
      
      boost::shared_ptr< CSimulation >         pSimulator ( &oSimulator, Null_Deleter() );
      boost::shared_ptr< ReconstructionSetup > pSetup     ( &oSetup, Null_Deleter() );
      boost::shared_ptr<GrainCostFunction>     pCost( new GrainCostFunction( pSetup, pSimulator ) );
      boost::shared_ptr< std::vector<SVoxel> > VoxelListPtr( & VoxelList, Null_Deleter() );

      
      SingleDetectorCost DetCost( pCost, VoxelListPtr, 0 );
      
      //---------------------------------------------------------------------------------
      //  Selecting subset of reflection vectors
      CUnitCell * pUnitCell = oSetup.SampleGeometry().GetStructure( VoxelList[0].nPhase );
      std::vector<CRecpVector> RecpVectors = pUnitCell->GetReflectionVectorList();
      
      std::cout << " started with " << RecpVectors.size() << " recp vectors " << std::endl;
      
      std::random_shuffle( RecpVectors.begin(), RecpVectors.end() );
      
      int NumRef = std::min( int( RecpVectors.size() ), 20 );
      
      std::vector<CRecpVector> RecpVectorSubset;
      RecpVectorSubset.insert( RecpVectorSubset.begin(), 
			       RecpVectors.begin(), RecpVectors.begin() + NumRef );

      pUnitCell->SetReflectionVectorList( RecpVectorSubset );
      //---------------------------------------------------------------------------------
      
      std::ofstream os( "Test.txt" );
      std::vector<double> v = DetCost.GetDetectorParameters( );
      
      std::vector<double> DerStepSize(6);
      
      DerStepSize[0] = 0.002;
      DerStepSize[1] = 0.002;
      DerStepSize[2] = 0.002;
      
      DerStepSize[3] = 0.005;
      DerStepSize[4] = 0.02;
      DerStepSize[5] = 0.02;

      DetCost.SetStepSize( DerStepSize );
      
      std::cout << "Initial cost before pertubation: " 
		<< pCost->GetAverageVoxelOverlapRatio( VoxelList ) << " " << DetCost( v ) << std::endl;
      for( int j = 0; j < v.size(); j ++ )
	std::cout << v[j] << " ";
      std::cout << std::endl;
      
      
      for( int i = -50; i < 50; i ++ )   // verying only two parameters
      {
	//----------------------------------------------
	//  perform optimization
	//----------------------------------------------
	time_t oStartTime, oStopTime;
	struct tm * pTimeObj;
	time ( &oStartTime );
	pTimeObj = localtime ( &oStartTime );

//------------------------------------------
	nlopt::opt opt( nlopt::LN_NELDERMEAD, v.size() );
//	nlopt::opt opt( nlopt::LD_LBFGS, v.size() );
	opt.set_max_objective( fn_wrapper, & DetCost );

	std::vector<double> lower_bound(6);
	std::vector<double> upper_bound(6);
	
	for( int j = 0; j < v.size(); j ++ )
	{
	  lower_bound[j] = v[j] * 0.9; 
	  upper_bound[j] = v[j] * 1.1; 
	}

	for( int j = 0; j < 3; j ++ )
	{
	  lower_bound[j] = v[j] * 0.99;
	  upper_bound[j] = v[j] * 1.01;
	}
	
	
	opt.set_lower_bounds( lower_bound );
	opt.set_upper_bounds( upper_bound );
	
	std::vector<double> Steps(6);
	for( int j = 0; j < Steps.size(); j ++ )
	  Steps[j] = ( upper_bound[j] - lower_bound[j] ) / double(202);
	

	std::vector<double> x = v;
	for( int j = 3; j < v.size(); j ++ )
	  x[j] = v[j] + double( i ) * Steps[j];
	

	

	std::vector<double> u = x;
    	opt.set_xtol_rel(1e-4);

	double minf;
	nlopt::result result = opt.optimize( x, minf );
//------------------------------------------
	time ( &oStopTime );

	double diff = 0;
	for( int j = 0; j < v.size(); j ++ )
	{
	  //  std::cout << u[j] << " " << x[j] << std::endl;
	  diff += std::fabs( x[j] - v[j] ) / v[j] ;
	}
	diff /= double( v.size() );



	os << DetCost( u ) << " " << DetCost( x ) 
	   << " " 
	   << ( 1.0 + 0.001 * Float( i ) )
	   << " " << diff
	   << " " << difftime( oStopTime, oStartTime ) << std::endl;
      }
      
      return 0;
    }
   
  };



}


#endif
