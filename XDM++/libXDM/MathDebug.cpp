/////////////////////////////////////////////////////////////////
//
//  File:    MathDebug.cpp
//  Author:  Frankie Li
//  e-mail:  sfli@cmu.edu
//
//
/////////////////////////////////////////////////////////////////
#include "MathDebug.h"


void PrintVector3(const SVector3 &c)
{
	cout << c.m_fX << " " << c.m_fY << " " << c.m_fZ << endl;
}


void PrintVector4(const SVector4 &c)
{
	cout << c.m_fX << " " 
			 << c.m_fY << " " 
			 << c.m_fZ << " "
			 << c.m_fW << " " << endl;
}

void PrintMatrix4x4(const SMatrix4x4 &m)
{
	for (UInt i = 0; i < 4; i ++){
		for (UInt j = 0; j < 4; j ++){
			cout << m.m[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void PrintRay (const CRay &c)
{


	
	const SVector3 &oStart = c.GetStart();
	const SVector3 &oDir  = c.GetDirection();
	
	
	cout << "reflected ray: " << endl; 
	cout << "oStart: " << oStart.m_fX << " " << oStart.m_fY 
			 << " " << oStart.m_fZ << endl;
	
	cout << "oDir: " << oDir.m_fX << " " << oDir.m_fY 
			 << " " << oDir.m_fZ << endl;
			
						
}
