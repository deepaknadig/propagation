/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2013, SOLUTT Corporation
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Author: Deepak Nadig Anantha <deepak@solutt.com>
 *         Kushal S P           <kushalsp007@gmail.com>
 */

/*
 * SUI propagation path loss model Logic
 * =====================================
 * The SUI model covers three terrain categories.
 * Category A represents the maximum path-loss category which is a hill terrain,
 * Category B represents an intermediate path-loss category, and
 * Category C represents the minimum path-loss category with mostly flat terrains.
 *
 * The empirical formulas for this model are as below:
 * The median path-loss for the SUI model can be generally written as
 *                   (1)  PLsui = A + 10*gamma*log10(d/d0) + s ;
 *
 * 
 * for d>do, where do=100m. The term A (Intercept) in the above equation is given by
 * 		    (2)  A = 20 * log10 ( 4*pi*d0/lambda) 
 * 
 * where lambda is the wavelength in meters. 
 * 
 * The path-loss exponent (gamma) is given by
 *                   (3) gamma = (a - (b*ht) + (c/ht)) + x * sigma_gamma ;  10m => ht => 80m
 * where, sigma_gamma is the std deviation of gamma.and x is a zero-mean gaussian variable N[0,1]
 * 
 * The parameters a,b and c depend on the terrain category and are defined below.
 *
 * Terrain Type A (Hilly/Moderate-to-heavy tree density): 		a = 4.6; b = 0.0075; c = 12.6; sigma_gamma = 0.57; mu_sigma = 10.6; sigma_sigma = 2.3;
 * Terrain Type B (Hilly/Flat Moderate-to-heavy tree density): 	a = 4.0; b = 0.0065; c = 17.1; sigma_gamma = 0.75; mu_sigma = 9.6;  sigma_sigma = 3.0;
 * Terrain Type C (Flat/Light tree density): 					a = 3.6; b = 0.005;  c = 20.0; sigma_gamma = 0.59; mu_sigma = 8.2;  sigma_sigma = 1.6;
 *
 * And "s" is the shadow fading component, given by:
 * 						(4) s =  y*sigma
 * where y is a zero-mean Gaussian variable of unit standard deviation; and sigma, the standard deviation of s;
 * 
 * Thus,
 * 					sigma = mu_sigma + z*sigma_sigma;
 * 
 * Corrective action taken for other frequency and receiver antenna heights:
 * PLsui = PLsui + PLdeltaf + PLdeltah + s;
 *
 * where PLdeltaf = 6 * log10 (frequency/2000);
 * for Category A and CategoryB, PLdeltah = -10.8 * log10 (Hr/2.0);
 * for Category C, PLdeltah = -20 * log10 (Hr/2.0);
 */

#include "ns3/propagation-loss-model.h"
#include "ns3/log.h"
#include "ns3/mobility-model.h"
#include "ns3/double.h"
#include "ns3/enum.h"
#include "ns3/pointer.h"
#include <cmath>
#include "sui-loss-model.h"
#include "ns3/random-variable-stream.h"
#define _USE_MATH_DEFINES

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("SUIPathLossModel");
NS_OBJECT_ENSURE_REGISTERED (SUIPathLossModel);

TypeId
SUIPathLossModel::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::SUIPathLossModel")

    .SetParent<PropagationLossModel> ()

    .AddConstructor<SUIPathLossModel> ()

    .AddAttribute ("MinDistance",
    				"The distance under which the propagation model refuses to give results (m). Default = 100m",
                    DoubleValue (100),
                    MakeDoubleAccessor (&SUIPathLossModel::SetMinDistance, &SUIPathLossModel::GetMinDistance),
                    MakeDoubleChecker<double> ())

   	.AddAttribute ("Frequency",
                   "The Frequency  (The frequency range is defined as 2000 MHz).",
                   DoubleValue (2000),
                   MakeDoubleAccessor (&SUIPathLossModel::m_frequency),
                   MakeDoubleChecker<double> ())

	.AddAttribute ("TxAntennaHeight",
				  "Height of the Transmitter Antenna (default is 45m).",
				  DoubleValue (45),
				  MakeDoubleAccessor (&SUIPathLossModel::m_txheight),
				  MakeDoubleChecker<double> ())

	.AddAttribute ("RxAntennaHeight",
				  "Height of the Reciever Antenna between 2m and 10m (default is 2m).",
				   DoubleValue (2),
				   MakeDoubleAccessor (&SUIPathLossModel::m_rxheight),
				   MakeDoubleChecker<double> ())
				   
	.AddAttribute ("Environment",
				  "Type of Terrain Category (default is CategoryA).",
				  EnumValue (CategoryA),
				  MakeEnumAccessor (&SUIPathLossModel::m_environment),
				  MakeEnumChecker (CategoryA, "CategoryA",
                                   CategoryB, "CategoryB",
                                   CategoryC, "CategoryC"))
	
	.AddAttribute ("EnableShadowing",
				  "Enable/Disable Shadowing (s), use 1/0 to enable/disable (default is 1).",
				   DoubleValue (1),
				   MakeDoubleAccessor (&SUIPathLossModel::m_shadowing),
				   MakeDoubleChecker<double> ());

  return tid;
}

SUIPathLossModel::SUIPathLossModel ()
{
}

void
SUIPathLossModel::SetMinDistance (double minDistance)
{
  m_minDistance = minDistance;  //Distance in meters.
}

double
SUIPathLossModel::GetMinDistance (void) const
{
  return m_minDistance;
}

void
SUIPathLossModel::SetFrequency (double frequency)
{
  m_frequency = frequency; // Frequency in GHz.
}

double
SUIPathLossModel::GetFrequency (void) const
{
  return m_frequency;
}

void
SUIPathLossModel::SetTxAntennaHeight (double Hb)
{
  m_txheight = Hb;
}

double
SUIPathLossModel::GetTxAntennaHeight (void)
{
  return m_txheight;
}

void
SUIPathLossModel::SetRxAntennaHeight (double Hm)
{
  m_rxheight = Hm;
}

double
SUIPathLossModel::GetRxAntennaHeight (void)
{
  return m_rxheight;
}

void
SUIPathLossModel::SetShadowing (double sh)
{
  m_shadowing = sh;
}

double
SUIPathLossModel::GetShadowing (void)
{
  return m_shadowing;
}

void
SUIPathLossModel::SetEnvironment (Environment env)
{
  m_environment = env;
}
SUIPathLossModel::Environment
SUIPathLossModel::GetEnvironment (void) const
{
  return m_environment;
}

double
SUIPathLossModel::GetLoss (Ptr<MobilityModel> x, Ptr<MobilityModel> y) const
{
	double mean = 0.0;
	double variance = 1.0;
	
	// Use NS_GLOBAL_VALUE="RngRun=20" from terminal to change the seed for RNG from default 1 to 20. 
	// Ex: $ NS_GLOBAL_VALUE="RngRun=20"  ./waf --run scratch/file-name
	Ptr<NormalRandomVariable> randx = CreateObject<NormalRandomVariable> ();
	Ptr<NormalRandomVariable> randy = CreateObject<NormalRandomVariable> ();
	Ptr<NormalRandomVariable> randz = CreateObject<NormalRandomVariable> ();

	randx->SetAttribute ("Mean", DoubleValue (mean));
	randx->SetAttribute ("Variance", DoubleValue (variance));
	
	randy->SetAttribute ("Mean", DoubleValue (mean));
	randy->SetAttribute ("Variance", DoubleValue (variance));
	
	randz->SetAttribute ("Mean", DoubleValue (mean));
	randz->SetAttribute ("Variance", DoubleValue (variance));
	
	double m_x = randx->GetValue ();
	double m_y = randy->GetValue ();
	double m_z = randz->GetValue ();

	double distance = x->GetDistanceFrom (y);
	double distance_m = distance; //  for distance in m
	if (distance_m < m_minDistance)
    {
      return 0.0;
    }
  
	double d0 = 100; // d0 is defined as 100m.
	
	double m_wavelength = 3e8 / (m_frequency * 1e6);
	double A = 20 * log10 ( 4*M_PI*d0/(m_wavelength));

	NS_LOG_DEBUG ("A  =" << A  << ", Wavelength = " << m_wavelength << ", pi = " << M_PI);

	double a;
	double b;
	double c;
	double sigma_gamma;
	double mu_sigma;
	double sigma_sigma;
	
	if (m_environment == CategoryA ) 
	{
	a = 4.6; b = 0.0075; c = 12.6; sigma_gamma = 0.57; mu_sigma = 10.6; sigma_sigma = 2.3;
  	} 
	else if (m_environment == CategoryB ) 
	{
	a = 4.0; b = 0.0065; c = 17.1; sigma_gamma = 0.75; mu_sigma = 9.6; sigma_sigma = 3.0;
	} 
	else 
	{
	a = 3.6; b = 0.005;  c = 20.0; sigma_gamma = 0.59; mu_sigma = 8.2; sigma_sigma = 1.6;
  	}

		// Enable/Disable Shadowing
		if (m_shadowing == 0)
		{
			m_x = 0;
			m_y = 0;
		}
		
		double m_gamma = a - (b*m_txheight) + (c/m_txheight) + (m_x * sigma_gamma);
		double s = m_y * (mu_sigma + (m_z * sigma_sigma));
		
		NS_LOG_DEBUG ("m_gamma =" << m_gamma << ",   a " << a << ",   b = " << b << ",   c = " << c );
		
		double PLsui = A + (10*m_gamma*log10(distance_m/d0)) + s;

		NS_LOG_DEBUG (" PL of SUI Model =" << PLsui << ",   distance = " << distance_m << ",   H m = " << m_rxheight << ",   H b = " << m_txheight << ",   Frequency = " << m_frequency);

	double PLdeltah;

	double PLdeltaf = 6 * log10 (m_frequency/2000);
	
	if ( (m_environment == CategoryA) || (m_environment ==  CategoryB) )
	{	
		PLdeltah = -10.8 * (log10 (m_rxheight/2.0));
	}
	else	
	{ 
		PLdeltah = -20 * (log10 (m_rxheight/2.0));
	}

		NS_LOG_DEBUG ("PL deltah =" << PLdeltah );
		NS_LOG_DEBUG ("PL deltaf =" << PLdeltaf );


	double loss_in_db = PLsui + PLdeltaf + PLdeltah ;

  	NS_LOG_DEBUG (" Path Loss = " << loss_in_db );

  	return (0 - loss_in_db);

}

double
SUIPathLossModel::DoCalcRxPower (double txPowerDbm, Ptr<MobilityModel> a, Ptr<MobilityModel> b) const
{
  return txPowerDbm + GetLoss (a, b);
}

int64_t
SUIPathLossModel::DoAssignStreams (int64_t stream)
{
  return 0;
}

}
