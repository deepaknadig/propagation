//SUI propagation path loss model Logic
//=========================================
//The SUI model covers three terrain categories. 
//Category A represents the maximum path-loss category which is a hilly terrain, 
//Category B represents an intermediate path-loss category, and 
//Category C represents the minimum path-loss category with mostly flat terrains.  
//
//The empical formulas for this model are as below: //
//
//The median path-loss for the SUI model can be generally written as
//
//                  (1)  PLsui = A + 10*u*log10(d/d0) ;
//
//for d>do, where do=100m. The term A in the above equation is given by 
//
//		    (2)  A = 20 * log10 ( 4*22*d0/lambda) ;
//
//where lambda is the wavelength in m. The path-loss exponent is given by
//
//                  (3) u = a - (b*ht) + (c/ht) ;
//
//in which the parameters a,b and c depend on the terrain category and are defined below.

//Terrain Type A: a = 4.6; b = 0.0075; c = 12.6;
//Terrain Type B: a = 4.0; b = 0.0065; c = 17.1;
//Terrain Type C: a = 3.6; b = 0.005;  c = 20.0;
//
//Corrective action taken for other frequency and reciever antenna heights:
//
// PLsui = PLsui + PLdeltaf + PLdeltah ;
//
// where PLdeltaf = 6 * log10 (frequency/2000);
//
// for CategoryA and CategoryB, PLdeltah = -10.8 * log10 (Hr/2.0);
//
// for CategoryC, PLdeltah = -20 * log10 (Hr/2.0);
//
//==========================================================================================


#include "ns3/propagation-loss-model.h"
#include "ns3/log.h"
#include "ns3/mobility-model.h"
#include "ns3/double.h"
#include "ns3/pointer.h"
#include <cmath>
#include "SUI-path-loss-model.h"

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
    				"The distance under which the propagation model refuses to give results (m). Default = 20m",
                    DoubleValue (200),
                    MakeDoubleAccessor (&SUIPathLossModel::SetMinDistance, &SUIPathLossModel::GetMinDistance),
                    MakeDoubleChecker<double> ())

    	.AddAttribute ("Frequency",
                   "The Frequency  (The frequency range is defined as 2000 MHz).",
                   DoubleValue (2000),
                   MakeDoubleAccessor (&SUIPathLossModel::m_frequency),
                   MakeDoubleChecker<double> ())


  	.AddAttribute ("Wavelength",
                   "The Wavelength  (The Wavelength range is defined in m ).",
                   DoubleValue (0.15),
                   MakeDoubleAccessor (&SUIPathLossModel::m_wavelength),
                   MakeDoubleChecker<double> ())
  
	.AddAttribute ("TxAntennaHeight",
				  "Height of the Transmitter Antenna (default is 50m).",
				  DoubleValue (50),
				  MakeDoubleAccessor (&SUIPathLossModel::m_txheight),
				  MakeDoubleChecker<double> ())

	.AddAttribute ("RxAntennaHeight",
				  "Height of the Reciever Antenna (default is 2m).",
				   DoubleValue (6),
				   MakeDoubleAccessor (&SUIPathLossModel::m_rxheight),
				   MakeDoubleChecker<double> ());

  return tid;
}

SUIPathLossModel::SUIPathLossModel ()
{
}

void
SUIPathLossModel::SetMinDistance (double minDistance)
{
  m_minDistance = minDistance/1000;  // /1000Distance in KM.
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
SUIPathLossModel::SetWavelength (double wavelength)
{
  m_wavelength = wavelength; // Wavelength is in m.
}

double
SUIPathLossModel::GetWavelength (void) const
{
  return m_wavelength;
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

  double distance = x->GetDistanceFrom (y);
  double distance_km = distance; //  for distance in km
  if (distance_km <= m_minDistance)
    {
      return 0.0;
    }
  
	
//parameters a,b and c depend on the terrain category and are defined below.
//Terrain Type A: a = 4.6; b = 0.0075; c = 12.6;
//Terrain Type B: a = 4.0; b = 0.0065; c = 17.1;
//Terrain Type C: a = 3.6; b = 0.005;  c = 20.0;

	double d0 = 100; // d0 is defined as 100m ie 0.1 in KM.

	double A = 20 * log10 ( 4*22*d0/(7*m_wavelength)) ;

	NS_LOG_DEBUG ("A  =" << A );

	double a;
	double b;
	double c;

	if (m_environment == CategoryA ) 
	{
	a = 4.6; b = 0.0075; c = 12.6;
  	} 
	else if (m_environment == CategoryB ) 
	{
	a = 4.0; b = 0.0065; c = 17.1;
	} 
	else 
	{
	a = 3.6; b = 0.005;  c = 20.0;
  	}

  
		double u = a - (b*m_txheight) + (c/m_txheight) ;

		NS_LOG_DEBUG ("u =" << u << ",   a " << a << ",   b = " << b << ",   c = " << c );

// SUI Path Loss Model Equation		

		double PLsui = A + (10*u*log10(distance_km/d0)) ;

		NS_LOG_DEBUG (" PL of SUI Model =" << PLsui << ",   distance = " << distance_km << ",   H m = " << m_rxheight << ",   H b = " << m_txheight << ",   Frequency = " << m_frequency);
 

//	Hb - Tx Antenna Height =  50m; 
//	Hm - Rx Antenna Height =  6m ; 
//	d0 = 100 m; 
//	f - Frequeny in GHz = 2 GHz;


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


// SUI Path Loss model equation for corrective action

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





