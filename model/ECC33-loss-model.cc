//ECC-33 propagation path loss model Logic
//=========================================
//	 PLecc33 = Afs + Abm - Gb - Gr -----------------------------------------------------------------------eq 1

//	Afs - Free Space Attenuation in dB
//	Abm - Basic Median path loss in dB
//	Gb - Transmitter Antenna height gain factor
//	Gr - Receiver Antenna height gain factor 

//	Afs = 92.4 + ( 20 * log10(d) ) + (20 * log10(f) ) --------------------------------------------- eq 2

//	Abm = 20.41 + ( 9.83 * log10(d) ) + ( 7.894 * log10(f) ) + ( 9.56 * (( log10(f))*( log10(f))) -- eq 3

//	Gb = log10(Hb/200)*(13.958 + (5.8*((log10(d))*(log10(d))) -------------------------------------- eq 4

//	For medium cities,

//	Gr = (42.57 + (13.7*log10(f))*((log10(Hr)) - 0.585 ) ------------------------------------------- eq 5

//	For large cities,

//	Gr = ( 0.759 * Hr ) - 1.892 -------------------------------------------------------------------- eq 6

//	Hb - Tx Antenna Height =  m; 
//	Hr - Rx Antenna Height; 
//	d - distance between Tx and Rx in Km; 
//	f - Frequeny in Ghz;

//===============================================================================================================

#include "ns3/propagation-loss-model.h"
#include "ns3/log.h"
#include "ns3/mobility-model.h"
#include "ns3/double.h"
#include "ns3/pointer.h"
#include <cmath>
#include "ECC33-path-loss-model.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("ECC33PathLossModel");
NS_OBJECT_ENSURE_REGISTERED (ECC33PathLossModel);

TypeId
ECC33PathLossModel::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::ECC33PathLossModel")

    .SetParent<PropagationLossModel> ()

    .AddConstructor<ECC33PathLossModel> ()

    .AddAttribute ("MinDistance",
    				"The distance under which the propagation model refuses to give results (m). Default = 20m",
                    DoubleValue (20),
                    MakeDoubleAccessor (&ECC33PathLossModel::SetMinDistance, &ECC33PathLossModel::GetMinDistance),
                    MakeDoubleChecker<double> ())

    .AddAttribute ("Frequency",
                   "The Frequency  (The frequency range is defined as 2 GHz).",
                   DoubleValue (2),
                   MakeDoubleAccessor (&ECC33PathLossModel::m_frequency),
                   MakeDoubleChecker<double> ())

     .AddAttribute ("TxAntennaHeight",
				  "Height of the Transmitter Antenna (default is 30m).",
				  DoubleValue (50),
				  MakeDoubleAccessor (&ECC33PathLossModel::m_txheight),
				  MakeDoubleChecker<double> ())

	  .AddAttribute ("RxAntennaHeight",
				  "Height of the Reciever Antenna (default is 6m).",
				   DoubleValue (2),
				   MakeDoubleAccessor (&ECC33PathLossModel::m_rxheight),
				   MakeDoubleChecker<double> ());

  return tid;
}

ECC33PathLossModel::ECC33PathLossModel ()
{
}

void
ECC33PathLossModel::SetMinDistance (double minDistance)
{
  m_minDistance = minDistance/1000; //Distance in KM.
}

double
ECC33PathLossModel::GetMinDistance (void) const
{
  return m_minDistance;
}

void
ECC33PathLossModel::SetFrequency (double frequency)
{
  m_frequency = frequency; // Frequency in MHz.
}

double
ECC33PathLossModel::GetFrequency (void) const
{
  return m_frequency;
}



void
ECC33PathLossModel::SetTxAntennaHeight (double Hb)
{
  m_txheight = Hb;
}

double
ECC33PathLossModel::GetTxAntennaHeight (void)
{
  return m_txheight;
}

void
ECC33PathLossModel::SetRxAntennaHeight (double Hr)
{
  m_rxheight = Hr;
}

double
ECC33PathLossModel::GetRxAntennaHeight (void)
{
  return m_rxheight;
}

void
ECC33PathLossModel::SetEnvironment (Environment env)
{
  m_environment = env;
}
ECC33PathLossModel::Environment
ECC33PathLossModel::GetEnvironment (void) const
{
  return m_environment;
}


double
ECC33PathLossModel::GetLoss (Ptr<MobilityModel> a, Ptr<MobilityModel> b) const
{

  double distance = a->GetDistanceFrom (b);
  double distance_km = distance / 1000;
  if (distance_km <= m_minDistance)
    {
      return 0.0;
    }
  
	double Afs = 92.4 + ( 20 * log10(distance_km) ) + ( 20 * log10(m_frequency) ) ;

NS_LOG_DEBUG ("A fs =" << Afs );

	double Abm = 20.41 + (9.83 * log10(distance_km)) + (7.894 * log10(m_frequency)) + (9.56 * ((log10(m_frequency))*(log10(m_frequency))));

NS_LOG_DEBUG ("A bm =" << Abm );

	double Gb = (log10(m_txheight/200)) * (13.958 + (5.8 * ( (log10(distance_km))*(log10(distance_km)) ) ) );


NS_LOG_DEBUG ("G b =" << Gb );
 
//	Afs - Free Space Attenuation in dB
//	Abm - Basic Median path loss in dB
//	Gb - Transmitter Antenna height gain factor
//	Gr - Receiver Antenna height gain factor 
//	Hb - Tx Antenna Height =  30m; 
//	Hr - Rx Antenna Height =  6m ; 
//	d - distance between Tx and Rx in meters = 5000m; 
//	f - Frequeny in MHz = 2000 MHz;

double Gr;

//	For medium cities,

if (m_environment == Suburban) {	

Gr = (42.57 + (13.7*log10(m_frequency))) * ((log10(m_rxheight)) - 0.585 ); 

NS_LOG_DEBUG ("G r =" << Gr );

}

//	For large cities,

else	{ 

Gr = ( 0.759 * m_rxheight ) - 1.892 ;

NS_LOG_DEBUG ("G r =" << Gr );
}

// ECC33 Path Loss model equation

double loss_in_db = Afs + Abm - Gb - Gr;

  NS_LOG_DEBUG ("dist =" << distance << ", Path Loss = " << loss_in_db << ", G r = " << Gr << ", G b = " << Gb << ", freq = " << m_frequency << ", Tx antenna height = " << m_txheight << ", Rx antenna height = " << m_rxheight << ", A fs = " << Afs << ", A bm = " << Abm);

  return (0 - loss_in_db);

}

double
ECC33PathLossModel::DoCalcRxPower (double txPowerDbm, Ptr<MobilityModel> a, Ptr<MobilityModel> b) const
{
  return txPowerDbm + GetLoss (a, b);
}

int64_t
ECC33PathLossModel::DoAssignStreams (int64_t stream)
{
  return 0;
}

}





