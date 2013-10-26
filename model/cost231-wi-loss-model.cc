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

#include "ns3/propagation-loss-model.h"
#include "ns3/log.h"
#include "ns3/mobility-model.h"
#include "ns3/double.h"
#include "ns3/pointer.h"
#include <cmath>
#include "cost231-wi-loss-model.h"

namespace ns3 {

NS_LOG_COMPONENT_DEFINE ("Cost231WILossModel");
NS_OBJECT_ENSURE_REGISTERED (Cost231WILossModel);

TypeId
Cost231WILossModel::GetTypeId (void)
{
  static TypeId tid = TypeId ("ns3::Cost231WILossModel")

    .SetParent<PropagationLossModel> ()

    .AddConstructor<Cost231WILossModel> ()

    .AddAttribute ("MinDistance",
    				"The distance under which the propagation model refuses to give results (m). Default = 20m",
                    DoubleValue (20),
                    MakeDoubleAccessor (&Cost231WILossModel::SetMinDistance, &Cost231WILossModel::GetMinDistance),
                    MakeDoubleChecker<double> ())

    .AddAttribute ("Width",
                   "The width of the road can range between 10 to 25m  (default is 10m).",
                   DoubleValue (10),
                   MakeDoubleAccessor (&Cost231WILossModel::m_width),
                   MakeDoubleChecker<double> ())

    .AddAttribute ("Frequency",
                   "The Frequency in MHz (default is 2000 MHz).",
                   DoubleValue (2000),
                   MakeDoubleAccessor (&Cost231WILossModel::m_frequency),
                   MakeDoubleChecker<double> ())

    .AddAttribute ("OrientationAngle",
                   "Orientation of the street w.r.t LoS (default is 90 degrees).",
                   DoubleValue (90.0),
                   MakeDoubleAccessor (&Cost231WILossModel::m_oriangle),
                   MakeDoubleChecker<double> ())

    .AddAttribute ("RoofHeight",
                   "Height of the building roof (default is 6m).",
                   DoubleValue (6),
                   MakeDoubleAccessor (&Cost231WILossModel::m_hroof),
                   MakeDoubleChecker<double> ())

    .AddAttribute ("MobileHeight",
				  "Height of the MS (default is 3m).",
				  DoubleValue (3),
				  MakeDoubleAccessor (&Cost231WILossModel::m_hmobile),
				  MakeDoubleChecker<double> ())

	  .AddAttribute ("BaseHeight",
						 "Height of the BS (default is 30m).",
						 DoubleValue (30),
						 MakeDoubleAccessor (&Cost231WILossModel::m_hbase),
						 MakeDoubleChecker<double> ());

  return tid;
}

Cost231WILossModel::Cost231WILossModel ()
{
}

void
Cost231WILossModel::SetMinDistance (double minDistance)
{
  m_minDistance = minDistance/1000; //Distance in KM.
}

double
Cost231WILossModel::GetMinDistance (void) const
{
  return m_minDistance;
}

void
Cost231WILossModel::SetFrequency (double frequency)
{
  m_frequency = frequency; // Frequency in MHz.
}

double
Cost231WILossModel::GetFrequency (void) const
{
  return m_frequency;
}

void
Cost231WILossModel::SetWidth (double width)
{
  m_width = width;
}

double
Cost231WILossModel::GetWidth (void)
{
  return m_width;
}

void
Cost231WILossModel::SetRoofHeight (double hroof)
{
  m_hroof = hroof;
}

double
Cost231WILossModel::GetRoofHeight (void)
{
  return m_hroof;
}

void
Cost231WILossModel::SetMobileHeight (double hmobile)
{
  m_hmobile = hmobile;
}

double
Cost231WILossModel::GetMobileHeight (void)
{
  return m_hmobile;
}

void
Cost231WILossModel::SetOrientationAngle (double oriangle)
{
  m_oriangle = oriangle; //Orientation Angle Phi in Degrees.
}

double
Cost231WILossModel::GetOrientationAngle (void)
{
  return m_oriangle;
}

void
Cost231WILossModel::SetBaseHeight (double hbase)
{
  m_hbase = hbase;
}

double
Cost231WILossModel::GetBaseHeight (void)
{
  return m_hbase;
}

void
Cost231WILossModel::SetEnvironment (Environment env)
{
  m_environment = env;
}
Cost231WILossModel::Environment
Cost231WILossModel::GetEnvironment (void) const
{
  return m_environment;
}


double
Cost231WILossModel::GetLoss (Ptr<MobilityModel> a, Ptr<MobilityModel> b) const
{

  double distance = a->GetDistanceFrom (b);
  double distance_km = distance / 1000;
  if (distance_km <= m_minDistance)
    {
      return 0.0;
    }
  // Calculation of L0 (free space loss)
  double L0 = 32.4 + 20 * log10(distance_km) + 20 * log10 (m_frequency);
  
  // Calculation of Lrts (roof top to street loss)

  double Lori;

  if ((0 <= m_oriangle) && (m_oriangle < 35)) {
	Lori = -10 + 0.354 * m_oriangle;
  } else if ((35 <= m_oriangle) && (m_oriangle < 55)) {
	Lori = 2.5 + 0.075 * (m_oriangle - 35);
  } else {
	Lori = 4.0 - 0.114 * (m_oriangle - 35);
  }
  
  double delta_hmobile = m_hroof - m_hmobile;
  
  double Lrts = -16.9 - (10 * log10(m_width)) + (10 * log10(m_frequency)) + (20 * log10(delta_hmobile)) + Lori;

  // Calculation of Lmsd (multiple screen diffraction loss)

  double delta_hbase = m_hbase - m_hroof;

  double Lbsh;

  if (m_hbase > m_hroof) {
	  Lbsh = -18 * (log10(1 + delta_hbase));
  } else {
	  Lbsh = 0;
  }

  double Ka;

  if (m_hbase > m_hroof) {
  	Ka = 54;
  } else if ((distance_km >= 0.5) && (m_hbase <= m_hroof)) {
  	Ka = 54 - 0.8 * delta_hbase;
    } else {
  	Ka = 54 - (1.6 * delta_hbase * distance_km);
  }

  double Kd;
  double Kf;
  double loss_in_db;

  if (m_hbase > m_hroof) {
	  Kd = 18;
  } else {
	Kd = 18 - (15 * (delta_hbase/m_hroof));
  }

  if (m_environment == Suburban) {
	  Kf = -4 + 0.7 * ((m_frequency/925) - 1);
  } else {
	  Kf = -4 + 1.5 * ((m_frequency/925) - 1);
  }
  double m_b = m_width * 2;
  double Lmsd = Lbsh + Ka + (Kd * log10(distance_km)) + (Kf * log10(m_frequency)) - (9 * log10(m_b));

  if ((Lrts + Lmsd) > 0) {
	loss_in_db = L0 + Lrts + Lmsd;
  } else {
	loss_in_db = L0;
  }


  NS_LOG_DEBUG ("dist =" << distance_km << ",   Lbsh = " << Lbsh << ", Path Loss = " << loss_in_db  << ",    m_b = " << m_b << ",   Kd = " << Kd  << ",   Lmsd = " << Lmsd << ",   Kf = " << Kf);

  return (0 - loss_in_db);

}

double
Cost231WILossModel::DoCalcRxPower (double txPowerDbm, Ptr<MobilityModel> a, Ptr<MobilityModel> b) const
{
  return txPowerDbm + GetLoss (a, b);
}

int64_t
Cost231WILossModel::DoAssignStreams (int64_t stream)
{
  return 0;
}

}
