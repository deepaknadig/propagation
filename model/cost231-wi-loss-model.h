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

#ifndef COST231_WI_LOSS_MODEL_H
#define COST231_WI_LOSS_MODEL_H

#include "ns3/nstime.h"
#include "ns3/propagation-loss-model.h"

namespace ns3 {

/**
 * \ingroup propagation
 *
 *  \brief The COST-231 Walfisch Ikegami Model is a combination of Walfisch and Ikegami models further developed by COST231 Project.
 * 
 * It is now called Empirical COST-Walfisch-Ikegami Model. The model considers
 * only the buildings in the vertical plane between the transmitter and the 
 * receiver. The accuracy of this empirical model is quite high because in
 * urban environments especially the propagation over the rooftops 
 * (multiple diffractions) is the most dominant part. Only wave guiding 
 * effects due to multiple reflections are not considered.  
 * 
 * The main parameters of the model are:
 * Frequency f (800...2000 MHz)
 * Height of the transmitter hTX (4...50 m)
 * Height of the receiver hRX (1...3 m)
 * Distance d between transmitter and receiver (20...5000 m)
 *
 */

class Cost231WILossModel : public PropagationLossModel
{

public:
  static TypeId GetTypeId (void);
  Cost231WILossModel ();
  enum Environment
  {
    Suburban, Urban
  };

  /**
   * \param a the mobility model of the source
   * \param b the mobility model of the destination
   * \returns the propagation loss (in dBm)
   */
  void SetFrequency (double frequency);
  double GetFrequency (void) const;

  void SetWidth (double width);
  double GetWidth (void);

  void SetRoofHeight (double hroof);
  double GetRoofHeight (void);

  void SetMobileHeight (double hmobile);
  double GetMobileHeight (void);

  void SetOrientationAngle (double oriangle);
  double GetOrientationAngle (void);

  void SetBaseHeight (double hbase);
  double GetBaseHeight (void);

  void SetEnvironment (Environment env);
  Environment GetEnvironment (void) const;

  double GetLoss (Ptr<MobilityModel> a, Ptr<MobilityModel> b) const;

  void SetMinDistance (double minDistance);
  double GetMinDistance (void) const;

private:
  virtual double DoCalcRxPower (double txPowerDbm, Ptr<MobilityModel> a, Ptr<MobilityModel> b) const;
  virtual int64_t DoAssignStreams (int64_t stream);

  double m_hroof; // in meter
  double m_hmobile; // in meter
  double m_hbase; // in meter
  double m_oriangle; // in degrees
  Environment m_environment;
  double m_minDistance; // in meter
  double m_frequency; // frequency in MHz
  double m_width; // width of the road in meters
  //double m_b; // building separation in meters
};

}

#endif /* COST231WIMODEL_H */
