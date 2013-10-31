/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2007,2008, 2009 INRIA, UDcast
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
 * Author: Mohamed Amine Ismail <amine.ismail@sophia.inria.fr>
 *                              <amine.ismail@udcast.com>
 */

#ifndef SUI_LOSS_MODEL_H
#define SUI_LOSS_MODEL_H

#include "ns3/nstime.h"
#include "ns3/propagation-loss-model.h"

namespace ns3 {

/**
 *  Frequency: 2000 MHz
 *  Tx Antenna Height: 50m
 *  Rx Antenna Height: 2m
 *  Tx-Rx Distance: 1 km
 *
 */

class SUIPathLossModel : public PropagationLossModel
{

public:
  static TypeId GetTypeId (void);
  SUIPathLossModel ();
  enum Environment
  {
    CategoryA, CategoryB, CategoryC
  };

  /**
   * \param a the mobility model of the source
   * \param b the mobility model of the destination
   * \returns the propagation loss (in dBm)
   */
  void SetFrequency (double frequency);
  double GetFrequency (void) const;

  void SetTxAntennaHeight (double Hb);
  double GetTxAntennaHeight (void);

  void SetRxAntennaHeight (double Hm);
  double GetRxAntennaHeight (void);

  void SetEnvironment (Environment env);
  Environment GetEnvironment (void) const;

  double GetLoss (Ptr<MobilityModel> a, Ptr<MobilityModel> b) const;

  void SetMinDistance (double minDistance);
  double GetMinDistance (void) const;

private:
  virtual double DoCalcRxPower (double txPowerDbm, Ptr<MobilityModel> a, Ptr<MobilityModel> b) const;
  virtual int64_t DoAssignStreams (int64_t stream);

  double m_txheight; // in meter
  double m_rxheight; // in meter
  Environment m_environment;
  double m_minDistance; // in meter
  double m_frequency; // frequency in GHz
  double m_wavelength; // wavelength in lambda
  
};

}

#endif /* SUIPATHLOSSMODEL_H */
