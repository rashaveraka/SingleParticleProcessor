// @(#)root/tmva $Id$
// Author: Helge Voss

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TActivationReLU                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Rectified linear unit function  for an ANN.                               *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Helge Voss                                                                *
 *                                                                                *
 * Copyright (c) 2005:                                                            *
 *      CERN, Switzerland                                                         *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/
  

#include <iostream>

#include "TFormula.h"
#include "TString.h"
#include "TMath.h"

#ifndef ROOT_TMVA_TActivationReLU
#include "TMVA/TActivationReLU.h"
#endif

ClassImp(TMVA::TActivationReLU)

////////////////////////////////////////////////////////////////////////////////
/// constructor for ReLU

TMVA::TActivationReLU::TActivationReLU()
{
   // sorry, I really don't know what I would possibly want to do here ;)

}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TMVA::TActivationReLU::~TActivationReLU()
{
}

////////////////////////////////////////////////////////////////////////////////
/// get expressions for the tanh and its derivative

TString TMVA::TActivationReLU::GetExpression()
{
   TString expr = "max(0,x)";

   return expr;
}

////////////////////////////////////////////////////////////////////////////////
/// writes the sigmoid activation function source code

void TMVA::TActivationReLU::MakeFunction( std::ostream& fout, const TString& fncName ) 
{
   fout << "double " << fncName << "(double x) const {" << std::endl;
   fout << "   // rectified linear unit" << std::endl;
   fout << "   return x>0 ? x : 0; " << std::endl;
   fout << "}" << std::endl;
}
