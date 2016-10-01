/* @(#)root/histpainter:$Id$ */

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class TPainter3dAlgorithms;
#pragma link C++ class TGraph2DPainter;
#pragma link C++ class TGraphPainter;
#pragma link C++ class THistPainter;
#pragma link C++ class TPaletteAxis+;

// needed since new class definition of TGraph2DPainter
#pragma extra_include "TGraph2D.h";
#pragma extra_include "TGraphDelaunay.h";
#pragma extra_include "TGraphDelaunay2D.h";


#endif