/// \file ROOT/THistDrawable.h
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-07-09
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_THistDrawable
#define ROOT7_THistDrawable

#include "ROOT/TCoopPtr.h"
#include "ROOT/TDrawable.h"
#include "ROOT/THistDrawOptions.h"
#include "ROOT/TLogger.h"

#include "TSystem.h"

#include <memory>

namespace ROOT {

template<int DIMENSIONS, class PRECISION> class THist;

namespace Internal {

template <int DIMENSION>
class THistPainterBase {
  static THistPainterBase<DIMENSION>* fgPainter;

protected:
  THistPainterBase() { fgPainter = this; }
  virtual ~THistPainterBase();

public:
  static THistPainterBase<DIMENSION>* GetPainter() {
    if (!fgPainter)
      gSystem->Load("libHistPainter");
    return fgPainter;
  }

  /// Paint a THist. All we need is access to its GetBinContent()
  virtual void Paint(TDrawable& obj, THistDrawOptions<DIMENSION> opts) = 0;
};

extern template class THistPainterBase<1>;
extern template class THistPainterBase<2>;
extern template class THistPainterBase<3>;

template <int DIMENSION, class PRECISION>
class THistDrawable final: public TDrawable {
private:
  TCoopPtr<THist<DIMENSION, PRECISION>> fHist;
  THistDrawOptions<DIMENSION> fOpts;

public:
  THistDrawable(TCoopPtr<THist<DIMENSION, PRECISION>> hist,
                THistDrawOptions<DIMENSION> opts): fHist(hist), fOpts(opts) {}

  ~THistDrawable() = default;

  /// Paint the histogram
  void Paint() final {
    THistPainterBase<DIMENSION>::GetPainter()->Paint(*this, fOpts);
  }
};

} // namespace Internal
} // namespace ROOT

#endif
