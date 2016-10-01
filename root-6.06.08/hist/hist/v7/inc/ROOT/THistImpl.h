/// \file ROOT/THistImpl.h
/// \ingroup Hist ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-03-23
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_THistImpl
#define ROOT7_THistImpl

#include <cctype>
#include "ROOT/RArrayView.h"
#include "ROOT/RTupleApply.h"

#include "ROOT/TAxis.h"

namespace ROOT {

namespace Hist {
/// Iterator over n dimensional axes - an array of n axis iterators.
template<int NDIM> using AxisIter_t = std::array<TAxisBase::const_iterator, NDIM>;
/// Range over n dimensional axes - a pair of arrays of n axis iterators.
template<int NDIM> using AxisIterRange_t = std::array<AxisIter_t<NDIM>, 2>;

/// Kinds of under- and overflow handling.
enum class EOverflow {
  kNoOverflow = 0x0, ///< Exclude under- and overflows
  kUnderflow = 0x1, ///< Include underflows
  kOverflow = 0x2, ///< Include overflows
  kUnderOver = 0x3, ///< Include both under- and overflows
};

inline bool operator&(EOverflow a, EOverflow b) {
  return static_cast<int>(a) & static_cast<int>(b);
}
}

namespace Detail {

/**
 \class THistImplPrecisionAgnosticBase
 Base class for THistImplBase that abstracts out the histogram's PRECISION.

 For operations such as painting a histogram, the PRECISION (type of the bin
 content) is not relevant; painting will cast the underlying bin type to double.
 To facilitate this, THistImplBase itself inherits from the
 THistImplPrecisionAgnosticBase interface.
 */
template <int DIMENSIONS>
class THistImplPrecisionAgnosticBase {
public:
  /// Type of the coordinate: a DIMENSIONS-dimensional array of doubles.
  using Coord_t = std::array<double, DIMENSIONS>;

  virtual ~THistImplPrecisionAgnosticBase() {}

  /// Number of dimensions of this histogram.
  constexpr int GetNDim() const { return DIMENSIONS; }
  /// Number of bins of this histogram, including all overflow and underflow
  /// bins. Simply the product of all axes' number of bins.
  virtual int GetNBins() const = 0;

  /// Given the coordinate `x`, determine the index of the bin.
  virtual int GetBinIndex(const Coord_t& x) const = 0;
  /// Given the coordinate `x`, determine the index of the bin, possibly growing
  /// axes for which `x` is out of range.
  virtual int GetBinIndexAndGrow(const Coord_t& x) = 0;

  /// Get the center in all dimensions of the bin with index `binidx`.
  virtual Coord_t GetBinCenter(int binidx) const = 0;
  /// Get the lower edge in all dimensions of the bin with index `binidx`.
  virtual Coord_t GetBinFrom(int binidx) const = 0;
  /// Get the upper edge in all dimensions of the bin with index `binidx`.
  virtual Coord_t GetBinTo(int binidx) const = 0;

  /// The bin's uncertainty. size() of the vector is a multiple of 2:
  /// several kinds of uncertainty, same number of entries for lower and upper.
  virtual std::vector<double> GetBinUncertainties(int binidx) const = 0;

  /// The bin content, cast to double.
  virtual double GetBinContentAsDouble(int binidx) const = 0;

  /// Get a TAxisView on axis with index iAxis.
  ///
  /// \param iAxis - index of the axis, must be 0 <= iAxis < DIMENSION
  virtual TAxisView GetAxis(int iAxis) const = 0;

  /// Get a Hist::AxisIterRange_t for the whole histogram, possibly
  /// restricting the range to non-overflow bins.
  ///
  /// \param withOverUnder - specifies for each dimension whether under and
  /// overflow should be included in the returned range.
  virtual Hist::AxisIterRange_t<DIMENSIONS>
    GetRange(const std::array<Hist::EOverflow, DIMENSIONS>& withOverUnder) const = 0;
};


/**
 \class THistImplBase
 Interface class for THistImpl.

 THistImpl is templated for a specific configuration of axes. To enable access
 through THist, THistImpl inherits from THistImplBase, exposing only dimension
 (`DIMENSION`) and bin type (`PRECISION`).
 */
template<int DIMENSIONS, class PRECISION>
class THistImplBase: public THistImplPrecisionAgnosticBase<DIMENSIONS> {
public:
  /// Type of a coordinate: an array of `DIMENSIONS` doubles.
  using Coord_t = typename THistImplPrecisionAgnosticBase<DIMENSIONS>::Coord_t;
  /// Type of the bin content (and thus weights).
  using Weight_t = PRECISION;
  /// Type of the Fill(x, w) function
  using FillFunc_t = void (THistImplBase::*)(const Coord_t& x, Weight_t w);

  /// Interface function to fill a vector or array of coordinates with
  /// corresponding weights.
  /// \note the size of `xN` and `weightN` must be the same!
  virtual void FillN(const std::array_view<Coord_t> xN,
                     const std::array_view<Weight_t> weightN) = 0;

  /// Interface function to fill a vector or array of coordinates.
  virtual void FillN(const std::array_view<Coord_t> xN) = 0;

  /// Retrieve the pointer to the overridden Fill(x, w) function.
  virtual FillFunc_t GetFillFunc() const = 0;


  /// Get the bin content (sum of weights) for bin index `binidx`.
  virtual PRECISION GetBinContent(int binidx) const = 0;

  /// Get the bin content (sum of weights) for bin index `binidx`, cast to
  /// double.
  double GetBinContentAsDouble(int binidx) const final {
    return (double) GetBinContent(binidx);
  }
};
} // namespace Detail


namespace Internal {
/** \name Histogram traits
    Helper traits for histogram operations.
 */
///\{

/// \name AxisTupleOperations
/// Template operations on axis tuple.
///@{
template <int IDX, class AXISTUPLE> struct TGetBinCount;

template <class AXES> struct TGetBinCount<0, AXES> {
  int operator()(const AXES& axes) const {
    return std::get<0>(axes).GetNBins();
  }
};


template <int I, class AXES>
struct TGetBinCount {
  int operator()(const AXES& axes) const {
    return std::get<I>(axes).GetNBins() * TGetBinCount<I - 1, AXES>()(axes);
  }
};


template <int IDX, class HISTIMPL, class AXES, bool GROW>
struct TGetBinIndex;

// Break recursion
template <class HISTIMPL, class AXES, bool GROW>
struct TGetBinIndex< -1, HISTIMPL, AXES, GROW> {
  int operator()(HISTIMPL*, const AXES&, const typename HISTIMPL::Coord_t&,
                 TAxisBase::EFindStatus& status) const {
    status = TAxisBase::EFindStatus::kValid;
    return 0;
  }
};

template <int I, class HISTIMPL, class AXES, bool GROW>
struct TGetBinIndex {
  int operator()(HISTIMPL* hist, const AXES& axes,
                 const typename HISTIMPL::Coord_t& x, TAxisBase::EFindStatus& status) const {
    int bin = std::get<I>(axes).FindBin(x[I]);
    if (GROW && std::get<I>(axes).CanGrow()
        && (bin < 0 || bin > std::get<I>(axes).GetNBinsNoOver())) {
      hist->GrowAxis(I, x[I]);
      status = TAxisBase::EFindStatus::kCanGrow;

      // Abort bin calculation; we don't care. Let THist::GetBinIndex() retry!
      return bin;
    }
    return bin + TGetBinIndex<I - 1, HISTIMPL, AXES, GROW>()(hist, axes, x, status)
                 * std::get<I>(axes).GetNBins();
  }
};


template<int I, class AXES> struct FillIterRange_t;

// Break recursion.
template<class AXES> struct FillIterRange_t<-1, AXES> {
  void operator()(Hist::AxisIterRange_t<std::tuple_size<AXES>::value>& /*range*/,
                  const AXES& /*axes*/,
                  const std::array<Hist::EOverflow, std::tuple_size<AXES>::value>& /*over*/) const {}
};

/** Fill `range` with begin() and end() of all axes, including under/overflow
  as specified by `over`.
*/
template<int I, class AXES>
struct FillIterRange_t {
  void operator()(Hist::AxisIterRange_t<std::tuple_size<AXES>::value> &range,
                  const AXES &axes,
                  const std::array<Hist::EOverflow, std::tuple_size<AXES>::value> &over) const {
    if (over[I] & Hist::EOverflow::kUnderflow)
      range[0][I] = std::get<I>(axes).begin_with_underflow();
    else
      range[0][I] = std::get<I>(axes).begin();
    if (over[I] & Hist::EOverflow::kOverflow)
      range[1][I] = std::get<I>(axes).end_with_overflow();
    else
      range[1][I] = std::get<I>(axes).end();
    FillIterRange_t<I - 1, AXES>()(range, axes, over);
  }
};



enum class EBinCoord {
  kBinFrom, ///< Get the lower bin edge
  kBinCenter, ///< Get the bin center
  kBinTo ///< Get the bin high edge
};

template<int I, class COORD, class AXES> struct FillBinCoord_t;

// Break recursion.
template<class COORD, class AXES> struct FillBinCoord_t<-1, COORD, AXES> {
  void operator()(COORD& /*coord*/, const AXES& /*axes*/, EBinCoord /*kind*/, int /*binidx*/) const {}
};

/** Fill `coord` with low bin edge or center or high bin edge of all axes.
*/
template<int I, class COORD, class AXES>
struct FillBinCoord_t {
  void operator()(COORD& coord, const AXES& axes, EBinCoord kind, int binidx) const {
    int axisbin = binidx % std::get<I>(axes).GetNBins();
    size_t coordidx = std::tuple_size<AXES>::value - I;
    switch (kind) {
      case EBinCoord::kBinFrom:
        coord[coordidx] = std::get<I>(axes).GetBinFrom(axisbin);
        break;
      case EBinCoord::kBinCenter:
        coord[coordidx] = std::get<I>(axes).GetBinCenter(axisbin);
        break;
      case EBinCoord::kBinTo:
        coord[coordidx] = std::get<I>(axes).GetBinTo(axisbin);
        break;
    }
    FillBinCoord_t<I - 1, COORD, AXES>()(coord, axes, kind,
                                         binidx / std::get<I>(axes).GetNBins());
  }
};



template <class... AXISCONFIG>
static std::array<TAxisView, sizeof...(AXISCONFIG)>
GetAxisView(const AXISCONFIG&...axes) noexcept {
  std::array<TAxisView, sizeof...(AXISCONFIG)> axisViews = {
     {TAxisView(axes)...}
  };
  return axisViews;
}

///\}
} // namespace Internal


template <int DIMENSIONS, class PRECISION> class THist;

namespace Detail {

template <int DIMENSIONS, class PRECISION, class STATISTICS, class... AXISCONFIG>
class THistImpl final: public THistImplBase<DIMENSIONS, PRECISION>,
   STATISTICS {
  static_assert(sizeof...(AXISCONFIG) == DIMENSIONS,
                "Number of axes must equal histogram dimension");
  friend class THist<DIMENSIONS, PRECISION>;

public:
  using ImplBase_t = THistImplBase<DIMENSIONS, PRECISION>;
  using Coord_t = typename ImplBase_t::Coord_t;
  using Weight_t = typename ImplBase_t::Weight_t;
  using typename ImplBase_t::FillFunc_t;
  template <int NDIM = DIMENSIONS> using AxisIterRange_t
    = typename Hist::AxisIterRange_t<NDIM>;

private:
  /// Get the number of bins in this histograms, including possible under- and
  /// overflow bins.
  int GetNBins() const final {
    return Internal::TGetBinCount<sizeof...(AXISCONFIG) - 1,
       decltype(fAxes)>()(fAxes);
  }

  /// Add `w` to the bin at index `bin`.
  void AddBinContent(int bin, Weight_t w) {
    fContent[bin] += w;
  }

  std::tuple<AXISCONFIG...> fAxes; ///< The histogram's axes
  std::vector<PRECISION> fContent; ///< The histogram's bin content

public:
  THistImpl(STATISTICS statConfig, AXISCONFIG... axisArgs);

  /// Retrieve the fill function for this histogram implementation, to prevent
  /// the virtual function call for high-frequency fills.
  FillFunc_t GetFillFunc() const final { return (FillFunc_t)&THistImpl::Fill; }

  /// Get the axes of this histogram.
  const std::tuple<AXISCONFIG...>& GetAxes() const { return fAxes; }

  /// Normalized axes access, converting the actual axis to TAxisConfig
  TAxisView GetAxis(int iAxis) const final {
    return std::apply(Internal::GetAxisView<AXISCONFIG...>, fAxes)[iAxis];
  }


  /// Gets the bin index for coordinate `x`; returns -1 if there is no such bin,
  /// e.g. for axes without over / underflow but coordinate out of range.
  int GetBinIndex(const Coord_t& x) const final {
    TAxisBase::EFindStatus status = TAxisBase::EFindStatus::kValid;
    int ret = Internal::TGetBinIndex<DIMENSIONS - 1, THistImpl,
       decltype(fAxes), false>()(nullptr, fAxes, x, status);
    if (status != TAxisBase::EFindStatus::kValid)
      return -1;
    return ret;
  }

  /// Gets the bin index for coordinate `x`, growing the axes as needed and
  /// possible. Returns -1 if there is no such bin,
  /// e.g. for axes without over / underflow but coordinate out of range.
  int GetBinIndexAndGrow(const Coord_t& x) final {
    TAxisBase::EFindStatus status = TAxisBase::EFindStatus::kCanGrow;
    int ret = - 1;
    while (status == TAxisBase::EFindStatus::kCanGrow) {
      ret = Internal::TGetBinIndex<DIMENSIONS - 1, THistImpl, decltype(fAxes), true>()
         (this, fAxes, x, status);
    }
    return ret;
  }

  /// Get the center coordinate of the bin.
  Coord_t GetBinCenter(int binidx) const final {
    using FillBinCoord_t
      = Internal::FillBinCoord_t<DIMENSIONS - 1, Coord_t, decltype(fAxes)>;
    Coord_t coord;
    FillBinCoord_t()(coord, fAxes, Internal::EBinCoord::kBinCenter, binidx);
    return coord;
  }

  /// Get the coordinate of the low limit of the bin.
  Coord_t GetBinFrom(int binidx) const final {
    using FillBinCoord_t = Internal::FillBinCoord_t<DIMENSIONS - 1, Coord_t, decltype(fAxes)>;
    Coord_t coord;
    FillBinCoord_t()(coord, fAxes, Internal::EBinCoord::kBinFrom, binidx);
    return coord;
  }

  /// Get the coordinate of the high limit of the bin.
  Coord_t GetBinTo(int binidx) const final {
    using FillBinCoord_t =  Internal::FillBinCoord_t<DIMENSIONS - 1, Coord_t, decltype(fAxes)>;
    Coord_t coord;
    FillBinCoord_t()(coord, fAxes, Internal::EBinCoord::kBinTo, binidx);
    return coord;
  }

  /// Fill an array of `weightN` to the bins specified by coordinates `xN`.
  /// For each element `i`, the weight `weightN[i]` will be added to the bin
  /// at the coordinate `xN[i]`
  /// \note `xN` and `weightN` must have the same size!
  void FillN(const std::array_view<Coord_t> xN,
             const std::array_view<Weight_t> weightN) final {
#ifndef NDEBUG
    if (xN.size() != weightN.size()) {
      R__ERROR_HERE("HIST") << "Not the same number of points and weights!";
      return;
    }
#endif

    for (int i = 0; i < xN.size(); ++i) {
      Fill(xN[i], weightN[i]);
      STATISTICS::Fill(xN[i], weightN[i]);
    }
  }

  /// Fill an array of `weightN` to the bins specified by coordinates `xN`.
  /// For each element `i`, the weight `weightN[i]` will be added to the bin
  /// at the coordinate `xN[i]`
  void FillN(const std::array_view<Coord_t> xN) final {
    for (int i = 0; i < xN.size(); ++i) {
      Fill(xN[i]);
      STATISTICS::Fill(xN[i]);
    }
  }

  /// Return the uncertainties for the given bin.
  std::vector<double> GetBinUncertainties(int binidx) const final {
    return STATISTICS::GetBinUncertainties(binidx, *this);
  }

  /// Add a single weight `w` to the bin at coordinate `x`.
  void Fill(const Coord_t& x, Weight_t w = 1.) {
    int bin = GetBinIndexAndGrow(x);
    if (bin >= 0)
      AddBinContent(bin, w);
    STATISTICS::Fill(x, w);
  }

  /// Get the content of the bin at position `x`.
  PRECISION GetBinContent(const Coord_t& x) const {
    int bin = GetBinIndex(x);
    if (bin >= 0)
      return GetBinContent(bin);
    return 0.;
  }


  /// Get the content of the bin at bin index `binidx`.
  PRECISION GetBinContent(int binidx) const final {
    return fContent[binidx];
  }

  /// Get the begin() and end() for each axis.
  ///
  ///\param[in] withOverUnder - Whether the begin and end should contain over-
  /// or underflow. Ignored if the axis does not support over- / underflow.
  AxisIterRange_t<DIMENSIONS>
     GetRange(const std::array<Hist::EOverflow, DIMENSIONS>& withOverUnder) const final {
    std::array<std::array<TAxisBase::const_iterator, DIMENSIONS>, 2> ret;
    Internal::FillIterRange_t<DIMENSIONS - 1, decltype(fAxes)>()(ret, fAxes, withOverUnder);
    return ret;
  }

  /// Grow the axis number `iAxis` to fit the coordinate `x`.
  ///
  /// The histogram (conceptually) combines pairs of bins along this axis until
  /// `x` is within the range of the axis.
  /// The axis must support growing for this to work (e.g. a `TAxisGrow`).
  void GrowAxis(int /*iAxis*/, double /*x*/) {
    // TODO: Implement GrowAxis()
  }
};


template <int DIMENSIONS, class PRECISION, class STATISTICS, class... AXISCONFIG>
THistImpl<DIMENSIONS, PRECISION, STATISTICS, AXISCONFIG...>::THistImpl(STATISTICS statConfig, AXISCONFIG... axisArgs):
  STATISTICS(statConfig), fAxes{axisArgs...}, fContent(GetNBins())
{}

#if 0
// In principle we can also have a runtime version of THistImpl, that does not
// contain a tuple of concrete axis types but a vector of `TAxisConfig`.
template <int DIMENSIONS, class PRECISION>
class THistImplRuntime: public THistImplBase<DIMENSIONS, PRECISION> {
public:
  THistImplRuntime(std::array<TAxisConfig, DIMENSIONS>&& axisCfg);
};
#endif

} // namespace Detail
} // namespace ROOT

#endif
