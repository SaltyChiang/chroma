#ifndef __SYSSOLVER_QSUNWAY_WILSON_PARAMS_H__
#define __SYSSOLVER_QSUNWAY_WILSON_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "actions/ferm/fermacts/wilson_fermact_params_w.h"
#include "actions/ferm/invert/qsunway_solvers/enum_qsunway_io.h"

#include <string>
#include "handle.h"

namespace Chroma
{

struct SysSolverQSUNWAYWilsonParams
{
  SysSolverQSUNWAYWilsonParams(XMLReader &xml, const std::string &path);
  SysSolverQSUNWAYWilsonParams()
  {
    solverType = MR;
    verboseP = false;
    // asymmetricP = false;           //< Use asymmetric version of the linear operator
    axialGaugeP = false;           //< Fix Axial Gauge?
    // SilentFailP = false;           //< If set to true ignore lack of convergence. Default is 'loud'
    // RsdToleranceFactor = Real(10); //< Tolerate if the solution achived is better (less) than rsdToleranceFactor*RsdTarget
    // tuneDslashP = false;           //< v0.3 autotune feature
    MROverParam = Real(1);
  };
  SysSolverQSUNWAYWilsonParams(const SysSolverQSUNWAYWilsonParams &p)
  {
    WilsonParams = p.WilsonParams;
    // AntiPeriodicT = p.AntiPeriodicT;
    MaxIter = p.MaxIter;
    RsdTarget = p.RsdTarget;
    // Delta = p.Delta;
    solverType = p.solverType;
    verboseP = p.verboseP;
    // asymmetricP = p.asymmetricP;
    axialGaugeP = p.axialGaugeP;
    // SilentFailP = p.SilentFailP;
    // RsdToleranceFactor = p.RsdToleranceFactor;
    // tuneDslashP = p.tuneDslashP;
    MROverParam = p.MROverParam;
  }

  WilsonFermActParams WilsonParams;
  // bool AntiPeriodicT;
  int MaxIter;
  Real RsdTarget;
  // Real Delta;
  QsunwaySolverType solverType;
  bool verboseP;
  // bool asymmetricP;
  bool axialGaugeP;
  // bool SilentFailP;
  // Real RsdToleranceFactor;
  // bool tuneDslashP;
  Real MROverParam;
};

void read(XMLReader &xml, const std::string &path, SysSolverQSUNWAYWilsonParams &p);

void write(XMLWriter &xml, const std::string &path,
           const SysSolverQSUNWAYWilsonParams &param);

} // namespace Chroma

#endif
