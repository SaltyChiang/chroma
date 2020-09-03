#ifndef __SYSSOLVER_QSUNWAY_CLOVER_PARAMS_H__
#define __SYSSOLVER_QSUNWAY_CLOVER_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"

#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/invert/qsunway_solvers/enum_qsunway_io.h"
// #include "actions/ferm/invert/qsunway_solvers/qsunway_gcr_params.h"
#include <string>
#include "handle.h"

namespace Chroma
{
  struct SysSolverQSUNWAYCloverParams
  {
    SysSolverQSUNWAYCloverParams(XMLReader &xml, const std::string &path);
    SysSolverQSUNWAYCloverParams()
    {
      solverType = CG;
      // asymmetricP = false; //< Use asymmetric version of the linear operator
      axialGaugeP = false; //< Fix Axial Gauge?
      // SilentFailP = false; //< If set to true ignore lack of convergence. Default is 'loud'
      // RsdToleranceFactor = Real(10); //< Tolerate if the solution achived is better (less) than rsdToleranceFactor*RsdTarget
      // tuneDslashP = false ; //< v0.3 autotune feature
      verboseP = false;
      MROverParam = Real(1);
      // backup_invP = false;
      // dump_on_failP = false;
    };

    SysSolverQSUNWAYCloverParams(const SysSolverQSUNWAYCloverParams &p)
    {
      CloverParams = p.CloverParams;
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
      // backup_invP = p.backup_invP;
      // backup_inv_param = p.backup_inv_param;
      // dump_on_failP = p.dump_on_failP;
      MROverParam = p.MROverParam;
    }

    CloverFermActParams CloverParams;
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
    // bool innerParamsP;
    Real MROverParam;

    // XML for Backup Solver
    bool backup_invP;
    GroupXML_t backup_inv_param;
    bool dump_on_failP;
  };

  void read(XMLReader &xml, const std::string &path, SysSolverQSUNWAYCloverParams &p);

  void write(XMLWriter &xml, const std::string &path,
             const SysSolverQSUNWAYCloverParams &param);

} // namespace Chroma

#endif
