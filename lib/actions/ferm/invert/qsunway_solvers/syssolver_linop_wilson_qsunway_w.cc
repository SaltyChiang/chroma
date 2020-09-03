/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/qsunway_solvers/syssolver_qsunway_wilson_params.h"
#include "actions/ferm/invert/qsunway_solvers/syssolver_linop_wilson_qsunway_w.h"
#include "io/aniso_io.h"

#include "handle.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/linop/lwldslash_w.h"
#include "meas/glue/mesplq.h"
// QSUNWAY Headers
extern "C"
{
#include <qsunway.h>
}

namespace Chroma
{
namespace LinOpSysSolverQSUNWAYWilsonEnv
{

//! Anonymous namespace
namespace
{
//! Name to be used
const std::string name("QSUNWAY_WILSON_INVERTER");

//! Local registration flag
bool registered = false;
} // namespace

LinOpSystemSolver<LatticeFermion> *createFerm(XMLReader &xml_in,
                                              const std::string &path,
                                              Handle<FermState<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>> state,

                                              Handle<LinearOperator<LatticeFermion>> A)
{
  return new LinOpSysSolverQSUNWAYWilson(A, state, SysSolverQSUNWAYWilsonParams(xml_in, path));
}

//! Register all the factories
bool registerAll()
{
  bool success = true;
  if (!registered)
  {
    success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
    registered = true;
  }
  return success;
}
} // namespace LinOpSysSolverQSUNWAYWilsonEnv

SystemSolverResults_t
LinOpSysSolverQSUNWAYWilson::qsunwayInvert(const T &chi_s,
                                           T &psi_s) const
{

  SystemSolverResults_t ret;

  void *spinorIn;

  T mod_chi;

  spinorIn = (void *)&(chi_s.elem(all.start()).elem(0).elem(0).real());

  void *spinorOut = (void *)&(psi_s.elem(all.start()).elem(0).elem(0).real());

  // Do the solve here
  StopWatch swatch1;
  swatch1.reset();
  swatch1.start();
  invertQsunway(spinorOut, spinorIn, (QsunwayInvertParam *)&qsunway_inv_param);
  swatch1.stop();

  // QDPIO::cout << "Cuda Space Required" << std::endl;
  // QDPIO::cout << "\t Spinor:" << qsunway_inv_param.spinorGiB << " GiB" << std::endl;
  // QDPIO::cout << "\t Gauge :" << q_gauge_param.gaugeGiB << " GiB" << std::endl;
  // QDPIO::cout << "QSUNWAY_" << solver_string << "_WILSON_SOLVER: time=" << qsunway_inv_param.secs << " s";
  // QDPIO::cout << "\tPerformance=" << qsunway_inv_param.gflops / qsunway_inv_param.secs << " GFLOPS";
  QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() << " s" << std::endl;

  ret.n_count = qsunway_inv_param.iter;

  return ret;
}

} // namespace Chroma
