/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/qsunway_solvers/syssolver_qsunway_clover_params.h"
#include "actions/ferm/invert/qsunway_solvers/syssolver_linop_clover_qsunway_w.h"
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
  namespace LinOpSysSolverQSUNWAYCloverEnv
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QSUNWAY_CLOVER_INVERTER");

      //! Local registration flag
      bool registered = false;
    } // namespace

    LinOpSystemSolver<LatticeFermion> *createFerm(XMLReader &xml_in,
                                                  const std::string &path,
                                                  Handle<FermState<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix>>> state,

                                                  Handle<LinearOperator<LatticeFermion>> A)
    {
      return new LinOpSysSolverQSUNWAYClover(A, state, SysSolverQSUNWAYCloverParams(xml_in, path));
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
  } // namespace LinOpSysSolverQSUNWAYCloverEnv

  SystemSolverResults_t
  LinOpSysSolverQSUNWAYClover::qsunwayInvert(const T &chi_s,
                                             T &psi_s) const
  {

    SystemSolverResults_t ret;

    void *spinorIn = (void *)&(chi_s.elem(rb[1].start()).elem(0).elem(0).real());
    void *spinorOut = (void *)&(psi_s.elem(rb[1].start()).elem(0).elem(0).real());

    // Do the solve here
    StopWatch swatch1;
    swatch1.reset();
    swatch1.start();
    invertQsunway(spinorOut, spinorIn, (QsunwayInvertParam *)&qsunway_inv_param);
    swatch1.stop();

    // QDPIO::cout << "Cuda Space Required" << std::endl;
    // QDPIO::cout << "\t Spinor:" << qsunway_inv_param.spinorGiB << " GiB" << std::endl;
    // QDPIO::cout << "\t Gauge :" << q_gauge_param.gaugeGiB << " GiB" << std::endl;
    // QDPIO::cout << "\t InvClover :" << qsunway_inv_param.cloverGiB << " GiB" << std::endl;
    // QDPIO::cout << "QSUNWAY_" << solver_string << "_CLOVER_SOLVER: time=" << qsunway_inv_param.secs << " s";
    // QDPIO::cout << "\tPerformance=" << qsunway_inv_param.gflops / qsunway_inv_param.secs << " GFLOPS";
    QDPIO::cout << "\tTotal Time (incl. load gauge)=" << swatch1.getTimeInSeconds() << " s" << std::endl;

    ret.n_count = qsunway_inv_param.iter;

    return ret;
  }

} // namespace Chroma
