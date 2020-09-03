// -*- C++ -*-
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by BiCGStab
 */

#ifndef __syssolver_linop_qsunway_clover_h__
#define __syssolver_linop_qsunway_clover_h__

#include "chroma_config.h"

#ifdef BUILD_QSUNWAY
extern "C"
{
#include <qsunway.h>
}

#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/fermbcs/simple_fermbc.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/invert/qsunway_solvers/syssolver_qsunway_clover_params.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "meas/gfix/temporal_gauge.h"
#include "io/aniso_io.h"
#include <string>

#include "util/gauge/reunit.h"

namespace Chroma
{

  //! Richardson system solver namespace
  namespace LinOpSysSolverQSUNWAYCloverEnv
  {
    //! Register the syssolver
    bool registerAll();
  } // namespace LinOpSysSolverQSUNWAYCloverEnv

  //! Solve a Clover Fermion System using the QSUNWAY inverter
  /*! \ingroup invert
 *** WARNING THIS SOLVER WORKS FOR Clover FERMIONS ONLY ***
   */

  class LinOpSysSolverQSUNWAYClover : public LinOpSystemSolver<LatticeFermion>
  {
  public:
    typedef LatticeFermion T;
    typedef LatticeColorMatrix U;
    typedef multi1d<LatticeColorMatrix> Q;

    typedef LatticeFermionF TF;
    typedef LatticeColorMatrixF UF;
    typedef multi1d<LatticeColorMatrixF> QF;

    typedef LatticeFermionF TD;
    typedef LatticeColorMatrixF UD;
    typedef multi1d<LatticeColorMatrixF> QD;

    typedef WordType<T>::Type_t REALT;
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverQSUNWAYClover(Handle<LinearOperator<T>> A_,
                                Handle<FermState<T, Q, Q>> state_,
                                const SysSolverQSUNWAYCloverParams &invParam_) : A(A_), invParam(invParam_), clov(new CloverTermT<T, U>::Type_t())
    {
      QDPIO::cout << "LinOpSysSolverQSUNWAYClover:" << std::endl;

      // FOLLOWING INITIALIZATION in test QSUNWAY program

      // 1) work out cpu_prec, cuda_prec, cuda_prec_sloppy
      int s = sizeof(WordType<T>::Type_t);
      if (s == 4)
      {
        cpu_prec = QSUNWAY_SINGLE_PRECISION;
      }
      else
      {
        cpu_prec = QSUNWAY_DOUBLE_PRECISION;
      }

      // 2) pull 'new; GAUGE and Invert params
      q_gauge_param = newQsunwayGaugeParam();
      qsunway_inv_param = newQsunwayInvertParam();

      // 3) set lattice size
      const multi1d<int> &latdims = Layout::subgridLattSize();

      if ((latdims[0] != 16) || (latdims[1] != 16) || (latdims[2] != 16) || (latdims[3] != 16))
      {
        QDPIO::cerr << "Subgrid lattice size for qsunway must be 16*16*16*16, please check your number of nodes and size fo the full grid." << std::endl;
        QDP_abort(1);
      }

      const multi1d<int> &griddims = Layout::getLogicalCoordFrom(Layout::numNodes() - 1);

      q_gauge_param.X[0] = griddims[3] + 1;
      q_gauge_param.X[1] = griddims[2] + 1;
      q_gauge_param.X[2] = griddims[1] + 1;
      q_gauge_param.X[3] = griddims[0] + 1;

      // 4) - deferred (anisotropy)

      // 5) - set QSUNWAY_WILSON_LINKS, QSUNWAY_GAUGE_ORDER
      // q_gauge_param.type = QSUNWAY_WILSON_LINKS;
      // q_gauge_param.gauge_order = QSUNWAY_QDP_GAUGE_ORDER; // gauge[mu], p

      // 6) - set t_boundary
      // Convention: BC has to be applied already
      // This flag just tells QSUNWAY that this is so,
      // so that QSUNWAY can take care in the reconstruct
      // if (invParam.AntiPeriodicT)
      // {
      //   q_gauge_param.t_boundary = QSUNWAY_ANTI_PERIODIC_T;
      // }
      // else
      // {
      //   q_gauge_param.t_boundary = QSUNWAY_PERIODIC_T;
      // }

      // Gauge fixing:

      // These are the links
      // They may be smeared and the BC's may be applied
      Q links_single(Nd);

      // Now downcast to single prec fields.
      for (int mu = 0; mu < Nd; mu++)
      {
        links_single[mu] = (state_->getLinks())[mu];
      }

      // GaugeFix
      if (invParam.axialGaugeP)
      {
        QDPIO::cout << "Fixing Temporal Gauge" << std::endl;
        temporalGauge(links_single, GFixMat, Nd - 1);
        for (int mu = 0; mu < Nd; mu++)
        {
          links_single[mu] = GFixMat * (state_->getLinks())[mu] * adj(shift(GFixMat, FORWARD, mu));
        }
        q_gauge_param.gauge_fix = QSUNWAY_GAUGE_FIXED_YES;
      }
      else
      {
        // No GaugeFix
        q_gauge_param.gauge_fix = QSUNWAY_GAUGE_FIXED_NO; // No Gfix yet
      }

      // deferred 4) Gauge Anisotropy
      const AnisoParam_t &aniso = invParam.CloverParams.anisoParam;
      if (aniso.anisoP)
      { // Anisotropic case
        Real gamma_f = aniso.xi_0 / aniso.nu;
        q_gauge_param.anisotropy = toDouble(gamma_f);
      }
      else
      {
        q_gauge_param.anisotropy = 1.0;
      }

      // MAKE FSTATE BEFORE RESCALING links_single
      // Because the clover term expects the unrescaled links...
      Handle<FermState<T, Q, Q>> fstate(new PeriodicFermState<T, Q, Q>(links_single));

      // if (aniso.anisoP)
      // { // Anisotropic case
      //   multi1d<Real> cf = makeFermCoeffs(aniso);
      //   for (int mu = 0; mu < Nd; mu++)
      //   {
      //     links_single[mu] *= cf[mu];
      //   }
      // }

      // Now onto the inv param:
      // Dslash type
      qsunway_inv_param.dslash_type = QSUNWAY_CLOVER_WILSON_DSLASH;

      // Invert type:
      switch (invParam.solverType)
      {
      case CG:
        qsunway_inv_param.inv_type = QSUNWAY_CG_INVERTER;
        solver_string = "CG";
        qsunway_inv_param.mr_over = 0.0;
        break;
      case BICGSTAB:
        qsunway_inv_param.inv_type = QSUNWAY_BICGSTAB_INVERTER;
        solver_string = "BICGSTAB";
        QDPIO::cerr << solver_string << " inverter is not implemented yet." << std::endl;
        QDP_abort(1);
        break;
      case GCR:
        qsunway_inv_param.inv_type = QSUNWAY_GCR_INVERTER;
        solver_string = "GCR";
        QDPIO::cerr << solver_string << " inverter is not implemented yet." << std::endl;
        QDP_abort(1);
        break;
      case MR:
        qsunway_inv_param.inv_type = QSUNWAY_MR_INVERTER;
        solver_string = "MR";
        qsunway_inv_param.mr_over = toDouble(invParam.MROverParam);
        break;
      default:
        QDPIO::cerr << "Unknown Solver type" << std::endl;
        QDP_abort(1);
        break;
      }

      // Mass

      // Fiendish idea from Ron. Set the kappa=1/2 and use
      // unmodified clover term, and ask for Kappa normalization
      // This should give us A - (1/2)D as the unpreconditioned operator
      // and probabl 1 - {1/4} A^{-1} D A^{-1} D as the preconditioned
      // op. Apart from the A_oo stuff on the antisymmetric we have
      // nothing to do...
      qsunway_inv_param.kappa = 0.5;

      // FIXME: If we want QSUNWAY to compute the clover coeff, we need to be able to deal
      // with awfuless of anisotropy
      // The value below is a dummy one.
      // qsunway_inv_param.clover_coeff = 1.0; // Dummy tree level value. Not used
      qsunway_inv_param.tol = toDouble(invParam.RsdTarget);
      qsunway_inv_param.maxiter = invParam.MaxIter;
      // qsunway_inv_param.reliable_delta = toDouble(invParam.Delta);

      // Solution type
      // qsunway_inv_param.solution_type = QSUNWAY_MATPC_SOLUTION;

      // Solve type
      // switch (invParam.solverType)
      // {
      // case CG:
      //   qsunway_inv_param.solve_type = QSUNWAY_NORMOP_PC_SOLVE;
      //   break;
      // case BICGSTAB:
      //   qsunway_inv_param.solve_type = QSUNWAY_DIRECT_PC_SOLVE;
      //   break;
      // case GCR:
      //   qsunway_inv_param.solve_type = QSUNWAY_DIRECT_PC_SOLVE;
      //   break;

      // case MR:
      //   qsunway_inv_param.solve_type = QSUNWAY_DIRECT_PC_SOLVE;
      //   break;

      // default:
      //   qsunway_inv_param.solve_type = QSUNWAY_NORMOP_PC_SOLVE;

      //   break;
      // }

      // if (invParam.asymmetricP)
      // {
      //   QDPIO::cout << "Using Asymmetric Linop: A_oo - D A^{-1}_ee D" << std::endl;
      //   qsunway_inv_param.matpc_type = QSUNWAY_MATPC_ODD_ODD_ASYMMETRIC;
      // }
      // else
      // {
      //   QDPIO::cout << "Using Symmetric Linop: 1 - A^{-1}_oo D A^{-1}_ee D" << std::endl;
      //   qsunway_inv_param.matpc_type = QSUNWAY_MATPC_ODD_ODD;
      // }

      qsunway_inv_param.dagger = QSUNWAY_DAG_NO;
      // qsunway_inv_param.mass_normalization = QSUNWAY_KAPPA_NORMALIZATION;

      // qsunway_inv_param.cpu_prec = cpu_prec;
      // qsunway_inv_param.cuda_prec = gpu_prec;
      // qsunway_inv_param.cuda_prec_sloppy = gpu_half_prec;
      // qsunway_inv_param.preserve_source = QSUNWAY_PRESERVE_SOURCE_NO;
      // qsunway_inv_param.gamma_basis = QSUNWAY_DEGRAND_ROSSI_GAMMA_BASIS;

      // qsunway_inv_param.dirac_order = QSUNWAY_DIRAC_ORDER;

      // // Autotuning
      // if (invParam.tuneDslashP)
      // {
      //   QDPIO::cout << "Enabling Dslash Autotuning" << std::endl;

      //   qsunway_inv_param.tune = QSUNWAY_TUNE_YES;
      // }
      // else
      // {
      //   QDPIO::cout << "Disabling Dslash Autotuning" << std::endl;

      //   qsunway_inv_param.tune = QSUNWAY_TUNE_NO;
      // }

      // if (invParam.innerParamsP)
      // {
      //   QDPIO::cout << "Setting inner solver params" << std::endl;
      //   // Dereference handle
      //   GCRInnerSolverParams ip = *(invParam.innerParams);
      // }
      // else
      // {
      //   QDPIO::cout << "Setting Precondition stuff to defaults for not using" << std::endl;
      //   qsunway_inv_param.inv_type_precondition = QSUNWAY_INVALID_INVERTER;
      //   qsunway_inv_param.tol_precondition = 1.0e-1;
      //   qsunway_inv_param.maxiter_precondition = 1000;
      //   qsunway_inv_param.verbosity_precondition = QSUNWAY_SILENT;
      //   qsunway_inv_param.gcrNkrylov = 1;
      // }

      // qsunway_inv_param.clover_order = QSUNWAY_PACKED_CLOVER_ORDER;

      if (invParam.verboseP)
      {
        qsunway_inv_param.verbosity = QSUNWAY_VERBOSE;
      }
      else
      {
        qsunway_inv_param.verbosity = QSUNWAY_SUMMARIZE;
      }

      // Set up the links
      void *gauge[4];

      for (int mu = 0; mu < Nd; mu++)
      {
        gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0, 0).real());
      }

      loadGaugeQsunway((void *)gauge, &q_gauge_param);

      //      Setup the clover term...
      QDPIO::cout << "Creating CloverTerm" << std::endl;
      clov->create(fstate, invParam_.CloverParams);

      void *cloverPtr = (void *)&(clov->getTriBuffer()[0].diag[0][0]);
      loadCloverQsunway(cloverPtr, &qsunway_inv_param);
    }

    //! Destructor is automatic
    ~LinOpSysSolverQSUNWAYClover()
    {
      QDPIO::cout << "Destructing" << std::endl;
      freeGaugeQsunway();
      freeCloverQsunway();
    }

    //! Return the subset on which the operator acts
    const Subset &subset() const { return A->subset(); }

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator()(T &psi, const T &chi) const
    {
      SystemSolverResults_t res;

      START_CODE();
      StopWatch swatch;
      swatch.start();

      //    T MdagChi;

      // This is a CGNE. So create new RHS
      //      (*A)(MdagChi, chi, MINUS);
      // Handle< LinearOperator<T> > MM(new MdagMLinOp<T>(A));
      if (invParam.axialGaugeP)
      {
        T g_chi, g_psi;

        // Gauge Fix source and initial guess
        QDPIO::cout << "Gauge Fixing source and initial guess" << std::endl;
        g_chi[rb[1]] = GFixMat * chi;
        g_psi[rb[1]] = GFixMat * psi;
        QDPIO::cout << "Solving" << std::endl;
        res = qsunwayInvert(g_chi,
                            g_psi);
        QDPIO::cout << "Untransforming solution." << std::endl;
        psi[rb[1]] = adj(GFixMat) * g_psi;
      }
      else
      {
        res = qsunwayInvert(chi,
                            psi);
      }

      swatch.stop();
      double time = swatch.getTimeInSeconds();

      {
        T r;
        r[A->subset()] = chi;
        T tmp;
        (*A)(tmp, psi, PLUS);
        r[A->subset()] -= tmp;
        res.resid = sqrt(norm2(r, A->subset()));
      }

      Double rel_resid = res.resid / sqrt(norm2(chi, A->subset()));

      QDPIO::cout << "QSUNWAY_" << solver_string << "_CLOVER_SOLVER: " << res.n_count << " iterations. Rsd = " << res.resid << " Relative Rsd = " << rel_resid << std::endl;

      // // Convergence Check/Blow Up
      // if (!invParam.SilentFailP)
      // {
      //   if (toBool(rel_resid > invParam.RsdToleranceFactor * invParam.RsdTarget))
      //   {
      //     QDPIO::cerr << "ERROR: QSUNWAY Solver residuum is outside tolerance: QSUNWAY resid=" << rel_resid << " Desired =" << invParam.RsdTarget << " Max Tolerated = " << invParam.RsdToleranceFactor * invParam.RsdTarget << std::endl;
      //     QDP_abort(1);
      //   }
      // }

      END_CODE();
      return res;
    }

  private:
    // Hide default constructor
    LinOpSysSolverQSUNWAYClover() {}

#if 1
    Q links_orig;
#endif

    U GFixMat;
    QsunwayPrecision_s cpu_prec;
    QsunwayPrecision_s gpu_prec;
    QsunwayPrecision_s gpu_half_prec;

    Handle<LinearOperator<T>> A;
    const SysSolverQSUNWAYCloverParams invParam;
    QsunwayGaugeParam q_gauge_param;
    QsunwayInvertParam qsunway_inv_param;

    Handle<CloverTermT<T, U>::Type_t> clov;

    SystemSolverResults_t qsunwayInvert(const T &chi_s,
                                        T &psi_s) const;

    std::string solver_string;
  };

} // namespace Chroma

#endif // BUILD_QSUNWAY
#endif
