// $Id: prec_asqtad_fermact_s.cc,v 1.2 2003-12-10 16:20:59 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */
// NEW $Id: asqtad_fermact_s.cc 2003/11/12 steve

#include "chromabase.h"
//#include "actions/ferm/linop/asqtad_linop_s.h"
//#include "actions/ferm/fermacts/asqtad_fermact_s.h"
//#include "actions/ferm/linop/lmdagm_s.h"

// #include "actions/ferm/linop/prec_asq_mdagm_s.h"
#include "actions/ferm/linop/prec_asqtad_linop_s.h"
#include "actions/ferm/fermacts/prec_asqtad_fermact_s.h"


//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the entire lattice
 *
 * \param u_fat, u_triple 	 fat7 and triple links    (Read)
 * \u has already had KS phases multiplied in.
 */
const EvenOddPrecLinearOperator<LatticeFermion>* 
EvenOddPrecAsqtadFermAct::linOp(const ConnectState& state_) const
{
  const AsqtadConnectState<LatticeFermion>& state = 
    dynamic_cast<const AsqtadConnectState<LatticeFermion>&>(state_);

  return new EvenOddPrecAsqtadLinOp(state.getFatLinks(), state.getTripleLinks(), Mass);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the checkerboarded lattice
 *
 * \param u_fat, u_triple 	 fat7 and triple links	       (Read)
 */
/*
const EvenOddPrecLinearOperator<LatticeFermion>* 
EvenOddPrecAsqtadFermAct::lMdagM(const ConnectState& state_) const
{
  const AsqtadConnectState<LatticeFermion>& state = 
    dynamic_cast<const AsqtadConnectState<LatticeFermion>&>(state_);
  
  return new PrecAsqtadMdagM(state.getFatLinks(), state.getTripleLinks(), Mass);
}
*/
