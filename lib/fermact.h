// -*- C++ -*-
// $Id: fermact.h,v 1.13 2004-12-07 17:11:50 bjoo Exp $

/*! @file
 * @brief Class structure for fermion actions
 */

#ifndef __fermact_h__
#define __fermact_h__

using namespace QDP;

#include "invtype.h"
#include "state.h"
#include "linearop.h"
#include "fermbc.h"

namespace Chroma
{
  //! Base class for quadratic matter actions (e.g., fermions)
  /*! @ingroup actions
   *
   * Supports creation and application for quadratic actions.
   * This is basically a foundry class with additional operations.
   *
   * The class holds info on the particulars of a bi-local action,
   * but it DOES NOT hold gauge fields. A specific dirac-operator
   * is a functional of the gauge fields, hence when a dirac-operator
   * is needed, it is created.
   *
   * The ConnectState holds gauge fields and whatever auxilliary info
   * is needed to create a specific dirac operator (linear operator)
   * on some background gauge field.
   *
   * The FermBC is the type of boundary conditions used for this action
   *
   * The linop and lmdagm functions create a linear operator on a 
 * fixed ConnectState
 *
 * The qprop solves for the full linear system undoing all preconditioning.
 * No direct functions are provided (or needed) to call a linear system
 * solver - that is a stand-alone function in the generic programming sense.
 *
 * dsdu provides the derivative of the action on a specific background
 * gauge field.
 *
 */
  template<typename T>
  class FermionAction
  {
  public:
    //! Given links, create the state needed for the linear operators
    /*! Default version uses a SimpleConnectState */
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u) const
      {
	multi1d<LatticeColorMatrix> u_tmp = u;
	getFermBC().modifyU(u_tmp);
	return new SimpleConnectState(u_tmp);
      }

    //! Given the links create a state with additional info held by the XMLReader
    /*! This is a default version, which ignores the reader and just calls the default
      createState() function */
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u,
					    XMLReader& reader,
					    const string& path) const
    {
      return createState(u);
    }

    //! Return the fermion BC object for this action
    /*! The user will supply the FermBC in a derived class */
    virtual const FermBC<T>& getFermBC() const = 0;

    //! Produce a linear operator for this action
    virtual const LinearOperator<T>* linOp(Handle<const ConnectState> state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    virtual const LinearOperator<T>* lMdagM(Handle<const ConnectState> state) const = 0;


    //! Produce the derivative of the Fermion Matrix with respect to 
    //! The Gauge field
    virtual const LinearOperator<T>* ldMdU(Handle<const ConnectState> state, const int mu) {
      QDPIO::cerr << "This is not yet implemented" << endl;
      QDP_abort(1);

      return 0x0;
    }

    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     *
     */
    virtual void qpropT(T& psi, 
			Handle<const ConnectState> state, 
			const T& chi, 
			const InvertParam_t& invParam,
			int& ncg_had) const;

    //! Compute quark propagator over simpler type
    /*! 
     * Solves  M.psi = chi
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     *
     */


    virtual void qprop(typename BaseType<T>::Type_t& psi, 
		       Handle<const ConnectState> state, 
		       const typename BaseType<T>::Type_t& chi, 
		       const InvertParam_t& invParam,
		       int& ncg_had) const;

    //! Compute dS_f/dU
    /*! Default version. Derived class should override this if needed. */
    virtual void dsdu(multi1d<LatticeColorMatrix>& result,
		      Handle<const ConnectState> state,
		      const T& psi) const
      {
	QDPIO::cerr << "FermionAction::dsdu not implemented" << endl;
	QDP_abort(1);
      }

  
    //! Virtual destructor to help with cleanup;
    virtual ~FermionAction() {}
  };


  //! Partial specialization for base class of quadratic matter actions (e.g., fermions)
  /*! @ingroup actions
   *
   * Supports creation and application for quadratic actions, specialized
   * to support arrays of matter fields
   */
  template<typename T>
  class FermionAction< multi1d<T> >
  {
  public:
    //! Expected length of array index
    virtual int size() const = 0;

    //! Given links, create the state needed for the linear operators
    /*! Default version uses a SimpleConnectState */
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u) const
      {
	multi1d<LatticeColorMatrix> u_tmp = u;
	getFermBC().modifyU(u_tmp);
	return new SimpleConnectState(u_tmp);
      }

    //! Given the links create a state with additional info held by the XMLReader
    /*! This is a default version, which ignores the reader and just calls the default
      createState() function */
    virtual const ConnectState* createState(const multi1d<LatticeColorMatrix>& u,
					    XMLReader& reader,
					    const string& path) const
    {
      return createState(u);
    }

    //! Return the fermion BC object for this action
    /*! The user will supply the FermBC in a derived class */
    virtual const FermBC< multi1d<T> >& getFermBC() const = 0;

    //! Produce a linear operator for this action
    virtual const LinearOperator< multi1d<T> >* linOp(Handle<const ConnectState> state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    virtual const LinearOperator< multi1d<T> >* lMdagM(Handle<const ConnectState> state) const = 0;

    virtual const LinearOperator< multi1d<T> >* ldMdU(Handle<const ConnectState> state, const int mu) {
      QDPIO::cerr << "This is not yet implemented" << endl;
      QDP_abort(1);

      return 0x0;
    }

    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     *
     */
    virtual void qpropT(multi1d<T>& psi, 
			Handle<const ConnectState> state, 
			const multi1d<T>& chi, 
			const InvertParam_t& invParam,
			int& ncg_had) const;

    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     *
     */
    virtual void qprop(T& psi, 
		       Handle<const ConnectState> state, 
		       const T& chi, 
		       const InvertParam_t& invParam,
		       int& ncg_had) const = 0;

    //! Compute dS_f/dU
    /*! Default version. Derived class should override this if needed. */
    virtual void dsdu(multi1d<LatticeColorMatrix>& result,
		      Handle<const ConnectState> state,
		      const multi1d<T>& psi) const
      {
	QDPIO::cerr << "FermionAction::dsdu not implemented" << endl;
	QDP_abort(1);
      }

    //! Virtual destructor to help with cleanup;
    virtual ~FermionAction() {}
  };


  //! Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Wilson-like fermion actions
   */

  template<typename T>
  class WilsonTypeFermAct : public FermionAction<T>
  {
  public:

    //! Produce a hermitian version of the linear operator
    virtual const LinearOperator<T>* gamma5HermLinOp(Handle<const ConnectState> state) const = 0;

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * Provides a default version
     *
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param xml_out  diagnostic output ( Modify )
     * \param state    gauge connection state ( Read )
     * \param invParam inverter parameters ( Read )
     * \param nonRelProp compute only a non-relativistic prop ( Read )
     * \param ncg_had  number of solver iterations ( Write )
     */
    virtual void quarkProp4(LatticePropagator& q_sol,   // Oops, need to make propagator type more general
			    XMLWriter& xml_out,
			    const LatticePropagator& q_src,
			    Handle<const ConnectState> state,
			    const InvertParam_t& invParam,
			    bool nonRelProp,
			    int& ncg_had);
  };

  //! Unpreconditioned Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Unpreconditioned like Wilson-like fermion actions
   */
  template<typename T>
  class UnprecWilsonTypeFermAct : public WilsonTypeFermAct<T>
  {
  public:
#if 0
    // This does not appear to be needed
    //! Compute quark propagator
    void qprop(T& psi, 
	       Handle<const ConnectState> state, 
	       const T& chi, 
	       const InvertParam_t& invParam,
	       int& ncg_had) const;
#endif
  };


  //! Even-odd preconditioned Wilson-like fermion actions
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions
   */
  template<typename T>
  class EvenOddPrecWilsonTypeFermAct : public WilsonTypeFermAct<T>
  {
  public:
    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecLinearOperator<T>* linOp(Handle<const ConnectState> state) const = 0;

    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * Provides a default version
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     */
    virtual void qpropT(T& psi, 
			Handle<const ConnectState> state, 
			const T& chi, 
			const InvertParam_t& invParam,
			int& ncg_had) const;

    //! Compute quark propagator over simpler type
    /*! 
     * Solves  M.psi = chi
     *
     * Provides a default version
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     */
    virtual void qprop(typename BaseType<T>::Type_t& psi, 
		       Handle<const ConnectState> state, 
		       const typename BaseType<T>::Type_t& chi, 
		       const InvertParam_t& invParam,
		       int& ncg_had) const;
  
  };


  //! Partial specialization of even-odd preconditioned Wilson-like fermion actions over arrays
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Wilson-like fermion actions.
   * Here, use arrays of matter fields.
 */
  template<typename T>
  class EvenOddPrecWilsonTypeFermAct< multi1d<T> > : public WilsonTypeFermAct< multi1d<T> >
  {
  public:


    //! Override to produce an even-odd prec. linear operator for this action
    /*! Covariant return rule - override base class function */
    virtual const EvenOddPrecLinearOperator< multi1d<T> >* linOp(Handle<const ConnectState> state) const = 0;

    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * Provides a default version
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     */
    virtual void qpropT(multi1d<T>& psi, 
			Handle<const ConnectState> state, 
			const multi1d<T>& chi, 
			const InvertParam_t& invParam,
			int& ncg_had) const;

    //! Compute quark propagator over simpler type
    /*! 
     * Solves  M.psi = chi
     *
     * Provides a default version
     *
     * \param psi      quark propagator ( Modify )
     * \param state    gauge connection state ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     */
    virtual void qprop(typename BaseType<T>::Type_t& psi, 
		       Handle<const ConnectState> state, 
		       const typename BaseType<T>::Type_t& chi, 
		       const InvertParam_t& invParam,
		       int& ncg_had) const = 0;
  
  };


  //! Staggered-like fermion actions
  /*! @ingroup actions
   *
   * Staggered-like fermion actions
   */
  template<typename T>
  class StaggeredTypeFermAct : public FermionAction<T>
  {
  public:
  };


  //! Even-odd preconditioned Staggered-like fermion actions
  /*! @ingroup actions
   *
   * Even-odd preconditioned like Staggered-like fermion actions
   */
  template<typename T>
  class EvenOddStaggeredTypeFermAct : public StaggeredTypeFermAct<T>
  {
  public:


    //! Provide a default version of qprop
    void qprop(LatticeStaggeredFermion& psi,
	       Handle<const ConnectState> state,
	       const LatticeStaggeredFermion& chi,
	       const InvertParam_t& invParam,
	       int& ncg_had);

    //! Return the quark mass
    virtual const Real  getQuarkMass() const = 0;

    //! Produce a linear operator for this action
    /*! NOTE: maybe this should be abstracted to a foundry class object */
    virtual const EvenOddLinearOperator<T>* linOp(Handle<const ConnectState> state) const = 0;

    //! Produce a linear operator M^dag.M for this action
    /*! NOTE: maybe this should be abstracted to a foundry class object */
    virtual const LinearOperator<T>* lMdagM(Handle<const ConnectState> state) const = 0;

    //! Virtual destructor to help with cleanup;
    virtual ~EvenOddStaggeredTypeFermAct() {}

  };

}

using namespace Chroma;

#endif
