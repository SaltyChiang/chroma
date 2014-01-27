// -*- C++ -*-
// $Id: lwldslash_base_w.h,v 3.4 2009-04-17 02:05:33 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#ifndef __lwldslash_base_h__
#define __lwldslash_base_h__

#include "linearop.h"

namespace Chroma 
{ 
  //! General Wilson-Dirac dslash
  /*!
   * \ingroup linop
   *
   * DSLASH
   *
   * This routine is specific to Wilson fermions!
   *
   * Description:
   *
   * This routine applies the operator D' to Psi, putting the result in Chi.
   *
   *	       Nd-1
   *	       ---
   *	       \
   *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
   *	       /    mu			  mu
   *	       ---
   *	       mu=0
   *
   *	             Nd-1
   *	             ---
   *	             \    +
   *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *	             /    mu			   mu
   *	             ---
   *	             mu=0
   *
   */
  template <typename T, typename P, typename Q>
  class WilsonDslashBase : 
    public DslashLinearOperator<T,P,Q>
  {
  public:

    //! No real need for cleanup here
    virtual ~WilsonDslashBase() {}

    //! Subset is all here
    const Subset& subset() const {return all;}

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign) const;
    virtual void derivAdd(P& ds_u,
		       const T& chi, const T& psi,
		       enum PlusMinus isign) const;

    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign, int cb) const ;
    virtual void derivAdd(P& ds_u,
		       const T& chi, const T& psi,
		       enum PlusMinus isign, int cb) const ;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  protected:
    //! Get the anisotropy parameters
    virtual const multi1d<Real>& getCoeffs() const = 0;
  };

  template<typename T>
  class HalfFermionType{};

  template<>
  class HalfFermionType<LatticeFermionF>  { 
  public:
    typedef LatticeHalfFermionF Type_t;
  };

  template<>
  class HalfFermionType<LatticeFermionD>  { 
  public:
    typedef LatticeHalfFermionD Type_t;
  };

  template<typename T, typename P, typename Q>
  void
  WilsonDslashBase<T,P,Q>::derivAdd(P& ds_u,
			  const T& chi, const T& psi,
			  enum PlusMinus isign) const
  {
    START_CODE();

    derivAdd(ds_u, chi, psi, isign, 0);
    derivAdd(ds_u, chi, psi, isign, 1);

    END_CODE();
  }


  //! Take deriv of D
  /*!
   * \param chi     left vector                                 (Read)
   * \param psi     right vector                                (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  template<typename T, typename P, typename Q>
  void
  WilsonDslashBase<T,P,Q>::deriv(P& ds_u,
			  const T& chi, const T& psi, 
			  enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    deriv(ds_u, chi, psi, isign, 0);
    derivAdd(ds_u, chi, psi, isign, 1);

    END_CODE();
  }


  //! Take deriv of D
  /*! \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$  */
  template<typename T, typename P, typename Q>
  void 
  WilsonDslashBase<T,P,Q>::deriv(P& ds_u,
				 const T& chi, const T& psi, 
				 enum PlusMinus isign, int cb) const
  {
    START_CODE();

    ds_u.resize(Nd);

    const multi1d<Real>& anisoWeights = getCoeffs();

		T temp_ferm1, temp_ferm2;
		typename HalfFermionType<T>::Type_t tmp_h;
		P temp_mat(1);

    for(int mu = 0; mu < Nd; ++mu) 
    {

      switch (isign) 
      {
      case PLUS:
      {
	// Undaggered: Minus Projectors
	switch(mu) 
	{ 
	case 0:
	  tmp_h[rb[1-cb]] = spinProjectDir0Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir0Minus(tmp_h);
	  break;
	case 1:
	  tmp_h[rb[1-cb]] = spinProjectDir1Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir1Minus(tmp_h);
	  break;
	case 2:
	  tmp_h[rb[1-cb]] = spinProjectDir2Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir2Minus(tmp_h);
	  break;
	case 3:
	  tmp_h[rb[1-cb]] = spinProjectDir3Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir3Minus(tmp_h);
	  break;
	default:
	  break;
	};
	
      }
      break;

      case MINUS:
      {
	// Daggered: Plus Projectors
	switch(mu) 
	{
	case 0:
	  tmp_h[rb[1-cb]] = spinProjectDir0Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir0Plus(tmp_h);
	  break;
	case 1:
	  tmp_h[rb[1-cb]] = spinProjectDir1Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir1Plus(tmp_h);
	  break;
	case 2:
	  tmp_h[rb[1-cb]] = spinProjectDir2Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir2Plus(tmp_h);
	  break;
	case 3:
	  tmp_h[rb[1-cb]] = spinProjectDir3Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir3Plus(tmp_h);
	  break;
	default:
	  break;
	};
      }
      break;

      default:
	QDP_error_exit("unknown case");
      }

      // QDP Shifts the whole darn thing anyhow
      temp_ferm2 = shift(temp_ferm1, FORWARD, mu);
      // This step supposedly optimised in QDP++
      (temp_mat[0])[rb[cb]] = traceSpin(outerProduct(temp_ferm2,chi));
    
      // Just do the bit we need.
      ds_u[mu][rb[cb]] = anisoWeights[mu] * temp_mat[0];
      ds_u[mu][rb[1-cb]] = zero;    
    }
    (*this).getFermBC().zero(ds_u);

    END_CODE();
  }

  template<typename T, typename P, typename Q>
  void
  WilsonDslashBase<T,P,Q>::derivAdd(P& ds_u,
				 const T& chi, const T& psi,
				 enum PlusMinus isign, int cb) const
  {
    START_CODE();

    const multi1d<Real>& anisoWeights = getCoeffs();

    T temp_ferm1, temp_ferm2;
    typename HalfFermionType<T>::Type_t tmp_h;
    P temp_mat(1);

    for(int mu = 0; mu < Nd; ++mu)
    {
      switch (isign)
      {
      case PLUS:
      {
	// Undaggered: Minus Projectors
	switch(mu)
	{
	case 0:
	  tmp_h[rb[1-cb]] = spinProjectDir0Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir0Minus(tmp_h);
	  break;
	case 1:
	  tmp_h[rb[1-cb]] = spinProjectDir1Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir1Minus(tmp_h);
	  break;
	case 2:
	  tmp_h[rb[1-cb]] = spinProjectDir2Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir2Minus(tmp_h);
	  break;
	case 3:
	  tmp_h[rb[1-cb]] = spinProjectDir3Minus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir3Minus(tmp_h);
	  break;
	default:
	  break;
	};

      }
      break;

      case MINUS:
      {
	// Daggered: Plus Projectors
	switch(mu)
	{
	case 0:
	  tmp_h[rb[1-cb]] = spinProjectDir0Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir0Plus(tmp_h);
	  break;
	case 1:
	  tmp_h[rb[1-cb]] = spinProjectDir1Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir1Plus(tmp_h);
	  break;
	case 2:
	  tmp_h[rb[1-cb]] = spinProjectDir2Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir2Plus(tmp_h);
	  break;
	case 3:
	  tmp_h[rb[1-cb]] = spinProjectDir3Plus(psi);
	  temp_ferm1[rb[1-cb]] = spinReconstructDir3Plus(tmp_h);
	  break;
	default:
	  break;
	};
      }
      break;

      default:
	QDP_error_exit("unknown case");
      }

      temp_ferm2 = shift(temp_ferm1, FORWARD, mu);
      (temp_mat[0])[rb[cb]] = traceSpin(outerProduct(temp_ferm2,chi));

      ds_u[mu][rb[cb]] += anisoWeights[mu] * temp_mat[0];
    }
    (*this).getFermBC().zero(ds_u);

    END_CODE();
  }


  //! Return flops performed by the operator()
  template<typename T, typename P, typename Q>
  unsigned long 
  WilsonDslashBase<T,P,Q>::nFlops() const {return 1320;}

#if 0

 class WilsonDslashBaseF : 
    public DslashLinearOperator<LatticeFermionF,
				multi1d<LatticeColorMatrixF>,
				multi1d<LatticeColorMatrixF> >
  {
  public:
    typedef LatticeFermionF T;
    typedef multi1d<LatticeColorMatrixF> P;
    typedef multi1d<LatticeColorMatrixF> Q;

    //! No real need for cleanup here
    virtual ~WilsonDslashBaseF() {}

    //! Subset is all here
    const Subset& subset() const {return all;}

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign) const {
      QDPIO::cout << "Not implemented" << endl;
      QDP_abort(1);
    }

    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
     */
    virtual void deriv(P& ds_u, 
		       const T& chi, const T& psi, 
		       enum PlusMinus isign, int cb) const {
      QDPIO::cout << "Not Implemented" << endl;
    }

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  protected:
    //! Get the anisotropy parameters
    virtual const multi1d<Real>& getCoeffs() const = 0;
  };

#endif

} // End Namespace Chroma


#endif
