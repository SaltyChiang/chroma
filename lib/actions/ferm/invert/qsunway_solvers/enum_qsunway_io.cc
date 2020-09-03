// -*- C++ -*-
/*! \file
 *  \brief QSUNWAY enum readers 
 */

#include "actions/ferm/invert/qsunway_solvers/enum_qsunway_io.h"
#include <string>

namespace Chroma
{

  namespace QsunwaySolverTypeEnv
  {
    bool registerAll(void)
    {
      bool success;
      success = theQsunwaySolverTypeMap::Instance().registerPair(std::string("CG"), CG);
      success &= theQsunwaySolverTypeMap::Instance().registerPair(std::string("BICGSTAB"), BICGSTAB);
      success &= theQsunwaySolverTypeMap::Instance().registerPair(std::string("GCR"), GCR);
      success &= theQsunwaySolverTypeMap::Instance().registerPair(std::string("MR"), MR);
      return success;
    }
    const std::string typeIDString = "QsunwaySolverType";
    bool regisered = registerAll();
  }; // namespace QsunwaySolverTypeEnv

  //! Read an QsunwaySolverType enum
  void read(XMLReader &xml_in, const std::string &path, QsunwaySolverType &t)
  {
    theQsunwaySolverTypeMap::Instance().read(QsunwaySolverTypeEnv::typeIDString, xml_in, path, t);
  }

  //! Write an QsunwaySolverType enum
  void write(XMLWriter &xml_out, const std::string &path, const QsunwaySolverType &t)
  {
    theQsunwaySolverTypeMap::Instance().write(QsunwaySolverTypeEnv::typeIDString, xml_out, path, t);
  }

  // namespace QsunwaySchwarzMethodEnv
  // {
  //   bool registerAll(void)
  //   {
  //     bool success;
  //     success = theQsunwaySchwarzMethodMap::Instance().registerPair(std::string("ADDITIVE_SCHWARZ"), ADDITIVE_SCHWARZ);
  //     success &= theQsunwaySchwarzMethodMap::Instance().registerPair(std::string("MULTIPLICATIVE_SCHWARZ"), MULTIPLICATIVE_SCHWARZ);
  //     return success;
  //   }
  //   const std::string typeIDString = "QsunwaySchwarzMethod";
  //   bool regisered = registerAll();
  // }; // namespace QsunwaySchwarzMethodEnv

  // //! Read an QsunwaySolverType enum
  // void read(XMLReader &xml_in, const std::string &path, QsunwaySchwarzMethod &t)
  // {
  //   theQsunwaySchwarzMethodMap::Instance().read(QsunwaySchwarzMethodEnv::typeIDString, xml_in, path, t);
  // }

  // //! Write an QsunwaySolverType enum
  // void write(XMLWriter &xml_out, const std::string &path, const QsunwaySchwarzMethod &t)
  // {
  //   theQsunwaySchwarzMethodMap::Instance().write(QsunwaySchwarzMethodEnv::typeIDString, xml_out, path, t);
  // }

} // namespace Chroma
