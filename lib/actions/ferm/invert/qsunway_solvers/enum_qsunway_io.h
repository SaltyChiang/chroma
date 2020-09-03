#ifndef enum_qsunway_io_h
#define enum_qsunway_io_h

#include "chromabase.h"
#include <string>
#include "singleton.h"
#include "io/enum_io/enum_type_map.h"

namespace Chroma
{
  /*!
   * Types and structures
   *
   * \ingroup io
   *
   * @{
   */

  //! Qsunway Solver type
  enum QsunwaySolverType
  {
    CG,
    BICGSTAB,
    GCR,
    MR
  };

  namespace QsunwaySolverTypeEnv
  {
    extern const std::string typeIDString;
    extern bool registered;
    bool registerAll(void); // Forward declaration
  }                         // namespace QsunwaySolverTypeEnv

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<QsunwaySolverType>> theQsunwaySolverTypeMap;

  // Reader and writer

  //! Read an QsunwaySolverType enum
  void read(XMLReader &r, const std::string &path, QsunwaySolverType &t);

  //! Write an QsunwaySolverType enum
  void write(XMLWriter &w, const std::string &path, const QsunwaySolverType &t);

  // enum QsunwaySchwarzMethod
  // {
  //   ADDITIVE_SCHWARZ,
  //   MULTIPLICATIVE_SCHWARZ
  // };

  // namespace QsunwaySchwarzMethodEnv
  // {
  //   extern const std::string typeIDString;
  //   extern bool registered;
  //   bool registerAll(void); // Forward declaration
  // }                         // namespace QsunwaySchwarzMethodEnv

  // // A singleton to hold the typemap
  // typedef SingletonHolder<EnumTypeMap<QsunwaySchwarzMethod>> theQsunwaySchwarzMethodMap;

  // // Reader and writer

  // //! Read an QsunwaySchwarzMethod enum
  // void read(XMLReader &r, const std::string &path, QsunwaySchwarzMethod &t);

  // //! Write an QsunwaySchwarzMethod enum
  // void write(XMLWriter &w, const std::string &path, const QsunwaySchwarzMethod &t);

  /*! @} */ // end of group io

} // namespace Chroma

#endif
