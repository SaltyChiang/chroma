#include "actions/ferm/invert/qsunway_solvers/syssolver_qsunway_wilson_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "chroma_config.h"

using namespace QDP;

namespace Chroma
{

  SysSolverQSUNWAYWilsonParams::SysSolverQSUNWAYWilsonParams(XMLReader &xml,
                                                             const std::string &path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "WilsonParams", WilsonParams);
    // read(paramtop, "AntiPeriodicT", AntiPeriodicT);

    // read(paramtop, "Delta", Delta);

    if (paramtop.count("SolverType") > 0)
    {
      read(paramtop, "SolverType", solverType);
    }
    else
    {
      solverType = MR;
    }

    // if (paramtop.count("AsymmetricLinop") > 0)
    // {
    //   read(paramtop, "AsymmetricLinop", asymmetricP);
    // }
    // else
    // {
    //   asymmetricP = false; // Symmetric is default
    // }

    if (paramtop.count("Verbose") > 0)
    {
      read(paramtop, "Verbose", verboseP);
    }
    else
    {
      verboseP = false;
    }

    if (paramtop.count("AxialGaugeFix") > 0)
    {
      read(paramtop, "AxialGaugeFix", axialGaugeP);
    }
    else
    {
      axialGaugeP = false;
    }

    // if (paramtop.count("SilentFail") > 0)
    // {
    //   read(paramtop, "SilentFail", SilentFailP);
    // }
    // else
    // {
    //   SilentFailP = false;
    // }

    // if (paramtop.count("RsdToleranceFactor") > 0)
    // {
    //   read(paramtop, "RsdToleranceFactor", RsdToleranceFactor);
    // }
    // else
    // {
    //   RsdToleranceFactor = Real(10); // Tolerate an order of magnitude difference by default.
    // }

    // if (paramtop.count("AutotuneDslash") > 0)
    // {
    //   read(paramtop, "AutotuneDslash", tuneDslashP);
    // }
    // else
    // {
    //   tuneDslashP = false;
    // }
    // QDPIO::cout << "tuneDslasP = " << tuneDslashP << std::endl;

    // if (paramtop.count("GCRInnerParams") > 0)
    // {
    //   innerParams = new GCRInnerSolverParams(paramtop, "./GCRInnerParams");
    //   innerParamsP = true;
    // }
    // else
    // {
    //   innerParamsP = false;
    // }

    if (paramtop.count("MROverParam") > 0)
    {
      read(paramtop, "MROverParam", MROverParam);
    }
    else
    {
      MROverParam = Real(1);
    }
  }

  void read(XMLReader &xml, const std::string &path,
            SysSolverQSUNWAYWilsonParams &p)
  {
    SysSolverQSUNWAYWilsonParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter &xml, const std::string &path,
             const SysSolverQSUNWAYWilsonParams &p)
  {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "WilsonParams", p.WilsonParams);
    // write(xml, "AntiPeriodicT", p.AntiPeriodicT);
    // write(xml, "Delta", p.Delta);
    write(xml, "SolverType", p.solverType);
    write(xml, "Verbose", p.verboseP);
    // write(xml, "AsymmetricLinop", p.asymmetricP);
    write(xml, "AxialGaugeFix", p.axialGaugeP);
    // write(xml, "SilentFail", p.SilentFailP);
    // write(xml, "RsdToleranceFactor", p.RsdToleranceFactor);
    // write(xml, "AutotuneDslash", p.tuneDslashP);
    // if (p.innerParamsP)
    // {
    //   write(xml, "GCRInnerParams", *(p.innerParams));
    // }
    write(xml, "MROverParam", p.MROverParam);

    pop(xml);
  }

} // namespace Chroma
