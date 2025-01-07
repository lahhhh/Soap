#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"
#include "Enrichment.h"

class DifferentialAnalysis
{
public:

    enum class DataType : int { Plain = 0, Gene = 1, Peak = 2, GeneActivity = 3, ChromVAR = 4 , Cicero = 5};

    G_CLASS_FUNCTION_DEFAULT(DifferentialAnalysis);

    DifferentialAnalysis(const CustomMatrix& df) : mat_(df) {};

    DataType data_type_{ DataType::Plain };

    CustomMatrix mat_;

    SOAP_SUBMODULES(Enrichment);

    G_SET_IDENTIFIER("DifferentialAnalysis");
};

