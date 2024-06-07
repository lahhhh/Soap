#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"
#include "Enrichment.h"

class DifferentialAnalysis
{
public:

    enum class DataType : int { Plain = 0, Gene = 1, Peak = 2, GeneActivity = 3, ChromVAR = 4 , Cicero = 5};

    G_CLASS_FUNCTION_DEFAULT(DifferentialAnalysis);

    template<typename Mat> requires (!std::same_as<std::decay_t<Mat>, DifferentialAnalysis>)
        DifferentialAnalysis(Mat&& df) : mat_(std::forward<Mat>(df)) {}

    DataType data_type_{ DataType::Plain };

    CustomMatrix mat_;

    SOAP_SUBMODULES(Enrichment);

    G_SET_IDENTIFIER("DifferentialAnalysis");
};

