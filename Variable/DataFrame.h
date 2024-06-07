#pragma once

#include "Identifier.h"

#include "CustomMatrix.h"

class DataFrame
{
public:

    enum class DataType : int { Plain = 0 };

    G_CLASS_FUNCTION_DEFAULT(DataFrame);

    template<typename Mat> requires (!std::same_as<std::decay_t<Mat>, DataFrame>)
    DataFrame(Mat&& df) : mat_(std::forward<Mat>(df)) {}

    DataType data_type_{ DataType::Plain };

    CustomMatrix mat_;

    G_SET_IDENTIFIER("DataFrame");
};
