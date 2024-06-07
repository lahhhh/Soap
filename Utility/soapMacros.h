#pragma once

#define SOAP_INLINE __forceinline

#define SOAP_EDITION "soap_1_0"

#include "soapFileName.h"
#include "soapMetadataName.h"
#include "soapVariableName.h"
#include "soapMacroFunction.h"
#include "SoapPhrase.h"

#define ITEM_IDENTIFIER_1 1997
#define ITEM_IDENTIFIER_2 1114

#define BARCODES_FILE_NAME_10X "barcodes.tsv.gz"
#define FEATURE_FILE_NAME_10X "features.tsv.gz"
#define MATRIX_FILE_NAME_10X "matrix.mtx.gz"

#define SOAP_DELIMITER "$%SOAP_DELIMITER%$"
#define SOAP_DELIMITER2 "$%SOAP_DELIMITER2%$"
#define SOAP_DELIMITER3 "$%SOAP_DELIMITER3%$"

#define SWITCH_ACCEPT "%SWITCH_ACCEPT%"
#define SWITCH_REJECT "%SWITCH_REJECT%"

#define SOAP_SUBMODULES(x) \
std::map<QString, x> x##s_;

#define DATA_SUBMODULES(x) this->data()->x##s_

#define SUBMODULES(x, mod) (x).mod##s_