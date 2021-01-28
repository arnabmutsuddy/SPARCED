#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <array>

#include "x_rdata.h"

namespace amici {
namespace model_SPARCED {

void total_cl_SPARCED(realtype *total_cl, const realtype *x_rdata){
    total_cl[0] = EGFR;
}

} // namespace amici
} // namespace model_SPARCED