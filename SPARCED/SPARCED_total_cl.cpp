#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <array>

#include "x_rdata.h"

namespace amici {
namespace model_SPARCED {

void total_cl_SPARCED(realtype *total_cl, const realtype *x_rdata){
    total_cl[0] = Ce_Cdk2_pRBp_E2F;
    total_cl[1] = Ce_Cdk2_pRBp;
    total_cl[2] = Cd_Cdk46_pRB_E2F;
    total_cl[3] = Cd_Cdk46_pRB;
    total_cl[4] = pRBpp_E2F;
}

} // namespace amici
} // namespace model_SPARCED