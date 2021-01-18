REML_Rm <- function(invsqrtW_,JYp_,loopsize_, patternlength_,rownumber_,ybetas_,etahat_,
tempmatR_,JXpi_,JXpp_,JXpx_,JXpdim_,JZpi_,JZpp_,JZpx_,JZpdim_,betacov_,C12_){
.Call( "REML_Rm_cpp",invsqrtW_,JYp_,loopsize_, patternlength_,rownumber_,ybetas_,etahat_,
tempmatR_,JXpi_,JXpp_,JXpx_,JXpdim_,JZpi_,JZpp_,JZpx_,JZpdim_,betacov_,C12_)
}
