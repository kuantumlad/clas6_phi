#include <vector>
#include <map>
#include <cmath>


bool hadronDeltaVzPass( int ExpOrSim, Float_t vx[], Float_t vy[], Float_t vz[], Float_t p[], Float_t cx[], Float_t cy[], Float_t cz[], UChar_t dc_sect[], int hadron_index, int e_index[] ){
  
  // using delta vertex cut to restrict hadrons produced in reaction close to reaction center and inside target.
  // cut is at +- 5.0 cm, the target positions of the rear and front windows.
  
  double hadron_corr_vz  = getCorrZ(ExpOrSim, vx[hadron_index], vy[hadron_index], vz[hadron_index], p[hadron_index]*cx[hadron_index], p[hadron_index]*cy[hadron_index], p[hadron_index]*cz[hadron_index], dc_sect[hadron_index]);
  double el_corr_vz = getCorrZ(ExpOrSim, vx[e_index[1]], vy[e_index[1]], vz[e_index[1]], p[e_index[1]]*cx[e_index[1]], p[e_index[1]]*cy[e_index[1]], p[e_index[1]]*cz[e_index[1]], dc_sect[e_index[1]]);
  double delta_vz = el_corr_vz - hadron_corr_vz;
  bool delta_vertex_cut = std::abs( el_corr_vz - hadron_corr_vz ) < 5.0; 
  
  if( delta_vertex_cut ) return true;
  
  return false;
}
