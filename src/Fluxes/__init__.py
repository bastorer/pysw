# Fluxes
#     This module provides the flux options
#     for the simulator class.


# Import the flux options
#from flux2 import flux2
#from flux2_vanish import flux2_vanish
#from flux2_vanish_v2 import flux2_vanish_v2
#from flux_1d_1L import flux_1d_1L
#from advection import advection
from FV_SW import fv_sw_flux, fv_sw_source
#from FV_SW_vanish import fv_sw_flux_vanish, fv_sw_source_vanish
#from FV_SW_Zhang import fv_sw_flux_zhang, fv_sw_source_zhang
#from FV_SW_linear import fv_sw_linear_flux, fv_sw_linear_source
from FV_Burger import fv_Burger_flux, fv_Burger_source
from FV_advection import fv_advection_flux, fv_advection_source
from SPECTRAL_SW import spectral_sw_flux, spectral_sw_source
