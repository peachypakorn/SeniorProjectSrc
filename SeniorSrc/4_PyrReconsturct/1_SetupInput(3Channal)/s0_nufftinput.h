#ifndef __S0__NUFFT__
#define __S0__NUFFT__
#include "../common_type.h"



typedef ap_axiu<112, 1, 1, 1> PACKETIN;
typedef ap_axiu<24, 1, 1, 1> PIXEL;

typedef ap_uint<112> packed_data;

void nufft_input_interface( hls::stream< packed_data >    &src_axi,
				 	   hls::stream<t_input_complex> sigViews[4][3],
					   hls::stream<t_disp_scalar>          dispFilters[4]);

#endif
