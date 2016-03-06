#ifndef __FORWARDFFT512__
#define __FORWARDFFT512__

#include "ap_fixed.h"
#include "hls_fft.h"

#include "../common_type.h"

#include <complex>


// Temporary Declaration
//typedef ap_fixed<17, 4> t_nufft_output_scalar;
//typedef std::complex< t_nufft_output_scalar > t_nufft_output_complex;


void s2_fft_512_stream3(
		hls::stream< t_nufft_output_complex > &in0,
		hls::stream< t_nufft_output_complex > &in1,
		hls::stream< t_nufft_output_complex > &in2,
		hls::stream< t_nufft_output_complex > &out0,
		hls::stream< t_nufft_output_complex > &out1,
		hls::stream< t_nufft_output_complex > &out2) ;


#endif
