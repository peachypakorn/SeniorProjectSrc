#include "../common_type.h"
#include "s0_nufftinput.h"
//D:\3DTV\Hardware\FPGA\AXIS\S0_pyrin\s0_nufftinput.h

#define C   3

template<typename T>
inline ap_fixed<16, 3> conv163( const T a) {
	ap_fixed<16,3> v;
	v.range(15,0) = a;
	return v;
}
template<typename T>
inline ap_fixed<16, 8> conv168( const T a) {
	ap_fixed<16,8> v;
	v.range(15,0) = a;
	return v;
}

void nufft_input_interface( hls::stream< packed_data >    &src_axi,
				 	   hls::stream<t_input_complex>        sigViews[4][3],
					   hls::stream<t_disp_scalar>          dispFilters[4]) {
	const t_disp_scalar  conV = 4.5;

	ap_fixed< 16,8> dA = 2.3;

	// Enable this if instantiated this alone
//#pragma HLS INTERFACE axis port=src_axi


#pragma HLS INTERFACE axis port=sigViews
#pragma HLS INTERFACE axis port=dispFilters

hls::stream<t_input_complex>  sigViewsInternal[4][3];
hls::stream<t_disp_scalar>    dispFiltersInternal[4];

//#pragma HLS data_pack variable=src_axi;
#pragma HLS data_pack variable=sigViews
#pragma HLS stream variable=sigViews depth=1

#pragma HLS data_pack variable=sigViewsInternal
#pragma HLS stream variable=sigViewsInternal depth=1520
#pragma HLS stream variable=dispFiltersInternal depth=1520


//#pragma HLS RESOURCE core=AXIS variable=src_axi   metadata="-bus_bundle INPUT_STREAM"
//#pragma HLS RESOURCE core=AXIS variable=sigViews  metadata="-bus_bundle OUTPUT_STREAM"
#pragma HLS DATAFLOW

	for(int i=0;i<1520;i++) {
#pragma HLS pipeline
		ap_uint<112>  packet = src_axi.read();
// Extracting PYR
				
		t_input_complex   pyrRGB[3];
#pragma HLS data_pack variable=pyrRGB
		/*
		pyrRGB[0].real().range(15,0) = conv163( packet.range(15,0)).range(15,0);
		pyrRGB[0].imag().range(15,0) = conv163( packet.range(31,16)).range(15,0);

		pyrRGB[1].real().range(15,0) = conv163( packet.range(47,32)).range(15,0);
		pyrRGB[1].imag().range(15,0) = conv163( packet.range(63 ,48)).range(15,0);

		pyrRGB[2].real().range(15,0) = conv163( packet.range(79,64)).range(15,0);
		pyrRGB[2].imag().range(15,0) = conv163( packet.range(95,80)).range(15,0);
*/
		ap_fixed<16,3> A[6];
		A[0].range(15,0) = packet.range(15,0);
		A[1].range(15,0) = packet.range(31,16);
		A[2].range(15,0) = packet.range(47,32);
		A[3].range(15,0) = packet.range(63 ,48);
		A[4].range(15,0) = packet.range(79,64);
		A[5].range(15,0) = packet.range(95,80);
		pyrRGB[0].real() = A[0];
		pyrRGB[0].imag()= A[1];

		pyrRGB[1].real()= A[2];
		pyrRGB[1].imag()= A[3];

		pyrRGB[2].real()= A[4];
		pyrRGB[2].imag()= A[5];


		// Extracting DISP
		ap_fixed<16,8> v;
		v.range(15,0) = packet.range( 111, 96);
		t_disp_scalar     depthScalarVal = v;

		for(int v=0;v<4;v++) {
#pragma HLS loop_unroll
			for(int c=0;c<C;c++) {
#pragma HLS loop_unroll
				sigViewsInternal[v][c].write(pyrRGB[c]);
			}
			t_disp_scalar displacement = ((conV - (v+1)) * dA - t_disp_scalar(0.5))*depthScalarVal;
			dispFiltersInternal[v].write(displacement);
		}
		//dispFilter.write(dispScalarVal);
	}
	
	for(int i=0;i<1520;i++) {
#pragma HLS pipeline
		for(int c=0;c<C;c++)
			sigViews[0][c].write(sigViewsInternal[0][c].read());
		dispFilters[0].write(dispFiltersInternal[0].read());
	}
	for(int i=0;i<1520;i++) {
#pragma HLS pipeline
		for(int c=0;c<C;c++)
			sigViews[1][c].write(sigViewsInternal[1][c].read());
		dispFilters[1].write(dispFiltersInternal[1].read());

	}

	for(int i=0;i<1520;i++) {
#pragma HLS pipeline
		for(int c=0;c<C;c++)
			sigViews[2][c].write(sigViewsInternal[2][c].read());
		dispFilters[2].write(dispFiltersInternal[2].read());

	}

	for(int i=0;i<1520;i++) {
#pragma HLS pipeline
		for(int c=0;c<C;c++)
			sigViews[3][c].write(sigViewsInternal[3][c].read());
		dispFilters[3].write(dispFiltersInternal[3].read());

	}

	
}
