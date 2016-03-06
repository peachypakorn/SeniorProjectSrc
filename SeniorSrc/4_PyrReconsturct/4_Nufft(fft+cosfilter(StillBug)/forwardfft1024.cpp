#include <hls_stream.h>
#include "forwardfft1024.h"

#include "ap_fixed.h"
#include "hls_fft.h"
#include "../nufft_cosfilter.h"

// configurable params
const char FFT_INPUT_WIDTH                     = 16;
const char FFT_OUTPUT_WIDTH                    = 32;
const char FFT_CONFIG_WIDTH                    = 8;
const char FFT_NFFT_MAX                        = 10;
const int  FFT_LENGTH                          = 1 << FFT_NFFT_MAX;


#include <complex>
using namespace std;


struct config1 : hls::ip_fft::params_t {
    static const unsigned ordering_opt = hls::ip_fft::natural_order;

    static const unsigned input_width = 16;
    static const unsigned output_width = 32;

    static const unsigned scaling_opt = hls::ip_fft::unscaled;

    static const unsigned config_width = FFT_CONFIG_WIDTH;
    static const unsigned max_nfft = 10;
        static const bool has_nfft = false;
};

typedef hls::ip_fft::config_t<config1> config_t;
typedef hls::ip_fft::status_t<config1> status_t;

typedef ap_fixed<FFT_INPUT_WIDTH,1> data_in_t;
typedef ap_fixed<FFT_OUTPUT_WIDTH,FFT_OUTPUT_WIDTH-FFT_INPUT_WIDTH+1> data_out_t;
typedef std::complex<data_in_t> cmpxDataIn;
typedef std::complex<data_out_t> cmpxDataOut;

void dummy_proc_fe(
    bool direction,
    config_t* config,
    cmpxDataIn in[FFT_LENGTH],
    cmpxDataIn out[FFT_LENGTH]);

void dummy_proc_be(
    status_t* status_in,
    bool* ovflo,
    cmpxDataOut in[FFT_LENGTH],
    cmpxDataOut out[FFT_LENGTH]);

void fft_top(
    bool direction,
    cmpxDataIn in[FFT_LENGTH],
    cmpxDataOut out[FFT_LENGTH],
    bool* ovflo);


void dummy_proc_fe(
    bool direction,
    config_t* config,
    cmpxDataIn in[FFT_LENGTH],
    cmpxDataIn out[FFT_LENGTH])
{
    int i;
    config->setDir(direction);
//    config->setSch(0x2AB);
    for (i=0; i< FFT_LENGTH; i++)
        out[i] = in[i];
}

void dummy_proc_be(
    status_t* status_in,
    bool* ovflo,
    cmpxDataOut in[FFT_LENGTH],
    cmpxDataOut out[FFT_LENGTH])
{
    int i;
    for (i=0; i< FFT_LENGTH; i++)
        out[i] = in[i];
//    *ovflo = status_in->getOvflo() & 0x1;
}

template < int unrelated_manner>
void fft_top(
    bool direction,
    complex<data_in_t> in[FFT_LENGTH],
    complex<data_out_t> out[FFT_LENGTH],
    bool* ovflo)
{
#pragma HLS interface ap_hs port=direction
#pragma HLS interface ap_fifo depth=1 port=ovflo
#pragma HLS interface ap_fifo depth=FFT_LENGTH port=in,out
#pragma HLS data_pack variable=in
#pragma HLS data_pack variable=out
#pragma HLS dataflow
    complex<data_in_t> xn[FFT_LENGTH];
    complex<data_out_t> xk[FFT_LENGTH];
    config_t fft_config;
    status_t fft_status;

    dummy_proc_fe(direction, &fft_config, in, xn);
    // FFT IP
    hls::fft<config1>(xn, xk, &fft_status, &fft_config);
    dummy_proc_be(&fft_status, ovflo, xk, out);
}


// Temporary Declaration
typedef ap_fixed<17, 4> t_nufft_output_scalar;
typedef std::complex< t_nufft_output_scalar > t_nufft_output_complex;






template < int unrelated_manner>
void fft_top_stream( hls::stream< t_nufft_output_complex > &in,
	                 hls::stream< t_nufft_output_complex > &out) {
#pragma HLS data_pack variable=in
#pragma HLS data_pack variable=out

#pragma HLS DATAFLOW
	complex<data_in_t> inM[FFT_LENGTH];
	complex<data_out_t> outM[FFT_LENGTH];
#pragma HLS data_pack variable=inM
#pragma HLS data_pack variable=outM
	for(int i=0;i<FFT_LENGTH;i++) {
#pragma HLS PIPELINE
		t_nufft_output_complex valIn = in.read();
		inM[i] = complex<data_in_t>(valIn.real(), valIn.imag());
	}
	bool direction;
	bool ovflo;
	fft_top<unrelated_manner>(true, inM, outM, &ovflo);

#ifndef __SYNTHESIS__
	{
		FILE * fo = fopen("fftOut.txt", "wb");
		for(int i=0;i<512;i++) {
			fprintf(fo, "%.8f %.8f\n", outM[i].real().to_float(), outM[i].imag().to_float());
		}
		fclose(fo);
	}
#endif

	for(int i=0;i<FFT_LENGTH;i++) {
#pragma HLS PIPELINE
		complex<data_out_t> valOut = outM[i];
		t_nufft_output_complex val(valOut.real(), valOut.imag());
		out.write(val);
	}
}

void s2_fft_1024_stream3(
		hls::stream< t_nufft_output_complex > in[3],
		hls::stream< t_nufft_output_complex > out[3]) {
#ifndef  NUFFTB
	#pragma HLS INTERFACE axis port=in
	#pragma HLS INTERFACE axis port=out
#endif

#pragma HLS data_pack variable=in
#pragma HLS data_pack variable=out

#pragma HLS DATAFLOW

	hls::stream< t_nufft_output_complex >  fftOut[3];
#pragma HLS data_pack variable=fftOut

	fft_top_stream<1>(in[0], fftOut[0]);
	fft_top_stream<2>(in[1], fftOut[1]);
	fft_top_stream<3>(in[2], fftOut[2]);

	//nufft_cosfilter<3, t_nufft_output_complex>(fftOut, out, 512, 2);
	nufft_cosfilter1024<3, t_nufft_output_complex>(fftOut, out);

}
