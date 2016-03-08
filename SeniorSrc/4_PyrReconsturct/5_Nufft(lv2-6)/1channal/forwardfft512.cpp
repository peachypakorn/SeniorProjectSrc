#include <hls_stream.h>
#include "forwardfft512.h"
#include "../nufft_cosfilter.h"
#include "hls_dsp.h"
//#include "hls_cmpy.h"

#include "ap_fixed.h"
#include "hls_fft.h"
//Configuration parameters for this instance.
//typedef hls::CmpyThreeMult ARCH;

// configurable params
const char FFT_INPUT_WIDTH                     = 16;
const char FFT_OUTPUT_WIDTH                    = 32;
const char FFT_CONFIG_WIDTH                    = 8;



#include <complex>
using namespace std;


struct config9 : hls::ip_fft::params_t {
    static const unsigned ordering_opt = hls::ip_fft::natural_order;
    static const unsigned config_width = 8;

    static const unsigned input_width = 16;
    static const unsigned output_width = 32;

    static const unsigned scaling_opt = hls::ip_fft::unscaled;

    static const unsigned max_nfft = 9;
    static const bool has_nfft = false;
    static const unsigned  channels = 1;

    static const unsigned arch_opt = hls::ip_fft::pipelined_streaming_io;
};


struct config8 : hls::ip_fft::params_t {
    static const unsigned ordering_opt = hls::ip_fft::natural_order;
    static const unsigned config_width = 8;

    static const unsigned input_width = 16;
    static const unsigned output_width = 32;

    static const unsigned scaling_opt = hls::ip_fft::unscaled;

    static const unsigned max_nfft = 8;
    static const bool has_nfft = false;
    static const unsigned  channels = 1;

    static const unsigned arch_opt = hls::ip_fft::pipelined_streaming_io;
};

struct config7 : hls::ip_fft::params_t {
    static const unsigned ordering_opt = hls::ip_fft::natural_order;
    static const unsigned config_width = 8;

    static const unsigned input_width = 16;
    static const unsigned output_width = 24;

    static const unsigned scaling_opt = hls::ip_fft::unscaled;

    static const unsigned max_nfft = 7;
    static const bool has_nfft = false;
    static const unsigned  channels = 1;

    static const unsigned arch_opt = hls::ip_fft::pipelined_streaming_io;
};
struct config6 : hls::ip_fft::params_t {
    static const unsigned ordering_opt = hls::ip_fft::natural_order;
    static const unsigned config_width = 8;

    static const unsigned input_width = 16;
    static const unsigned output_width = 24;

    static const unsigned scaling_opt = hls::ip_fft::unscaled;

    static const unsigned max_nfft = 6;
    static const bool has_nfft = false;
    static const unsigned  channels = 1;

    static const unsigned arch_opt = hls::ip_fft::pipelined_streaming_io;
};

//typedef hls::ip_fft::config_t<config9> config_t;
//typedef hls::ip_fft::status_t<config9> status_t;

typedef ap_fixed<FFT_INPUT_WIDTH,1> data_in_t;
typedef ap_fixed<FFT_OUTPUT_WIDTH,FFT_OUTPUT_WIDTH-FFT_INPUT_WIDTH+1> data_out_t;
typedef std::complex<data_in_t> cmpxDataIn;
typedef std::complex<data_out_t> cmpxDataOut;



template< typename config_f, typename Tin, int FFT_LENGTH>
void dummy_proc_fe(
    bool direction,
	hls::ip_fft::config_t<config_f> * config,
    Tin in[FFT_LENGTH],
    Tin out[FFT_LENGTH])
{
    int i;
    config->setDir(direction);
    if(FFT_LENGTH == 256)config->setNfft(8);
    else if(FFT_LENGTH == 256)config->setNfft(7);
    else if(FFT_LENGTH == 256)config->setNfft(6);
    else config->setNfft(9);
   // config->setNfft(log_size);
//    config->setSch(0x2AB);
    for (i=0; i< FFT_LENGTH; i++)
        out[i] = in[i];
}

template< typename config_f, typename Tout, int FFT_LENGTH>
void dummy_proc_be(
		hls::ip_fft::status_t<config_f> * status_in,
    bool* ovflo,
	Tout in[FFT_LENGTH],
	Tout out[FFT_LENGTH])
{
    int i;
    for (i=0; i< FFT_LENGTH; i++)
        out[i] = in[i];
//    *ovflo = status_in->getOvflo() & 0x1;
}

template < typename config_f, int FFT_LENGTH, typename Tin, typename Tout, int unrelated_manner>
void fft_top(
    bool direction,
    Tin in[FFT_LENGTH],
    Tout out[FFT_LENGTH],
    bool* ovflo)
{
//#pragma HLS inline

#pragma HLS interface ap_hs port=direction
#pragma HLS interface ap_fifo depth=1 port=ovflo
#pragma HLS interface ap_fifo depth=FFT_LENGTH port=in,out
#pragma HLS data_pack variable=in
#pragma HLS data_pack variable=out
#pragma HLS dataflow
    Tin  xn[FFT_LENGTH];
    Tout xk[FFT_LENGTH];

#pragma HLS data_pack variable=xn
#pragma HLS data_pack variable=xk
//#pragma HLS STREAM  variable=xn depth=512
//#pragma HLS STREAM  variable=xk depth=512
    hls::ip_fft::config_t<config_f> fft_config;
    hls::ip_fft::status_t<config_f> fft_status;

#pragma HLS STREAM  variable=fft_config depth=1
    dummy_proc_fe<config_f, Tin, FFT_LENGTH>(direction, &fft_config, in, xn);
    // FFT IP
    hls::fft<config_f>(xn, xk, &fft_status, &fft_config);
    dummy_proc_be<config_f, Tout, FFT_LENGTH>(&fft_status, ovflo, xk, out);
}


// Temporary Declaration
//typedef ap_fixed<17, 4> t_nufft_output_scalar;
//typedef std::complex< t_nufft_output_scalar > t_nufft_output_complex;

void nufft_cosfilter1(hls::stream< t_nufft_output_complex > nufftIn[1],
					 hls::stream< t_nufft_output_complex > nufftOut[1],
					 const int nL, const int m) {

#pragma HLS data_pack variable=nufftIn
#pragma HLS data_pack variable=nufftOut

#pragma HLS DATAFLOW
	nufft_cosfilter<1>(nufftIn, nufftOut, nL, m);
}


void nufft_cosfilter3(hls::stream< t_nufft_output_complex > nufftIn[3],
					 hls::stream< t_nufft_output_complex > nufftOut[3],
					 const int nL, const int m) {

#pragma HLS data_pack variable=nufftIn
#pragma HLS data_pack variable=nufftOut
#pragma HLS DATAFLOW
	nufft_cosfilter<3>(nufftIn, nufftOut, nL, m);
}


template < typename config_f, int DEPTHIN, int FFT_LENGTH, typename Tin, typename Tout, int unrelated_manner>
void fft_top_stream( hls::stream< t_nufft_output_complex > &in,
	                 hls::stream< t_nufft_output_complex > &out) {
#pragma HLS inline
#pragma HLS data_pack variable=in
#pragma HLS data_pack variable=out

#pragma HLS DATAFLOW
	Tin  inM[FFT_LENGTH];
	Tout outM[FFT_LENGTH];
#pragma HLS data_pack variable=inM
#pragma HLS data_pack variable=outM
#pragma HLS STREAM variable=DEPTHIN
	for(int i=0;i<FFT_LENGTH;i++) {
#pragma HLS PIPELINE
		t_nufft_output_complex valIn = in.read();
		inM[i] = complex<data_in_t>(valIn.real(), valIn.imag());
	}
	bool direction;
	bool ovflo;
	fft_top<config_f, FFT_LENGTH, Tin, Tout, unrelated_manner>(true, inM, outM, &ovflo);
	for(int i=0;i<FFT_LENGTH;i++) {
#pragma HLS PIPELINE
		Tout valOut = outM[i];
		t_nufft_output_complex val(valOut.real(), valOut.imag());
		out.write(val);
	}
#ifndef __SYNTHESIS__
	{
		if(FFT_LENGTH==512){
		FILE * fo = fopen("fftOut.txt", "wb");
		for(int i=0;i<512;i++) {
			fprintf(fo, "%.8f %.8f\n", outM[i].real().to_double(), outM[i].imag().to_float());
		}
		fclose(fo);
		}
	}
#endif
}

template < typename T, int C>
void splitStream(hls::stream<T> &strIn,
			     hls::stream<T> strOut[C], int length) {
#pragma HLS inline
#pragma HLS DATAFLOW

	for(int i=0;i<length*C;i++) {
#pragma HLS PIPELINE
		T v = strIn.read();
		if (i>=0 && i<length) {
			strOut[0].write(v);
		}
		if (i>=length && i<length*2) {
			strOut[1].write(v);
		}
		if (i>=(length*2) && i<length*3) {
			strOut[2].write(v);
		}
	}

}

void s2_fft_512_stream3(
		hls::stream< t_nufft_output_complex > &in0,
		hls::stream< t_nufft_output_complex > &in1,
		hls::stream< t_nufft_output_complex > &in2,
		hls::stream< t_nufft_output_complex > &out0,
		hls::stream< t_nufft_output_complex > &out1,
		hls::stream< t_nufft_output_complex > &out2) {
#pragma HLS inline off

#pragma HLS INTERFACE axis port=in0
#pragma HLS INTERFACE axis port=in1
#pragma HLS INTERFACE axis port=in2
#pragma HLS INTERFACE axis port=out0
#pragma HLS INTERFACE axis port=out1
#pragma HLS INTERFACE axis port=out2

#pragma HLS data_pack variable=in0
#pragma HLS data_pack variable=out0
#pragma HLS data_pack variable=in1
#pragma HLS data_pack variable=out1
#pragma HLS data_pack variable=in2
#pragma HLS data_pack variable=out2

#pragma HLS DATAFLOW

	// Assuming that the buffer is actually on the other side
	hls::stream< t_nufft_output_complex > FFT512InSplit[3];
	hls::stream< t_nufft_output_complex > FFT256InSplit[3];
	hls::stream< t_nufft_output_complex > FFT128InSplit[3];
	hls::stream< t_nufft_output_complex > FFT64InSplit[3];

#pragma HLS STREAM variable=FFT512InSplit depth=512
#pragma HLS STREAM variable=FFT256InSplit depth=256
#pragma HLS STREAM variable=FFT128InSplit depth=128
#pragma HLS STREAM variable=FFT64InSplit depth=128

#pragma HLS DATA_PACK variable=FFT512InSplit
#pragma HLS DATA_PACK variable=FFT256InSplit
#pragma HLS DATA_PACK variable=FFT128InSplit
#pragma HLS DATA_PACK variable=FFT64InSplit


	hls::stream< t_nufft_output_complex > FFT512Out;
	hls::stream< t_nufft_output_complex > FFT256Out;
	hls::stream< t_nufft_output_complex > FFT128Out;
	hls::stream< t_nufft_output_complex > FFT64Out;
#pragma HLS DATA_PACK variable=FFT512Out
#pragma HLS DATA_PACK variable=FFT256Out
#pragma HLS DATA_PACK variable=FFT128Out
#pragma HLS DATA_PACK variable=FFT64Out

	hls::stream< t_nufft_output_complex > nFFT32Out;
	hls::stream< t_nufft_output_complex > nFFT256Out;
	hls::stream< t_nufft_output_complex > nFFT128Out;
	hls::stream< t_nufft_output_complex > nFFT64Out;
	#pragma HLS DATA_PACK variable=nFFT32Out
	#pragma HLS DATA_PACK variable=nFFT256Out
	#pragma HLS DATA_PACK variable=nFFT128Out
	#pragma HLS DATA_PACK variable=nFFT64Out

#pragma HLS STREAM variable=nFFT256Out depth=768
#pragma HLS STREAM variable=nFFT128Out depth=384
#pragma HLS STREAM variable=nFFT64Out depth=192
#pragma HLS STREAM variable=nFFT32Out depth=96


	hls::stream< t_nufft_output_complex > nFFT32InSplit[3];
	hls::stream< t_nufft_output_complex > nFFT256InSplit[3];
	hls::stream< t_nufft_output_complex > nFFT128InSplit[3];
	hls::stream< t_nufft_output_complex > nFFT64InSplit[3];
#pragma HLS DATA_PACK variable=nFFT32InSplit
#pragma HLS DATA_PACK variable=nFFT256InSplit
#pragma HLS DATA_PACK variable=nFFT128InSplit
#pragma HLS DATA_PACK variable=nFFT64InSplit


#pragma HLS STREAM variable=nFFT256InSplit depth=256
#pragma HLS STREAM variable=nFFT128InSplit depth=128
#pragma HLS STREAM variable=nFFT64InSplit depth=64
#pragma HLS STREAM variable=nFFT32InSplit depth=32


	LoopSplittingData:
	for(int i=0;i<512+256+128+64;i++) {
#pragma HLS PIPELINE rewind
		t_nufft_output_complex vin0 = in0.read();
		t_nufft_output_complex vin1 = in1.read();
		t_nufft_output_complex vin2 = in2.read();

		if (i<512) {
			FFT512InSplit[0].write(vin0);
			FFT512InSplit[1].write(vin1);
			FFT512InSplit[2].write(vin2);
		}
		if (i>=512 && i < 768) {
			FFT256InSplit[0].write(vin0);
			FFT256InSplit[1].write(vin1);
			FFT256InSplit[2].write(vin2);
		}
		if (i>=768 && i< 896) {
			FFT128InSplit[0].write(vin0);
			FFT128InSplit[1].write(vin1);
			FFT128InSplit[2].write(vin2);
		}
		if (i>=896&& i< 960) {
			FFT64InSplit[0].write(vin0);
			FFT64InSplit[1].write(vin1);
			FFT64InSplit[2].write(vin2);
		}
	}

	for(int c=0;c<3;c++) {
#pragma HLS DATAFLOW
		fft_top_stream< config9, 1536, 512, cmpxDataIn, cmpxDataOut, 0> ( FFT512InSplit[c], FFT512Out);
		nufft_cosfilter_single<t_nufft_output_complex>(FFT512Out, nFFT256Out, 256, 2);
	}
	for(int c=0;c<3;c++) {
#pragma HLS DATAFLOW
		fft_top_stream< config8, 768, 256, cmpxDataIn, cmpxDataOut, 0>( FFT256InSplit[c], FFT256Out);
		nufft_cosfilter_single<t_nufft_output_complex>(FFT256Out, nFFT128Out, 128, 2);
	}

	typedef std::complex< ap_fixed<24, (24-16+1) > >   config7type;

	//template < typename config_f, int DEPTHIN, int FFT_LENGTH, typename Tin, typename Tout, int unrelated_manner>
	for(int c=0;c<3;c++) {
#pragma HLS DATAFLOW
		fft_top_stream< config7, 384, 128, cmpxDataIn, config7type, 0> ( FFT128InSplit[c], FFT128Out);
		nufft_cosfilter_single<t_nufft_output_complex>(FFT128Out, nFFT64Out, 64, 2);
	}
	for(int c=0;c<3;c++) {
#pragma HLS DATAFLOW
		fft_top_stream< config6, 192, 64, cmpxDataIn, config7type, 0> ( FFT64InSplit[c], FFT64Out);
		nufft_cosfilter_single<t_nufft_output_complex>(FFT64Out, nFFT32Out, 32, 2);
	}

	splitStream<t_nufft_output_complex,3>(nFFT256Out,nFFT256InSplit, 256 );
	splitStream<t_nufft_output_complex,3>(nFFT128Out,nFFT128InSplit, 128 );
	splitStream<t_nufft_output_complex,3>(nFFT64Out ,nFFT64InSplit, 64 );
	splitStream<t_nufft_output_complex,3>(nFFT32Out ,nFFT32InSplit, 32 );
	// Test Data Out
	LoopDataOut:
	for(int i=0;i<(256+128+64+32);i++) {
#pragma HLS pipeline
		t_nufft_output_complex  v0;
		t_nufft_output_complex  v1;
		t_nufft_output_complex  v2;
		if (i<256) {
			v0 = nFFT256InSplit[0].read();
			v1 = nFFT256InSplit[1].read();
			v2 = nFFT256InSplit[2].read();
		}
		if (i>=256 && i < 384) {
			v0 = nFFT128InSplit[0].read();
			v1 = nFFT128InSplit[1].read();
			v2 = nFFT128InSplit[2].read();
		}
		if (i>=384 && i< 448) {
			v0 = nFFT64InSplit[0].read();
			v1 = nFFT64InSplit[1].read();
			v2 = nFFT64InSplit[2].read();
		}
		if (i>=448&& i< 480) {
			v0 = nFFT32InSplit[0].read();
			v1 = nFFT32InSplit[1].read();
			v2 = nFFT32InSplit[2].read();
		}
		out0.write(v0);
		out1.write(v1);
		out2.write(v2);
	}




}

void s2_fft_512_stream(
		hls::stream< t_nufft_output_complex > &in0,

		hls::stream< t_nufft_output_complex > &out0) {
#pragma HLS inline off

#pragma HLS INTERFACE axis port=in0

#pragma HLS INTERFACE axis port=out0


#pragma HLS data_pack variable=in0
#pragma HLS data_pack variable=out0


#pragma HLS DATAFLOW

	// Assuming that the buffer is actually on the other side
	hls::stream< t_nufft_output_complex > FFT512InSplit;
	hls::stream< t_nufft_output_complex > FFT256InSplit;
	hls::stream< t_nufft_output_complex > FFT128InSplit;
	hls::stream< t_nufft_output_complex > FFT64InSplit;

#pragma HLS STREAM variable=FFT512InSplit depth=512
#pragma HLS STREAM variable=FFT256InSplit depth=256
#pragma HLS STREAM variable=FFT128InSplit depth=128
#pragma HLS STREAM variable=FFT64InSplit depth=128

#pragma HLS DATA_PACK variable=FFT512InSplit
#pragma HLS DATA_PACK variable=FFT256InSplit
#pragma HLS DATA_PACK variable=FFT128InSplit
#pragma HLS DATA_PACK variable=FFT64InSplit


	hls::stream< t_nufft_output_complex > FFT512Out;
	hls::stream< t_nufft_output_complex > FFT256Out;
	hls::stream< t_nufft_output_complex > FFT128Out;
	hls::stream< t_nufft_output_complex > FFT64Out;
#pragma HLS DATA_PACK variable=FFT512Out
#pragma HLS DATA_PACK variable=FFT256Out
#pragma HLS DATA_PACK variable=FFT128Out
#pragma HLS DATA_PACK variable=FFT64Out

	hls::stream< t_nufft_output_complex > nFFT32Out;
	hls::stream< t_nufft_output_complex > nFFT256Out;
	hls::stream< t_nufft_output_complex > nFFT128Out;
	hls::stream< t_nufft_output_complex > nFFT64Out;
	#pragma HLS DATA_PACK variable=nFFT32Out
	#pragma HLS DATA_PACK variable=nFFT256Out
	#pragma HLS DATA_PACK variable=nFFT128Out
	#pragma HLS DATA_PACK variable=nFFT64Out

//#pragma HLS STREAM variable=nFFT256Out depth=768
//#pragma HLS STREAM variable=nFFT128Out depth=384
//#pragma HLS STREAM variable=nFFT64Out depth=192
//#pragma HLS STREAM variable=nFFT32Out depth=96

	#pragma HLS STREAM variable=nFFT256Out depth=256
	#pragma HLS STREAM variable=nFFT128Out depth=128
	#pragma HLS STREAM variable=nFFT64Out depth=64
	#pragma HLS STREAM variable=nFFT32Out depth=32

//	hls::stream< t_nufft_output_complex > nFFT32InSplit;
//	hls::stream< t_nufft_output_complex > nFFT256InSplit;
//	hls::stream< t_nufft_output_complex > nFFT128InSplit;
//	hls::stream< t_nufft_output_complex > nFFT64InSplit;
//#pragma HLS DATA_PACK variable=nFFT32InSplit
//#pragma HLS DATA_PACK variable=nFFT256InSplit
//#pragma HLS DATA_PACK variable=nFFT128InSplit
//#pragma HLS DATA_PACK variable=nFFT64InSplit

//
//#pragma HLS STREAM variable=nFFT256InSplit depth=256
//#pragma HLS STREAM variable=nFFT128InSplit depth=128
//#pragma HLS STREAM variable=nFFT64InSplit depth=64
//#pragma HLS STREAM variable=nFFT32InSplit depth=32


	LoopSplittingData:
	for(int i=0;i<512+256+128+64;i++) {
#pragma HLS PIPELINE rewind
		t_nufft_output_complex vin0 = in0.read();
		//t_nufft_output_complex vin1 = in1.read();
		//t_nufft_output_complex vin2 = in2.read();

		if (i<512) {
			FFT512InSplit.write(vin0);
			//FFT512InSplit[1].write(vin1);
			//FFT512InSplit[2].write(vin2);
		}
		if (i>=512 && i < 768) {
			FFT256InSplit.write(vin0);
			//FFT256InSplit[1].write(vin1);
			//FFT256InSplit[2].write(vin2);
		}
		if (i>=768 && i< 896) {
			FFT128InSplit.write(vin0);
			//FFT128InSplit[1].write(vin1);
			//FFT128InSplit[2].write(vin2);
		}
		if (i>=896&& i< 960) {
			FFT64InSplit.write(vin0);
			//FFT64InSplit[1].write(vin1);
			//FFT64InSplit[2].write(vin2);
		}
	}


		fft_top_stream< config9, 1536, 512, cmpxDataIn, cmpxDataOut, 0> ( FFT512InSplit, FFT512Out);
		nufft_cosfilter_single<t_nufft_output_complex>(FFT512Out, nFFT256Out, 256, 2);


		fft_top_stream< config9, 768, 256, cmpxDataIn, cmpxDataOut, 0>( FFT256InSplit, FFT256Out);
		nufft_cosfilter_single<t_nufft_output_complex>(FFT256Out, nFFT128Out, 128, 2);


	typedef std::complex< ap_fixed<24, (24-16+1) > >   config7type;

	//template < typename config_f, int DEPTHIN, int FFT_LENGTH, typename Tin, typename Tout, int unrelated_manner>

		fft_top_stream< config9, 384, 128, cmpxDataIn, cmpxDataOut, 0> ( FFT128InSplit, FFT128Out);
		nufft_cosfilter_single<t_nufft_output_complex>(FFT128Out, nFFT64Out, 64, 2);


		fft_top_stream< config9, 192, 64, cmpxDataIn, cmpxDataOut, 0> ( FFT64InSplit, FFT64Out);
		nufft_cosfilter_single<t_nufft_output_complex>(FFT64Out, nFFT32Out, 32, 2);


//	splitStream<t_nufft_output_complex,3>(nFFT256Out,nFFT256InSplit, 256 );
//	splitStream<t_nufft_output_complex,3>(nFFT128Out,nFFT128InSplit, 128 );
//	splitStream<t_nufft_output_complex,3>(nFFT64Out ,nFFT64InSplit, 64 );
//	splitStream<t_nufft_output_complex,3>(nFFT32Out ,nFFT32InSplit, 32 );
	// Test Data Out
	LoopDataOut:
	for(int i=0;i<(256+128+64+32);i++) {
#pragma HLS pipeline
		t_nufft_output_complex  v0;
		//t_nufft_output_complex  v1;
		//t_nufft_output_complex  v2;
		if (i<256) {
			v0 = nFFT256Out.read();
			//v1 = nFFT256InSplit[1].read();
			//v2 = nFFT256InSplit[2].read();
		}
		if (i>=256 && i < 384) {
			v0 = nFFT128Out.read();
			//v1 = nFFT128InSplit[1].read();
			//v2 = nFFT128InSplit[2].read();
		}
		if (i>=384 && i< 448) {
			v0 = nFFT64Out.read();
			//v1 = nFFT64InSplit[1].read();
			//v2 = nFFT64InSplit[2].read();
		}
		if (i>=448&& i< 480) {
			v0 = nFFT32Out.read();
			//v1 = nFFT32InSplit[1].read();
			//v2 = nFFT32InSplit[2].read();
		}
		out0.write(v0);
		//out1.write(v1);
		//out2.write(v2);
	}

}
