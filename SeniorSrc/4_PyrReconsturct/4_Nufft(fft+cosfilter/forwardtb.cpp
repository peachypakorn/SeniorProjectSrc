
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include "forwardfft1024.h"
#include <stdio.h>

#include <fstream>
#include <string>
#include <sstream>

using namespace std;
// Temporary Declaration
typedef ap_fixed<17, 4> t_nufft_output_scalar;
typedef std::complex< t_nufft_output_scalar > t_nufft_output_complex;


#define BUF_SIZE 64
int main()
{
	static hls::stream< t_nufft_output_complex > xn_input[3];
		static hls::stream< t_nufft_output_complex> xk_output[3];
		int ndat = 1024;
		//t_nufft_output_complex testinput[ndat];

//		int ndat = 256 ;
		FILE * fi = fopen("inputMem.txt", " wb");
			for(int j=0;j<4;j++)
		{

			for(int i=0;i<ndat;i++) {

				t_nufft_output_complex val ((float)( rand()^0xFFFF) /65535, (float)( rand()^0xFFFF) /65535);
				if(j==0){
					fprintf(fi,"%.8f %.8f \n",val.real().to_float(),val.imag().to_float());
				//testinput[i] = val;
				}
				xn_input[0].write( val);
				xn_input[1].write( val);
				xn_input[2].write( val);
			}

			s2_fft_1024_stream3( xn_input, xk_output);
			FILE * fo = fopen("pyr_out.txt", " wb");

			for(int i=0;i<ndat/2;i++) {
				t_nufft_output_complex valOut1 = xk_output[0].read();
				t_nufft_output_complex valOut2 =xk_output[1].read();
				t_nufft_output_complex valOut3 =xk_output[2].read();
				fprintf(fo, "%.8f %.8f  ", valOut1.real().to_float(), valOut1.imag().to_float());
				fprintf(fo, "%.8f %.8f  ", valOut2.real().to_float(), valOut2.imag().to_float());
				fprintf(fo, "%.8f %.8f \n", valOut3.real().to_float(), valOut3.imag().to_float());
				std::cout << valOut1 << valOut2 << valOut3<< std::endl;
			}

			fclose(fo);
		}
			fclose(fi);
			return 0;
}
