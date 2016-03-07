#include "median_filter.h"
#include <stdio.h>
#include <string.h>

int main() {
int width = 1520;
int hight = 10;
hls::stream<data_t> output("test");
hls::stream<data_t> src("test2");

int i ;

for(int j = 0;j<10 ;j++){
for (i = 0; i < 1520; i++) {

	data_t val;
	val = (t_output_scalar)temp[i];
	//val.last =0;
	src.write(val);
	printf("%d  %f \n",i,temp[i]);
}
}
median_strm(width,hight,src,output);

FILE * fo = fopen("mdf_out.txt", " wb");
	for( i=0;i<15200;i++)

	{
		data_t val = output.read();
		//data_t val = 0.789;
		//std::cout << val << std::endl;
		fprintf(fo, "%f  \n", val.to_double());
		//printf("%f \n",val.data.to_double());
		//printf("%d \n",i);
	}
	fclose(fo);
	return 0;
}

