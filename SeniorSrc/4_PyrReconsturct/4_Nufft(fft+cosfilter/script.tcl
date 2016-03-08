############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 2015 Xilinx Inc. All rights reserved.
############################################################
open_project S0_FFT1024
set_top s2_fft_1024_stream3
add_files S0_FFT1024/forwardfft1024.cpp
add_files S0_FFT1024/forwardfft1024.h
add_files nufft_cosfilter.h
add_files -tb S0_FFT1024/forwardtb.cpp
open_solution "solution3"
set_part {xc7z020clg484-1}
create_clock -period 10 -name default
config_dataflow -default_channel fifo -fifo_depth 1
config_schedule -effort high -verbose
#source "./S0_FFT1024/solution3/directives.tcl"
csim_design -O
csynth_design
cosim_design -O -trace_level all
export_design -evaluate verilog -format ip_catalog -vendor "pitchaya.com" -display_name "S2_FFT_1024"
