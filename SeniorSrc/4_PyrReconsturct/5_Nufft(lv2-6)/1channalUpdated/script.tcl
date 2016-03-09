############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 2015 Xilinx Inc. All rights reserved.
############################################################
open_project nufft512
set_top s2_fft_512_stream
add_files nufft512/forwardfft512.cpp
add_files nufft512/forwardfft512.h
add_files nufft_cosfilter.h
add_files nufft512/temp.txt
add_files -tb nufft512/forwardfft512_tb.cpp
open_solution "solution2"
set_part {xc7z020clg484-1}
create_clock -period 10 -name default
config_rtl -encoding onehot -reset all -reset_level high
config_dataflow -default_channel fifo -fifo_depth 1
#source "./nufft512/solution2/directives.tcl"
csim_design -O
csynth_design
cosim_design -O -trace_level all
export_design -format ip_catalog
