#!/usr/bin/sh

paste x_coord.dat y_coord.dat pmf_1d.dat > pmf_.dat
sort -g pmf_.dat > pmf_xy.dat

paste x_coord.dat z_coord.dat pmf_1d.dat > pmf_.dat
sort -g pmf_.dat > pmf_xz.dat

paste y_coord.dat z_coord.dat pmf_1d.dat > pmf_.dat
sort -g pmf_.dat > pmf_yz.dat

paste x_coord.dat y_coord.dat z_coord.dat pmf_1d.dat > pmf_.dat
sort -g pmf_.dat > pmf_xyz.dat
