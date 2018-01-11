#!/bin/bash
FOLDER=/home/dlecina/daniel/scripts/msm
FOLDER=.
rm discretized/*
rm matrix/*
rm pcca/*

$FOLDER/mm_1_3_cluster
$FOLDER/mm_1_3_assign
$FOLDER/mm_1_3_connectivity
$FOLDER/mm_1_3_timescales
