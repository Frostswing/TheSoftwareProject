#!/bin/bash

mkdir val_mem_leeks
valgrind ./symnmf sym inputs/spheres.txt > val_mem_leeks/symc_spheres_res.txt
echo
echo
valgrind ./symnmf ddg inputs/spheres.txt > val_mem_leeks/ddgc_spheres_res.txt
echo
echo
valgrind ./symnmf norm inputs/spheres.txt > val_mem_leeks/normc_spheres_res.txt
echo
echo
valgrind ./symnmf sym inputs/iris.txt > val_mem_leeks/symc_iris_res.txt
echo
echo
valgrind ./symnmf ddg inputs/iris.txt > val_mem_leeks/ddgc_iris_res.txt
echo
echo
valgrind ./symnmf norm inputs/iris.txt > val_mem_leeks/normc_iris_res.txt
echo
echo
valgrind ./symnmf sym inputs/drive1.txt > val_mem_leeks/symc_drive1_res.txt
echo
echo
valgrind ./symnmf ddg inputs/drive1.txt > val_mem_leeks/ddgc_drive1_res.txt
echo
echo
valgrind ./symnmf norm inputs/drive1.txt > val_mem_leeks/normc_drive1_res.txt
echo
echo
valgrind ./symnmf sym inputs/drive2.txt > val_mem_leeks/symc_drive2_res.txt
echo
echo
valgrind ./symnmf ddg inputs/drive2.txt > val_mem_leeks/ddgc_drive2_res.txt
echo
echo
valgrind ./symnmf norm inputs/drive2.txt > val_mem_leeks/normc_drive2_res.txt
echo
echo

rm -r val_mem_leeks

