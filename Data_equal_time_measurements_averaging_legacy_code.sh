#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

dtau=0.05
run_no=25

for U in 10.00
do
   for Mu in $(seq -10.00 0.20 10.00) $(seq -20.00 0.50 20.00) 
   do
       for L in $(seq 10 2 42)
       do           
       
        python3 ./Data_equal_time_measurements_connected_correlators_neighbor_averaged_normal_average_legacy_code.py 10 $U $Mu $L $dtau $run_no &> connected_neighbor_averaged_normal_average_output.out

        echo "U_$U, Mu_$Mu, L_$L"
       
        done
   done
done



