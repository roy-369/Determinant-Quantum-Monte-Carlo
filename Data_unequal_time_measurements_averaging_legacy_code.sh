#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

N=12
Dtau=0.025
Run_no=20

for U in 8.0
do
      for Mu in $(seq 0.00 0.20 7.00) 
      do
             for L in 10 12 14 16 18 20 30 40 50 60 70 80
             do       
                 python3 ./Data_unequal_time_measurements_average_legacy_code.py $N $U $Mu $L $Dtau $Run_no &> un_eq_average_output.out                  
             done
      echo "U_$U, Mu_$Mu, L_$L"
      done

     # echo "U_$U, Mu_$Mu, L_$L"
done



