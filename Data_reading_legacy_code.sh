#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

N=12
dtau=0.025
run_no=25

for U in 8.0 
do
      for Mu in $(seq 0.00 0.20 7.00) 
      do
             for L in 10 12 14 16 18 20 30 40 40 50 60 70 80 #20 30 40
             do
           
                   python3 ./Data_read_legacy_code_equal_time.py $N $U $Mu $L $dtau $run_no &> eq_time_data_read_output.out
                   python3 ./Data_read_legacy_code_unequal_time.py $N $U $Mu $L $dtau $run_no &> un_eq_time_data_read_output.out                    
             done
      echo "U_$U, Mu_$Mu, L_$L"
      done
      #echo "U_$U, Mu_$Mu, L_$L"
done



