#!/bin/bash

# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

module load gnu/9.1.0


for U in 14.0
do
   mkdir -p U_$U 
   cd U_$U || exit

   for Mu in $(seq -3.00 0.50 -3.00)
   do
        mkdir -p Mu_$Mu
        cd Mu_$Mu || exit

        for trot_slices in $(seq 20 20 20)
        do
               mkdir -p L_$trot_slices
               cd L_$trot_slices || exit

                
               counter=0
               for s_val in $(seq 0 1 2)
               do
                      
	        	
                        mkdir -p Realization_$counter
                        cd Realization_$counter || exit
                                   

mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#SBATCH --job-name=l$trot_slices-Mu$Mu
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH -t 48:00:00                    
#SBATCH --partition trivedi
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL                        
#SBATCH --mail-user=roy.369@osu.edu                             
hostname
#SBATCH --no-requeue
module load gnu/9.1.0
cd \$SLURM_SUBMIT_DIR
time

cp ../../../../zeeman2db.f ./zeeman2db.f
cp ../../../../twod.in ./twod.in
cp ../../../../paramp.dat ./paramp.dat
sed -i 's/Slices/'$trot_slices'/g' ./twod.in
sed -i 's/Uval/'$U'/g' ./twod.in
sed -i 's/CHEM/'$Mu'/g' ./twod.in
sed -i 's/Realize/'$counter'/g' ./twod.in
sed -i 's/Slices/'$trot_slices'/g' ./paramp.dat


date +'%N' | tr -d 0 > seed
date +'%N' | tr -d 0 >> seed
gfortran -O3 -march=native zeeman2db.f -llapack -lblas -o bandmott
./bandmott &> output.out

time
EOF
)"
			 echo "$rawjob" &> job.bat
			 sbatch job.bat
			 counter=$(echo "$counter+1" | bc -l)
     
                         cd ..
                 done
                 cd ..
        done
        cd ..
        echo "Trot_Slice_$trot_slices, Mu_$Mu, Realization_$counter"
    done
    cd ..

done

