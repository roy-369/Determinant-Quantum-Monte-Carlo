#!/bin/bash

# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

module load intel

Dtau=0.025

for U in 10.0
do
    mkdir -p U_$U
    cd U_$U

    for t_prime in -0.25
    do
        mkdir -p t_prime_$t_prime
        cd t_prime_$t_prime

        for Mu in $(seq 0.50 0.50 5.00)
        do
            mkdir -p Mu_$Mu
            cd Mu_$Mu || exit

            for trot_slices in $(seq 40 20 80)
            do
                 mkdir -p dtau_${Dtau}_L_$trot_slices
                 cd dtau_${Dtau}_L_$trot_slices || exit

                 for counter in $(seq 0 1 20)
                 do
                      
#	        	echo "$counter"
                     mkdir -p Realization_$counter
                     cd Realization_$counter || exit

           

mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#SBATCH --job-name=l$trot_slices-Mu$Mu
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH -t 48:00:00                    
#SBATCH --mem=48G
#SBATCH --mail-type=FAIL                        
#SBATCH --mail-user=roy.369@osu.edu                             
hostname
#SBATCH --no-requeue
module load intel
cd \$SLURM_SUBMIT_DIR
time
source ~/.bashrc
source ~/.bash_profile

cp ../../../../../square_temp.geom ./square.geom
cp ../../../../../ggeom ./ggeom
cp ../../../../../sq_temp.in ./sq_curr.in

sed -i 's/Uval/'$U'/g' square.geom
sed -i 's/t1/'$t_prime'/g' square.geom
sed -i 's/Uval/'$U'/g' sq_curr.in
sed -i 's/delt/'$Dtau'/g' sq_curr.in
sed -i 's/t1/'$t_prime'/g' sq_curr.in
sed -i 's/CHEM/'$Mu'/g' sq_curr.in
sed -i 's/Slices/'$trot_slices'/g' sq_curr.in
sed -i 's/realize/'$counter'/g' sq_curr.in

./ggeom sq_curr.in &> output.out

time
EOF
)"
			 echo "$rawjob" &> job.bat
			 sbatch job.bat
                         cd ..
                   done
                   cd ..
            done
            cd ..
       done
       cd ..
       echo "tprime_$t_prime, Trot_Slice_$trot_slices, Mu_$Mu, Realization_$counter"
   done
   cd ..   
  
done

