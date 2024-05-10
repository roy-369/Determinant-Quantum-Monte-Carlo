#!/bin/bash
Ntausk=10
chi_meas=1
dent_meas=1
warm_sweep=100
meas_sweep=500
Nwrap=4
Difflim=0.0001
Errrat=0.001
Doauto=0
Orthlen=12
Eorth=1.d-10
Dopair=1
Num_pair=1
Num_try=0
Strt_type=0
t=1
Dtau=0.05

for U in 1.0
do 
    mkdir -p U_$U
    cd U_$U || exit

    for Mu in 0.0
    do
        mkdir -p Mu_$Mu
        cd Mu_$Mu || exit

        for V in 1.0
        do
            mkdir -p V_$V
            cd V_$V || exit

            for trot_slices in 20
            do
                echo "Dtau $Dtau"
                mkdir -p dtau_${Dtau}_L_${trot_slices}
                cd dtau_${Dtau}_L_${trot_slices} || exit

                for counter in 1 2
                do
                    mkdir -p Realization_$counter
                    cd Realization_$counter || exit

mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#SBATCH --job-name=l$trot_slices-Mu$Mu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 8:00:00                    
#SBATCH --partition trivedi
#SBATCH --exclude=u012,u013,u014
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL                        
#SBATCH --mail-user=roy.369@osu.edu                             
hostname
#SBATCH --no-requeue
module load gnu/9.1.0
cd \$SLURM_SUBMIT_DIR
time

outname=z8l${trot_slices}u${U}dt${Dtau}mu${Mu}V${V}r${counter}
job_file=job8l${trot_slices}u${U}dt${Dtau}mu${Mu}V${V}r${counter}.dat

cp ../../../../../negu2009.f ./negu2009.f
cp ../../../../../param-negu2009.dat ./param-negu2009.dat
sed -i'' -e 's/Slices/'$trot_slices'/g' ./param-negu2009.dat
                    

echo $outname > $job_file
echo $Ntausk $chi_meas $dent_meas >> $job_file
echo $RANDOM >> $job_file
echo $t $Mu $V $Dtau >> $job_file
echo $warm_sweep $meas_sweep >> $job_file
echo $U >> $job_file
echo $Nwrap $Difflim $Errrat >> $job_file
echo $Doauto $Orthlen $Eorth $Dopair $Num_pair >> $job_file
echo $Num_try >> $job_file
echo $Strt_type >> $job_file


gfortran -O3 -march=native negu2009.f -llapack -lblas -o bandmott
./bandmott < $job_file 2>&1 | tee output.log
sleep 10
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
      done
cd ..
done 
                    
