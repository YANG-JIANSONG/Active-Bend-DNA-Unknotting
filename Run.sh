total_step=500
save_step=1000
jump_step=500
# L=300
for type in 0 31;do
   rm -f -r $type
   mkdir  $type
   cd $type
   cp -r ../source source

for L in 300 ; do
      # echo Ls$a    
      rm -f -r L$L
    #   mkdir Lp$Lp
      cp -r ./source L$L
      cd L$L
      cp -r ../source source
    # Lp=1
    for a in  $(seq 1 10 );do

        cp -r ./source A$a
        echo A$a
        cd A$a
                > output_index_0921.txt
                > output_all_save_0921.txt
                > trj_reco_0921.xyz
        cp -r ../source source
        sed -i '$d' in.py && echo "Run_lammps(lp=17,seed=$a,L=$L,type=$type,degree=150)" >> in.py
        nohup python in.py &
        cd ..
done
        # bash If_right_click_Run.sh
        cd ..
done
        cd ..
done
