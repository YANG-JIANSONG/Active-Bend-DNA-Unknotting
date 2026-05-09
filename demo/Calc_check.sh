for type in 0 31;do

   cd $type
   cp -r ../source/check.py source/check.py 
   cp -r ../source/merge_overall_tri.sh source/merge_overall_tri.sh
    cd L300
    cp -r ../source/check.py check.py 
    sed -i '$d' check.py && echo "check($type)" >> check.py
    cp -r ../source/merge_overall_tri.sh merge_overall_tri.sh
    # bash merge_overall_tri.sh
    > output_all_save_check.txt
    > trj_reco_check.xyz

    nohup python check.py &
    cd ..
    cd ..
done
    
