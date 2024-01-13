#!/bin/bash

filelist="training_data.txt"
matrix="gcc_submatrix.txt"

for test in ILAKFLHWL-ila1.5men LGYGFVNYI-a6.3pwp LLLGIGILV.7zuc RMFPNAPYL-a7b2.6rsy; do
    for ii in 1 2; do 
        echo "/home/kevin/data/maap/md.a2/SLLMWITQC-1g4.2p5e/hlaonly/analysis/mdigest/t1/$matrix" > $filelist

        for training in ILAKFLHWL-ila1.5men LGYGFVNYI-a6.3pwp LLLGIGILV.7zuc RMFPNAPYL-a7b2.6rsy; do
            if [ $test != $training ]; then
                echo "/home/kevin/data/maap/md.a2/$training/hlaonly/analysis/mdigest/t1/$matrix" >> $filelist
                echo "/home/kevin/data/maap/md.a2/$training/hlaonly/analysis/mdigest/t2/$matrix" >> $filelist
            fi
        done

        echo "Testing $test"
        echo "Training set contains: "
        cat $filelist
        python svm_skf.py --filelist $filelist --test_data "/home/kevin/data/maap/md.a2/$test/hlaonly/analysis/mdigest/t$ii/$matrix" > "svm_result_${test}_t${ii}.txt" 2>&1
    done
done
