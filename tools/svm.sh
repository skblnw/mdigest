#!/bin/bash

filelist = "training_data.txt"

for test in a b c; do
    echo "path/to/txt" > $filelist

    for training in a b c; do
        if [ $test != $training ]; then
            echo "path/to/txt" >> $filelist
        fi
    done

    python svm.py --filelist $filelist --test_data "/path/to/txt" > "svm_result_$test.txt" 2>&1
done
