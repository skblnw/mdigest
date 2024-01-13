for ii in ILAKFLHWL-ila1.5men LGYGFVNYI-a6.3pwp LLLGIGILV.7zuc RMFPNAPYL-a7b2.6rsy SLLMWITQC-1g4.2p5e; do
    (cd "$ii/hlaonly/analysis/mdigest/t1" && rm gcc_submatrix.txt distances_submatrix.txt && python ../run.py > "runlog.txt" 2>&1 &)
    (cd "$ii/hlaonly/analysis/mdigest/t2" && rm gcc_submatrix.txt distances_submatrix.txt && python ../run.py > "runlog.txt" 2>&1 &)
done
