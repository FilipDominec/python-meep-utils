for r in `seq 500 500 5000`
do 
	mpirun -np 2 ../../scatter.py model=PlasmonicDimers radius=${r}e-12 comment=dimer cellsizey=10e-9
done
