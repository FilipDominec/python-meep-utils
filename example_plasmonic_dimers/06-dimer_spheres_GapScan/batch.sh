for r in 250 # `seq 500 500 5000`
do 
	mpirun -np 2 ../../scatter.py model=PlasmonicDimers gap=${r}e-12 comment=dimer cellsizey=10e-9 radius=3e-9
done
