for r in 0 125 250 500 1000 2000
do 
	mpirun -np 2 ../../scatter.py model=PlasmonicDimers gap=${r}e-12 comment=dimer cellsizey=5e-9 radius=1e-9 resolution=250p
done
