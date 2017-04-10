for r in 125 250 500 # 1000 2000
do 
	mpirun -np 2 ../../scatter.py model=PlasmonicDimers gap=${r}e-12 comment=dimer cellsizey=10e-9 radius=3e-9 resolution=125p
done
