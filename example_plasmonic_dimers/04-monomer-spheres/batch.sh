for r in `seq 500 500 5000`
do 
	../../scatter.py model=PlasmonicDimers radius=${r}e-12 comment=singlesphere cellsizey=10e-9
done
