for x in `seq 20 8 96`; do np ../../scatter.py model=Fishnet xholesize=${x}u yholesize=${x}u cornerradius=10u slabcdist=40u slabthick=4u; ../../effparam.py;  done
