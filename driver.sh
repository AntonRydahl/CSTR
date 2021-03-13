if [ $# -ne 2 ]
then
    echo "Incompatible number of arguments supplied. Please run:"
    echo "	./driver.sh <number of simulations> <number of threads>"
    echo "For instance"
    echo "	./driver.sh 100 8"
    echo "will run 100 simulations using 8 threads."
    exit 0
fi

make clean
make
export OMP_NUM_THREADS=$2
./project $1
export MYHOME=`pwd`
echo $MYHOME
VAR1="/driver.m"
VAR2="$MYHOME$VAR1"
echo "Executing file on path:"
echo $VAR2
matlab -nodesktop -nosplash -sd $MYHOME -r "driver($1);cd figures;fig_2_pdf();exit;"
