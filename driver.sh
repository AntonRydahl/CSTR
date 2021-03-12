if [ $# -ne 1 ]
then
    echo "Incompatible number of arguments supplied. Please run:"
    echo "./driver.sh <number of simulations>"
    echo "For instance:"
    echo "./driver.sh 10"
    exit 0
fi

make clean
make
export OMP_NUM_THREADS=20
./project $1
export MYHOME=`pwd`
echo $MYHOME
VAR1="/driver.m"
VAR2="$MYHOME$VAR1"
echo "Executing file on path:"
echo $VAR2
matlab -nodesktop -sd $MYHOME -r "driver($1);pause(10);exit;"

rm *.txt
