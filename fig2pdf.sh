export MYHOME=`pwd`
echo $MYHOME
VAR1="/fig_2_pdf.m"
VAR2="$MYHOME$VAR1"
echo "Executing file on path:"
echo $VAR2
matlab -nosplash -nodesktop -sd $MYHOME -r "addpath(genpath('../../Matlab_plot_tools'));fig_2_pdf();exit;"

exit 0
