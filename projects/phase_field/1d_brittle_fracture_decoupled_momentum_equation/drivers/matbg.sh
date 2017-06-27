#! /bin/bash
module load matlab/r2012a
# Clear the DISPLAY.
unset DISPLAY  
# unset DISPLAY for some shells

# Call MATLAB with the appropriate input and output,
# make it immune to hangups and quits using ''nohup'',
# and run it in the background.
nohup nice -10 matlab -nodisplay -nodesktop -nojvm -nosplash < $1 > $2 &
