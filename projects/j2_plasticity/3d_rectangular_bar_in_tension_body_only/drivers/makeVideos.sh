#! /bin/bash

cd ../results/numRefinements0
mkdir -p videos
cd plots
mencoder mf://plot_displacement_time*.png -mf w=2880:h=1620:fps=60:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o ../videos/plot_displacement.avi
mencoder mf://plot_displacement_topdown_time*.png -mf w=2880:h=1620:fps=60:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o ../videos/plot_displacement_topdown.avi

cd ../../numRefinements1
mkdir -p videos
cd plots
mencoder mf://plot_displacement_time*.png -mf w=2880:h=1620:fps=60:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o ../videos/plot_displacement.avi
mencoder mf://plot_displacement_topdown_time*.png -mf w=2880:h=1620:fps=60:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o ../videos/plot_displacement_topdown.avi
