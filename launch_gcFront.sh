#!/bin/bash
export PATH=$PATH:/home/alvaro/Matlab_2022/bin/

nohup matlab -nodisplay -nodesktop -r "run('/home/alvaro/github/ConsortiumEngineeringTool/jupyter_notebook/run_gcFront.m')" > /home/alvaro/github/ConsortiumEngineeringTool/jupyter_notebook/gcfront.txt &