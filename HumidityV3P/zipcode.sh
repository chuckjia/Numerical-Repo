#!/bin/bash

cd ~/Documents/Workspace/Git/Numerical-Repo/
printf "Check if result folder is cleared!\nEnter to proceed: " 
read act
if [[ $act == "Y" ||  $act = "y" ]]
then
zip -r -q ~/Downloads/humiditycode.zip HumidityV3 HumidityV3_MATLAB
echo "Successfully ziped. Results are in: ~/Downloads/humiditycode.zip"
else
echo "Aborted. Nothing is created."
fi
