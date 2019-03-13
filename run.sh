# script to run the python program, save the output to a log file and save the log file. 

conda activate DNA

#pythonfile=$1
pythonfile="main.py"
date_str="$(date '+%Y.%m.%d_%H.%M.%S')"
log_str="${date_str}_log.log"
python_str="${date_str}_${pythonfile}"

cd src
cp $pythonfile "../logs/${python_str}"
#python -u ${pythonfile} 2>&1 |  tee "CM_logs/${log_str}"  & disown
#python -u ${pythonfile} 2>&1 | grep -vP --line-buffered '(?!.*total= )(?=\[CV\])^.*$' | grep -vP --line-buffered '\[LightGBM\] \[Warning\].*will be ignored\. Current value' | tee "CM_logs/${log_str}" & disown

python -u ${pythonfile} 2>&1 | tee "../logs/${log_str}" & disown
cd ..

rm out.log
ln -s "logs/${log_str}" out.log
