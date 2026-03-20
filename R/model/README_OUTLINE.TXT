##################################################
#  README model call. Juan 2025 model
#
##################################################


1. Generate input files
A Generate input files in Anders format
# source(file = 'clear.R')
# source(file = 'common/model_input.R')
# ALL()
# => creates files in model/input directory e.g.  "sirAntibioticsModelWords_<mode>.csv" and "modelOutput<n>.csv", "combinations<n>.csv"
B Convert "sirAntibioticsModelWords_<mode>.csv" files into juan format
# Output: input/<mode>/"sirAntibioticsModelWordsJuan_<n>_<patientindex>.csv"


2. Run juans model (python) in ubuntu
#cd /root/dire/program/model/bin
#conda activate juan
#python run_all_models_cli.py 
#    which calls "run_model_cpu.py" which is a CPU port of juans code 
#    located in /root/dire/program/model/Confidence-based-Prediction-of-Antibiotic-Resistance
# Output: output/<version>/<mode>/predsraw/"sirAntibioticsModelWordsJuan_<n>_<patientindex>.csv"

3. Run conformal predictions (RScript.exe) in anaconda prompt on windows
# cd git/dire/model/bin
# python run_all_conformal_cli.py
# Output: output/<version>/<mode>/predscp/"sirAntibioticsModelWordsJuan_<n>_<patientindex>_cp.csv"

4. Post-process conformal predictions files
A. To same format as when running anders model (modelOutput_answer<n>.csv, modelOutput_preds<n>.csv, modelOutput_preds<n>.csv )
# source(file = 'clear.R')
# source(file = 'model/convert_to_preds.R')
# ALL()
# Output: output/<version>/<mode>/predsalt/modelOutput_preds_<n>.csv etc.
B To merged long data sets
# Output output/<version>/<mode>/predsalt/modelOutput_long_<n>.csv etc.
# Output output/<version>/<mode>/predsalt/modelOutput_long_merged.rds etc.

5. Make statistics from the output/<version>/<mode>/predsalt/modelOutput_long_merged.rds files.
#source(file = 'clear.R')
#source(file = 'model/prediction_statistics.R')
#ALL()
#Output: output/<version>/<mode>/statistics/<prefix>allstats_<sum/sample/word>_<normal/strict>.csv
#Output: output/<version>/common/statistics/<prefix>allstats_<sum/sample/word>.rds
# The latter had "mode" and "cpmode" as columns. Long formats.








# source(file = 'clear.R')
# source(file = 'model/prediction_statistics.R')

  
