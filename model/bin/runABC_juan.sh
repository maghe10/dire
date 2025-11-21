#!/bin/bash

# Loop from 4 to 13
for i in {4..13}
do
  echo "Running model C with number_of_choices = $i"

  python ./runmodel_juan.py \
    --number_of_choices "$i" \
    --significant_level 0.10 \
    --output_folder "/root/dire/data/Analyser/processed/R/model/output/temp/Mode-A/preds" \
    --load_sir_csvfile "sirAntibioticsModelWordsJuan2020_12_Mode-A.csv"

  echo "Completed run C for number_of_choices = $i"
  echo "-----------------------------------------"
done


for i in {4..13}
do
  echo "Running model B with number_of_choices = $i"

  python ./runmodel_juan.py \
    --number_of_choices "$i" \
    --significant_level 0.10 \
    --output_folder "/root/dire/data/Analyser/processed/R/model/output/temp/Mode-B/preds" \
    --load_sir_csvfile "sirAntibioticsModelWordsJuan2020_12_Mode-B.csv"

  echo "Completed run C for number_of_choices = $i"
  echo "-----------------------------------------"
done

for i in {4..13}
do
  echo "Running model C with number_of_choices = $i"

  python ./runmodel_juan.py \
    --number_of_choices "$i" \
    --significant_level 0.10 \
    --output_folder "/root/dire/data/Analyser/processed/R/model/output/temp/Mode-C/preds" \
    --load_sir_csvfile "sirAntibioticsModelWordsJuan2020_12_Mode-C.csv"

  echo "Completed run C for number_of_choices = $i"
  echo "-----------------------------------------"
done
