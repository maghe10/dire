mkdir -p /root/dire/program/model/dire/saved_models
cp /root/dire/program/model/data/conformal_set_random.p /root/dire/program/model/dire/saved_models
unzip /root/dire/program/model/data/bert_updated.zip -d /root/dire/program/model/dire/saved_models

mkdir -p /root/dire/program/model/dire/data
cp /root/dire/program/model/data/vocabulary.pth /root/dire/program/model/dire/data

conda env create -f /root/dire/program/model/dire/environment_ubuntu_eval.yml

conda activate dml


