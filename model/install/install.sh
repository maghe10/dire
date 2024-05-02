mkdir -p /root/sequencing/dire/model/dire/saved_models
cp /root/sequencing/dire/model/data/conformal_set_random.p /root/sequencing/dire/model/dire/saved_models
unzip /root/sequencing/dire/model/data/bert_updated.zip -d /root/sequencing/dire/model/dire/saved_models

mkdir -p /root/sequencing/dire/model/dire/data
cp /root/sequencing/dire/model/data/vocabulary.pth /root/sequencing/dire/model/dire/data

conda env create -f /root/sequencing/dire/model/dire/environment_ubuntu_eval.yml

conda activate dml




