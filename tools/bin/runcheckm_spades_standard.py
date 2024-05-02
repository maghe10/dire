#checkm data setRoot /root/sequencing/checkm/data

#checkm lineage_wf -r --max_memory 16.GB -t 1 -x fna /root/sequencing/checkm/data/test_data /root/sequencing/checkm/data/output
#--reduced-tree


#(M) > checkm tree <bin folder> <output folder>
#(R) > checkm tree_qa <output folder>
#(M) > checkm lineage_set <output folder> <marker file>
#(M) > checkm analyze <marker file> <bin folder> <output folder>
#(M) > checkm qa <marker file> <output folder>

#checkm tree -r -t 40 -x fna /root/sequencing/checkm/data/test_data /root/sequencing/checkm/data/output

import os

os.chdir("/root/sequencing/tools/bin")
os.system("mkdir -p /root/sequencing/intermediate/checkm_spades_standard/")
os.system("checkm lineage_wf -r -t 1 -x fasta /root/sequencing/dire/spades_standard/assembly /root/sequencing/intermediate/checkm_spades_standard/")
os.system("checkm qa /root/sequencing/intermediate/checkm_spades_standard/lineage.ms /root/sequencing/intermediate/checkm_spades_standard/ --tab_table -f /root/sequencing/intermediate/checkm_spades_standard/qa_spades_standard.tsv")
