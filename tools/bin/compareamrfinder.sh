cd /root/sequencing/storage
awk '{ print $6 }' amrfinder/sample2.tsv > amrfinder_2_genes.txt
awk '{ print $6 }' amrfinder_optimal/sample2.tsv > amrfinderoptimal_2_genes.txt
sort amrfinder_2_genes.txt > amrfinder_2_genes_sortes.txt
sort amrfinderoptimal_2_genes.txt > amrfinderoptimal_2_genes_sortes.txt
diff amrfinderoptimal_2_genes_sortes.txt amrfinder_2_genes_sortes.txt

