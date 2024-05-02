cd /root/sequencing/storage
find trimmed_fastqc/ -name *.gz -exec ls -lv {} \; > trimmedfastq_ls.txt
find trimmed/ -name *.gz -exec ls -lv {} \;  > trimmed_ls.txt
awk '{ print $5 }' trimmedfastq_ls.txt > trimmedfastq_sizes.txt
awk '{ print $5 }' trimmed_ls.txt > trimmed_sizes.txt
sort trimmedfastq_sizes.txt > trimmedfastq_sorted.txt
sort trimmed_sizes.txt > trimmed_sorted.txt
diff trimmedfastq_sorted.txt trimmed_sorted.txt
