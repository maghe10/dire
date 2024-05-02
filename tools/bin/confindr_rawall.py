from confindr_src import confindr


from confindr_src import confindr

# Find read files.
paired_reads = confindr.find_paired_reads('/root/sequencing/in/reads', forward_id='_1', reverse_id='_2')
print(paired_reads)
# Run confindr. This assumes that you have already downloaded the databases. If you haven't,
# you can run confindr.check_for_databases_and_download(database_location='/root/sequencing/confindr/databases', tmpdir='/root/sequencing/confindr/tmp')
for pair in paired_reads:
    confindr.find_contamination(pair=pair,
                                forward_id='_1', # change if yours is different
                                threads=4,
                                base_cutoff=3, 
                                output_folder='/root/sequencing/intermediate/confindrraw_private',
                                databases_folder='/root/sequencing/confindr/databases')

