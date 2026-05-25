from confindr_src import confindr

# Find read files.
paired_reads = confindr.find_paired_reads('/root/dire/data/Illumina/trimmed', forward_id='_1', reverse_id='_2')
# Run confindr. This assumes that you have already downloaded the databases. If you haven't,
# you can run confindr.check_for_databases_and_download(database_location='path/where/you/want/to/download, tmpdir='a/tmp/dir')
for pair in paired_reads:
    print(pair)
    confindr.find_contamination(pair=pair,
                                forward_id='_1', # change if yours is different
                                threads=4,
				base_cutoff=3,	 
                                output_folder='/root/sequencing/intermediate/confindr',
                                databases_folder='/root/.confindr_db')