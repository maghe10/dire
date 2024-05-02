wget https://figshare.com/ndownloader/files/41228577 -O test_samples.tar.gz && \
    tar -xzvf test_samples.tar.gz && \
    rm test_samples.tar.gz
confindr -i test_samples -o test_out
confindr -i test_samples -o test_out --rmlst