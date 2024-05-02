#find /root/sequencing/storage/unicycler/ -type f -name "unicycler.log" -exec grep "best$" '{}' ';'
find /root/sequencing/storage/unicycler/ -type f -name "unicycler.log" -exec grep -Hn "best$" '{}' ';'

find /root/sequencing/storage/unicycler/ -type f -name "unicycler.log" -exec grep -Hn "best$" '{}' ';' | grep -v 127

