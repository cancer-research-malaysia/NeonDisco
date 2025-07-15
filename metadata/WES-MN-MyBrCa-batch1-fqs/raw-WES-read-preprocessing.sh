# get a list of read1 fastqs off s3

aws s3 ls s3://crm.sequencing.raw.data.sharing/batch1/SLX-15056/ | grep SLX-15056.i709_i517 | awk '{print $NF}'|
grep 'r_1.fq.gz$' > s3_extracted_SD0033_MNw_R1_fastqs.txt

aws s3 ls s3://crm.sequencing.raw.data.sharing/batch1/SLX-15056/ | grep SLX-15056.i709_i517 | awk '{print $NF}'|
grep 'r_2.fq.gz$' > s3_extracted_SD0033_MNw_R2_fastqs.txt

# download read 1 files and read 2 files off s3 separately based on the txt list

while read -r line; do echo "$line"; aws s3 cp s3://crm.sequencing.raw.data.sharing/batch1/SLX-15056/"$line" .; done < s3_extracted_SD0033_MNw_R1_fastqs.txt

while read -r line; do echo "$line"; aws s3 cp s3://crm.sequencing.raw.data.sharing/batch1/SLX-15056/"$line" .; done < s3_extracted_SD0033_MNw_R2_fastqs.txt

# then merge the split read 1 and read 2 files into R1 and R2 files using cat

while read -r line; do echo "$line"; done < s3_extracted_SD0033_MNw_R1_fastqs.txt | xargs cat > SD0033_WES-MNw_R1.fq.gz

while read -r line; do echo "$line"; done < s3_extracted_SD0033_MNw_R2_fastqs.txt | xargs cat > SD0033_WES-MNw_R2.fq.gz






