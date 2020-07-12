Here are the pipeline for updating blast database. We will need to do this every Jan, Apr, Jul, Oct. It is done this year, so next update is Jan 2013.

0.	Add 7zip to path:
a.	/mnt/cluster/tools/p7zip_9.20.1/bin

1.	Back up the old nr and virus database:

mv /mnt/san/cluster/xdeng/blastdb/nr /mnt/san/cluster/xdeng/blastdb/nr_today
mv /mnt/san/cluster/xdeng/blastdb/virus /mnt/san/cluster/xdeng/blastdb/virus_today

2.	Download the virus genome and nr to current directory /mnt/san/cluster/xdeng/blastdb/

wget ftp://ftp.ncbi.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
3.	Update taxonomy:
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_prot.zip
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxcat.zip
7za e taxcat.zip
7za e gi_taxid_prot.zip

4.	Create new empty directory
a.	mkdir /mnt/san/cluster/xdeng/blastdb/nr
b.	mkdir /mnt/san/cluster/xdeng/blastdb/virus
5.	Process the downloaded files 
                python nr_virus.py
6.	Run makeblastdb
a.	cd nr 
b.	makeblastdb -in non_virusNR.fa -dbtype prot -parse_seqids -out nr
c.	cd virus
d.	makeblastdb -in virus.fa -dbtype prot -parse_seqids -out virus
7.	Change mode
a.	sudo chmod 777 * -R


======================================================================================================
Instruction for running virus discovery pipeline
0 (set up one time only):
                setup path for the following tools
                /mnt/cluster/xdeng/script
                /usr/bin/
                /mnt/cluster/tools
1. download the data into directory:
                /mnt/cluster/xdeng/Eric/NewProject/fastq
        wget -r -l1 --no-parent -A.gz --user='chiulab' --password='!wantdata' vddc.ucsf.edu/~miseq/miseq_data/PT_121407_actual/
                change the mode to be rwe for all files in the project directory
2. create sample file inside fastq folder:
                ls -1 *.gz >samples.txt
3. go to project directory: 
                cd ..
                run readseeds2.py to generate pipeline script
                sudo chmod 777 * -R
                  #sh pipeline_run.sh
4. Now you see several scripts created: run them in order
                sh polyA.sh
                   sh clone.sh
                sh trim.sh
                sh bowtieNT.sh  # for RNA seq
                sh soap.sh
                sh combineContig_reads.sh

                //start this for overnight run
                //use top.sh and kill.sh to monitor progress
                sh blast_virus.sh
                sh blast_virus_parser.sh

                //start this for overnight run
                //use top.sh and kill.sh to monitor progress

                sh blast_nr.sh
                sh blast_nr_filter.sh
                sh blast_output_sort.sh

Also reinstall blast software every Jan and July or if there is a new version
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
