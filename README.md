# PacBio_seq_analysis

#Use genbank file of the template to get the barcode and library site annotation: 
1) Open genbank file in txt, note down features and locations to edit as follows
2) Modify genbank file annotation with edit_gb_file to replace misc_feature with barcode. Simply open genbank as text and see which feature you need to edit. If it is 2nd, you will use index 1 as counting starts from 0
example, this is the function definition def change_feature_type(template, i, strtloc, endloc, new_type):
template is gb file, i is the index of the feature, strtloc is the start of the feature as mentioned in gb file, endloc is end, new_type is the new annotation/new feature type
3) Open gb_feature_edit_parameters and input any number of changes you want to make

#Minimap2 to generate SAM file:

1) install minimap2 via pip
2) for the sam file, use minimap2 to create consensus sequence. Use fasta file for template
3) Run minimap2 after editing features: minimap2 --MD -Lax map-pb template_file.fasta fastq_data.fastq > file_to_create.sam

#Specify template, sam file locations:
1) Input the location for genbank template and also the sam file
2) Add bam file and bai file location (only the directory for these files has to exist, the code will generate these files)

#Functions in the main file:

1) sort_sam = True       #sort sam file and create bam file
   index_file = True     #index bam file

   Be sure to install pysam 

2) barcode_filter = True #Check qscore in PBQscore, all barcode reads with Q-score >= 30 are considered

   Specify the total number of reads (tot) to parse through in the PBQscore file and length of the barcode (bclen)

3) filtered_Reads_file = True       #sort and index filtered barcode reads

4) aaRS_filter = True               #filter_aaRS reads
   qcutoff = 20 #qscore cutoff for aaRS reads

   Be sure to specify the qcutoff which is the minimum Qscore you want

5) mutations = True                 #generate a dataframe of aaRS mutations and corresponding barcodes

   Specify the sequence of the synthetase at the amino acid as well as DNA level (MaPylRS AA and MaPylRS DNA respectively)
    
6) bc_mut = True                   #generate a dataframe for number of barcodes per unique mutant
   mut_bc = True                   #generate a dataframe for number of unique mutants per barcode

   Generate graphs using matplotlib, make sure matplotlib is installed



