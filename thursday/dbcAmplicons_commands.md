Running the dbcAmplicons pipeline
===============================================

**ALL of this should only be done in an interactive session on the cluster**

srun -t 1-0 -c 4 -n 1 --mem 10000 --reservation workshop --pty /bin/bash

The goal of today is to process raw Illumina sequence reads to abudance tables for the 16sV1-V3 amplicon. To do so 1) we need have all the software installed, and 2) The Illumina sequence data in our project folder (mca_example). We will first need to prepare the input metadata files: barcodes, primers, and samples. Amplicon processing with dbcAmplicons includes the following steps: preprocessing, join, classify and abundances

<img src="thursday/Workflow.png" width="200">

Over the next week finish processing the rest of the amplicons.

**1\.** First lets validate our install and environment

Change directory into the workshops space

	cd
	cd mca_example
	ls --color

you should see 3 directories: bin, Illumina_Reads and src

Lets verify the software is accessible

	dbcVersionReport.sh

you should see the version info for dbcAmplicons, flash2 and RDP.

Next lets make a metadata directory and transfer our barcode, primer and sample sheet to their

	mkdir metadata

We can pull down the barcode and primer table from github

	cd metadata
	wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-September-Microbial-Community-Analysis-Workshop/master/metadata/dbcBarcodeLookupTable.txt
	wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2017-September-Microbial-Community-Analysis-Workshop/master/metadata/PrimerTable.txt	

You will also need to transfer the sample sheet you created into the metadata directory 

Once all the metadata tables are in the metadata folder, lets go back to the main workshop folder

	cd ..

If all this is correct, we are ready to begin.

**1\.** First lets validate our input files:

	dbcAmplicons validate -h
	dbcAmplicons validate -B metadata/dbcBarcodeLookupTable.txt -P metadata/PrimerTable.txt -S metadata/workshopSamplesheet.txt

If there are any errors, fix them (can do so in nano) and validate again.

Once we are sure our input files are ok

**2\.** Lets Preprocess the data 

	dbcAmplicons preprocess -h

First lets 'test' preprocessing, by only running the first 'batch' of reads

	dbcAmplicons preprocess -B metadata/dbcBarcodeLookupTable.txt -P metadata/PrimerTable.txt -S metadata/workshopSamplesheet.txt -O Slashpile.intermediate -1 Illumina_Reads/Slashpile_only_R1.fastq.gz --test > preprocess.log

View preprocess.log and the file Identified_barcodes.txt, make sure the results make sense

	cat preprocess.log
	cat Slashpile.intermediate/Identified_Barcodes.txt

Lets see what it looks like when you get the primer orientation incorrect. Try running the above but add the parameter, --I1 asis. The look at the output again.

Now run all reads

	dbcAmplicons preprocess -B metadata/dbcBarcodeLookupTable.txt -P metadata/PrimerTable.txt -S metadata/workshopSamplesheet.txt -O Slashpile.intermediate -1 Illumina_Reads/Slashpile_only_R1.fastq.gz > preprocess.log

Again view the output to make sure it makes sense

	cat preprocess.log
	cat Slashpile.intermediate/Identified_Barcodes.txt

Finnally, look at the output in the Slashpile.intermediate folder, how many subfolders are there? What do these coorespond to? What is inside each folder? View a few reads in one of the files.

**3\.** Lets merge/join the read pairs 

	dbcAmplicons join -h

	dbcAmplicons join -t 4 -O Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3 -1 Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3_R1.fastq.gz > join-16sV1V3.log
	cat join-16sV1V3.log

Try changing the parameter --max-mismatch-density, first to 0.1, then to 0.5, how do they differ.

If you prefer the default, run join again with defaults.

**4\.** Classify the merged reads using RDP

	dbcAmplicons classify -h

	dbcAmplicons classify -p 4 --gene 16srrna -U Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3.extendedFrags.fastq.gz -O Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3

classify produces a fixrank file, view the first 6 lines of the output file

	head Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3.fixrank

For the other amplicons what parameters in classify would need to be changed?

During the next week try running classify on the original preprocessed reads (pairs), how do the results compare to the overlapped set?

---

**5\.** Produce Abundance tables

	mkdir Slashpile.results

	dbcAmplicons abundance -h

	dbcAmplicons abundance -S metadata/workshopSamplesheet.txt -O Slashpile.results/16sV1V3 -F Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3.fixrank --biom > abundance.16sV1V3.log
	cat abundance.16sV1V3.log

Try changing the -r parameter see what changes? what about the -t parameter? Once done playing rerun the above to get the final biom file for the next phase of analysis.

**6\.** Split Reads by samples, for downstream processing in another application (post preprocessing/merging), or for submission to the SRA.

	splitReadsBySample.py -h
	splitReadsBySample.py -O SplitBySample/16sV1V3 -1 Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3_R1.fastq.gz -2 Slashpile.intermediate/16sV1V3/Slashpile-16sV1V3_R2.fastq.gz

---

**7\.** Write out the softare versions to a file in the folder.

	dbcVersionReport.sh &> > VersionInfo.txt

**8\.** Now perform 'join' and then remainder of the pipeline on the other amplicons: 16sV3V4, ITS1, ITS2, LSU

Post joining evaluate the results.

Given the results, are their any changes you would make to downstream for any of the amplicons?

Once join is complete for all amplicons, setup and run the remainder of the pipeline for each amplicon.

**Advanced** Set up a slurm script to sbatch the job on the cluster.

---

**FYI** Preprocessing Pairs with inline BC, Mills lab protocol

	preproc_Pair_with_inlineBC.py -h
