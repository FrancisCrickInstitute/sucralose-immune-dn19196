module purge
module load Nextflow/19.10.0
module load Singularity/2.6.0-foss-2016b

nextflow run nf-core/ampliseq \
	 -profile crick \
	 -resume \
	 --reads "fastq" \
	 --FW_primer CCTACGGGNGGCWGCAG \
	 --RV_primer GACTACHVGGGTATCTAATCC \
	 -c ampliseq.config \
	 --trunclenf 280 --trunclenr 210
#    --metadata "data/Metadata.tsv"
