data=out/test_s.trimmed.truncated.gz
ref=data/endo/horse_chrom31.fa
bwa index ${ref}
bwa aln -l 15000 -t 10 ${ref} ${data} > out/aln.sai
bwa samse ${ref} out/aln.sai ${data} | samtools sort -O BAM -@ 20 -T test -  > out/aln.bam
