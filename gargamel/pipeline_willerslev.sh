
#OUTPREFIX="out/test"
#mkdir -p $(dirname $OUTPREFIX)
#software/gargammel/gargammel.pl -n 1000000 -se --comp 0,0,1 -damage 0.03,0.4,0.01,0.3 -o ${OUTPREFIX} IN_TRIMMED/
### -n 1000 is the number of reads simulated
### -se -> singl-end. for paired-end just remove this argument
### --comp 0,0,1 sample p(exogenous|microbial) == 0, p(contamination) == 0, p(hostDNA) == 1
### -damage 0.03,0.4,0.01,0.3 v: nick frequency l: length of overhanging ends (geometric parameter) d: prob. of deamination of Cs in double-stranded parts s: prob. of deamination of Cs in single-stranded parts

print_time_stamps() {
	dt=$(echo "$2 - $1" | bc)
	days=$(echo "$dt/86400" | bc)
	dt2=$(echo "$dt%86400" | bc)
	hours=$(echo "$dt2/3600" | bc)
	dt3=$(echo "$dt2%3600" | bc)
	minutes=$(echo "$dt3/60" | bc)
	seconds=$(echo "$dt3%60" | bc)
	printf "Total runtime: %d:%02d:%02d:%02.4f\n" $days $hours $minutes $seconds
}


[ -z "$1" ] && echo "No argument supplied" && exit 1


INFILE_W_EXT=$1
L=${#INFILE_W_EXT}-6
INFILE=${INFILE_W_EXT:0:${L}}

ref=${2:-"horse_chrom31.fa"}
OUTDIR="./out/${INFILE}"

printf "\n"
printf "Pipeline: Trimming, aligning and extracting the relevant columns of ${INFILE_W_EXT}.\n"
printf "Using ${ref} as reference and outputting to ${OUTDIR}."


mkdir -p ${OUTDIR}
OUT_TRIMMED="${OUTDIR}/${INFILE}.trimmed"

time0=$(date +%s.%N)


# trimming$
printf "\n\n"
printf "Removing Adapters: \n\n"
AdapterRemoval --file1 ${INFILE_W_EXT} --threads 3  --basename ${OUT_TRIMMED} --trimns --trimqualities --gzip


IN_TRIMMED="${OUT_TRIMMED}.truncated.gz"

OUTFILE="${OUTDIR}/${INFILE}"
OUT_SAI="${OUTFILE}.sai"
OUT_BAM="${OUTFILE}.bam"
OUT_TXT="${OUTDIR}/gargamel_${INFILE}.txt"
mapDamage_dir="${OUTDIR}/mapDamage"


# creating index of reference file
printf "\n\n"
printf "Creating index of ref file: \n\n"
bwa index ${ref}

# align trimmed read to reference using 20 threads
printf "\n\n"
printf "Aligning trimmed reads: \n\n"
bwa aln -l 15000 -t 20 ${ref} ${IN_TRIMMED} > ${OUT_SAI}

# create bam file
#bwa samse ${ref} ${OUT_SAI} ${IN_TRIMMED} | samtools sort -O BAM -@ 20 -T test -  > ${OUT_BAM}
printf "\n\n"
printf "Creating bam file: \n\n"
bwa samse ${ref} ${OUT_SAI} ${IN_TRIMMED} | samtools view -Sb -@ 20 > ${OUT_BAM}

# extract correct columns of file
printf "\n\n"
printf "Extracting correct columns of file. "
samtools view -F 4 ${OUT_BAM} | cut -f2,6,10,12-100 > ${OUT_TXT}

printf "\n\n"
print_time_stamps $time0 $(date +%s.%N)

#mapDamge of simulated gargamel sample
#mapDamage -i ${OUT_BAM} -r ${ref} -d ./mapDamage_gargamel/
printf "\n\n"
printf "Running mapdamage: \n\n"
mapDamage --no-stats -i ${OUT_BAM} -r ${ref} -d ${mapDamage_dir}

printf "\n\n"
print_time_stamps $time0 $(date +%s.%N)

printf "\n\n"
printf "Deleting temporary files. \n\n"
shopt -s extglob
rm !(*.gz|*.fa|*.sh|*.fq|out)
cd $OUTDIR
rm !(*.txt|*.bam|mapDamage)
cd ../..
shopt -u extglob
