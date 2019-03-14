#!/bin/bash
if [ $# != 3 ]
then
        echo "USAGE: $0 [FASTA file path] [Read length] [Absolute path to bowtie2 index + prefix]"
        exit 0;
fi

fasta_fp=$1;
l_read=$2;
bt2_prefix=$3;

if [ ! -e $fasta_fp ]
then
	echo "Could not find fasta file ${fasta_fp}"
	exit 0;
fi

echo "Creating temporary files directory"
rm -f -r temp*
mkdir temp

MUSIC -preprocess_FASTA ${fasta_fp} temp

if [ ! -e temp/chr_ids.txt ]
then
	echo "Could not find chromosome id's file path."
	exit 0;
fi

chr_ids=(`cat temp/chr_ids.txt`);

# Process all the chromosomes.
echo "Writing mapping commands."
rm -f temp_map_reads.csh
for chr_id in ${chr_ids[@]}
do
	mkdir temp/${chr_id}
	echo "MUSIC -fragment_sequence_2_stdout temp/${chr_id}.bin ${l_read} read | bowtie2 -x ${bt2_prefix} -k 5 -p 4 -f -S /dev/stdout -U - | MUSIC -preprocess SAM stdin temp/${chr_id}" >> temp_map_reads.csh
done
chmod 755 temp_map_reads.csh

# Merge all the reads per chromosome.
echo "Writing mapped read processing commands."
rm -f temp_process_mapping.csh
for chr_id in ${chr_ids[@]}
do
	echo "find temp -name '${chr_id}_mapped_reads.txt' | xargs -Ifiles cat files | MUSIC -get_multimapability_signal_per_mapped_reads stdin ${chr_id}.bin ${l_read}" >> temp_process_mapping.csh
done
chmod 755 temp_process_mapping.csh

echo "Done. Run the scripts named temp_map_reads.csh then temp_process_mapping.csh."

