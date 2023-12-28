#!/bin/bash


usage() {
    echo "Usage: make binary file from the bedgraph file. "
    echo " 	Normalized value: bins per million mapped reads (BPM) , depth (nFrag or ReadsInTss) and cell_num (per 100 cells)."
    echo "	bed_to_bdg.sh [-i bedgraph] [-o out_dir] [-f cutoff]  [-m mark] [-t cell_type/cluster]"
    echo "		      [-s bin_size] [-c chrom_size]"
    echo "Description:"
    echo "	-i bedgraph, the input bedgraph file path."
    echo "    	-o out_dir, the output dir name."
    echo "    	-f cutoff, used to binarye."
    echo "    	-m mark, the ChIPseq mark name, example H3K27me3."
    echo "    	-t cell_type / cluster name,"
    echo "    	-s bin_size, should the same with bedgraph file. Default:500."
    echo "    	-c chrom_size, the chroms size file."
    exit -1
}
 


bin_size=500

while getopts 'h:i:o:f:m:t:s:c:' opts
do
	case $opts in
		i) bdg="$OPTARG";;
		o) out_dir="$OPTARG";;
		f) cutoff="$OPTARG";;
		m) mark="$OPTARG";;
		t) cell_type="$OPTARG";;
		s) bin_size="$OPTARG";;
		c) chrom_size="$OPTARG";;
		h) usage;;
        	?) usage;;
    	esac
done


function binary_one_chr(){
	# $1 input_bdg, should the same with bin_size
	# $2 out, out file path
	# $3 mark
	# $4 cell type / cluster
	# $5 chr
	# $6 chrom_size file
	# $7 bin_size
	bdg=$1;out=$2;mark=$3;cell_type=$4;chr=$5;chrom_size=$6;bin_size=$7
	grep -w $chr $chrom_size|bedtools makewindows -g - -w $bin_size|\
	bedtools intersect -a - -b $bdg -wao|\
	awk -v cut=$cutoff '{if ($5 != -1 && $7 >cut){v=1}else{v=0};print v}'|\
	sed "1i$cell_type\t$chr\n$mark" > $out
}

if [ ! -d "$out_dir" ]
then
	mkdir -p $out_dir
fi

chroms=(`cut -f 1 $chrom_size|grep -v random|grep -v chrUn|grep -v hap|grep -v chrM|grep -v chrY|tr '\n' ' '`)

for chr in ${chroms[@]}
do
	o=${out_dir}/${cell_type}_${chr}"_binary.txt" 
	binary_one_chr $bdg $o $mark $cell_type $chr $chrom_size $bin_size &
done
wait

