#!/bin/bash

usage() {
    echo "Usage: convert bed to bedgraph file. "
    echo " 	Normalized value: bins per million mapped reads (BPM) , depth (nFrag or ReadsInTss) and cell_num_scale)."
    echo "	bed_to_bdg.sh [-i bed] [-o out_name] [-s bin_size] [-c chrom_size] [-l black_list]"
    echo "		      [-n norm_method] [-g 
    
    
    
    
    ] [-n cell_num_scale] [-b bin_merge]"
    echo "Description:"
    echo "	-i bed, the input bed file path."
    echo "    	-o out_name, the output file name."
    echo "    	-s bin_size, the bin size to bin genome. Default:100."
    echo "    	-c chrom_size, the chroms size file."
    echo "    	-l black_list, the blacklist file."
    echo "    	-m norm_method, the normalized method: nFrag or ReadsInTss. Default: ReadsInTss."
    echo "    	-g refgene, the refGene.txt file path."
    echo "    	-n cell_num_scale, the depth of cell number. Default: 1."
    echo "	-b bin_merge, keep the orginal bin(bin) or merge the bins into regions with "
    echo "	   same value and remove 0 bins(region). Default: region."
    exit -1
}
 


bin_size=100
norm_method="ReadsInTss"
cell_num_scale=1
bin_merge="region"
black_list=FALSE

while getopts 'h:i:o:s:c:l:m:g:n:b:' opts
do
	case $opts in
		i) bed="$OPTARG";;
		o) out_name="$OPTARG";;
		s) bin_size="$OPTARG";;
		c) chrom_size="$OPTARG";;
		l) black_list="$OPTARG";;
		m) norm_method="$OPTARG";;
		g) refgene="$OPTARG";;
		n) cell_num_scale="$OPTARG";;
		b) bin_merge="$OPTARG";;
		h) usage;;
        	?) usage;;
    	esac
done

out_dir=`cd "$(dirname $out_name)" && pwd`
out_name=`basename $out_name`
out_name=${out_dir}/$out_name
bdg=${out_name}.bdg
bw=${out_name}.bw

if [ $norm_method == "nFrag" ]
then
	seqDepth=`wc -l $bed|tr ' ' '\t'|cut -f 1`
elif [ $norm_method == "ReadsInTss" ]
then
    if [ $refgene == *.txt ]
    then
        seqDepth=`awk 'OFS="\t"{if ($4 == "+"){tss=$5}else{tss=$6};s=tss-2000;e=tss+2000;if (s<0){s=1};print $3,s,e}' $refgene|sort -k1,1 -k2,2n -S 10G --parallel=10 -u|bedtools intersect -a $bed -b - -wa|sort -k1,1 -k2,2n -S 10G --parallel=10 -u|wc -l|tr ' ' '\t'`
    elif [ $refgene == *.gtf ]
    then
        seqDepth=`awk 'OFS="\t"{if ($7 == "+"){tss=$4}else{tss=$5};s=tss-2000;e=tss+2000;if (s<0){s=1};print $1,s,e}' $refgene|sort -k1,1 -k2,2n -S 10G --parallel=10 -u|bedtools intersect -a $bed -b - -wa|sort -k1,1 -k2,2n -S 10G --parallel=10 -u|wc -l|tr ' ' '\t'`
    fi
fi


scale_factor=`echo "1000000 / $seqDepth" | bc -l`
scale_factor=`echo "scale=5;$scale_factor / $cell_num_scale" |bc -l`
echo "scale_factor:"$scale_factor";seqDepth: "$seqDepth";norm_method: "$norm_method

tmp1=${bdg}.tmp1
tmp2=${bdg}.tmp2
bedtools makewindows -g $chrom_size -w $bin_size|\
bedtools coverage -a - -b $bed |cut -f 1-4|\
perl -slane 'if ($F[3]!=0){$a=sprintf("%.2f",$F[3]*$s);$"="\t";print "@F[0..2]\t$a"}' -- -s=$scale_factor > $tmp1

if [ -z $black_list ]
then
	bedtools intersect -a $tmp1 -b $black_list -v > $tmp2
	rm $tmp1
else
	tmp2=$tmp1
fi

if [ $bin_merge == "region" ]
then
	bedtools groupby -i $tmp2 -g 1,4 -c 2,3 -o min,max|\
	awk -v OFS='\t' '{print $1, $3, $4, $2}'|\
	sort -k1,1 -k2,2n -S 10G --parallel=10 > $bdg
else
	sort -k1,1 -k2,2n -S 10G --parallel=10 $tmp2 > $bdg
fi

rm $tmp2

#covert bedgraph to bigwig
bedGraphToBigWig $bdg $chrom_size $bw
echo "bed_to_bdg done!"


#-----------------------------------------------------------------------------------------------------------------------------






