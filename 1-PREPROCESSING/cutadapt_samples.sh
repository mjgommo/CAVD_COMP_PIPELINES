accounting="LAB_JLP"
refpf=""
cutpf="cutadapt"
outdir="/data3/Project_BeatrizMartinez_RNASeqIntegration/Processed.Reads.2014"
mkdir -p $sgelogsdir

readsdir="$outdir/../Raw.Reads.AOValvs.2014"

sgelogsdir="$outdir/sge_logs"
mkdir -p $sgelogsdir

minlength=30
minoverlap=7
adaptor='AR00#=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
adaptorb='AR00#b=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
adaptor2='AR000=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
adaptor2b='AR000b=GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
trimlength=15

echo "cutadapt version"
echo `cutadapt --version`

tail='.gz'

ghjid="-h"

	 

echo "Adaptors1:$adaptor $adaptorb"
date +"%F %X"
echo ""

for f in `ls $readsdir/*$tail`
do

	pf=$(basename ${f%$tail})
	pf=${pf%.fq}
	pf=${pf%.fastq}
	sname=${pf%__*}
	sname=${sname##*.}
	echo "$sname $pf"

	hjid=$ghjid

	jid=`echo "fastqc $f" \
	| qsub -terse -P prod -N "Fastqc_$sname" -l h_vmem=10G -v HOLD_PIPELINE=1 $hjid -A $accounting \
	-o $sgelogsdir/Fastqc_$sname.log -e $sgelogsdir/Fastqc_$sname.err`
	echo "JID:$jid. Name:Fastqc_${sname}. COMMAND:fastqc $f"
#	hjid="-hold_jid $jid"
	
	jid=`echo "cutadapt  -m $minlength --match-read-wildcards -O$minoverlap -o $outdir/$pf.$cutpf.fq.gz --too-short-output=$outdir/$pf.$cutpf.too_short.fq.gz -a $adaptor -a $adaptorb $f > $outdir/$pf.$cutpf.log 2>&1" \
	| qsub -P prod -terse -N "Cutadapt_$sname" -l h_vmem=3G -v HOLD_PIPELINE=1 $hjid -A $accounting \
	-o $sgelogsdir/Cutadapt_$sname.log -e $sgelogsdir/Cutadapt_$sname.err`
	echo "JID:$jid. Name:Cutadapt_${sname}. COMMAND:cutadapt -m $minlength --match-read-wildcards -O$minoverlap -o $outdir/$pf.$cutpf.fq.gz --too-short-output=$outdir/$pf.$cutpf.too_short.fq.gz -a $adaptor -a $adaptorb $f > $outdir/$pf.$cutpf.log 2>&1"
	hjid="-hold_jid $jid"

	jid=`echo "fastqc $outdir/$pf.$cutpf.fq.gz" \
	| qsub -terse -P prod -N "Fastqc_$sname" -l h_vmem=10G -v HOLD_PIPELINE=1 $hjid -A $accounting \
	-o $sgelogsdir/Fastqc_$sname.log -e $sgelogsdir/Fastqc_$sname.err`
	echo "JID:$jid. Name:Fastqc_${sname}. COMMAND:fastqc $outdir/$pf.$cutpf.fq.gz"
	hjid="-hold_jid $jid"

done

