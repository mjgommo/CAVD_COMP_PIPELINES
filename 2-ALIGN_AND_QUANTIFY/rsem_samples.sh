accounting="LAB_JLP"
refpf="GRCh38.78"
rsempf="rsem"
gtf="/data3/genomes/Homo_sapiens/rsem/Homo_sapiens.GRCh38.78.wBioTypes.primary_assembly.rsem"
rsemdir="/data3/Project_BeatrizMartinez_RNASeqIntegration/Aligned.Reads.2012"
mkdir -p $rsemdir
readsdir="$rsemdir/../Processed.Reads.2012"
sgelogsdir="$rsemdir/sge_logs"
tail='cutadapt.fq.gz'

ths=6

echo "Reference:$gtf"
date +"%F %X"
echo ""
allsamplesgenes=""
allsamplestrans=""
mkdir -p $sgelogsdir
holdalljid=""

ghjid="-h"

#for lsname in `ls $readsdir/*$tail | grep -oP "\.C[\w]+\." | sort | uniq`
for lsname in `ls $readsdir/*$tail`
do
  sname=$(basename ${lsname%.$tail})
  sname=${sname#\/}
  sname=${sname#.}
  sname=${sname%.}
  sname=${sname%__}
  
  allsamplesgenes=$allsamplesgenes" $rsemdir/$sname.$refpf.$rsempf.genes.results"
  allsamplestrans=$allsamplestrans" $rsemdir/$sname.$refpf.$rsempf.isoforms.results"

  if [ ! -e $rsemdir/$sname.$refpf.$rsempf.genes.results ]
   then

        lfq=$rsemdir/$sname".fq"

	hjid=$ghjid

	echo "Sample:$sname   File:$lfq"
	jid=`echo "zcat $readsdir/*$sname*$tail > $lfq" \
	| qsub -P prod -terse -N "GetFq_$sname" -v HOLD_PIPELINE=1 $hjid -A $accounting \
	-o $sgelogsdir/GetFq_$sname.log -e $sgelogsdir/GetFq_$sname.err`
	echo "JID:$jid. Name:GetFq_$sname. COMMAND:zcat $readsdir/*$sname*$tail > $lfq"
        ls $readsdir/*$sname*$tail
#	hjid="-hold_jid $jid"
	
        jid=`echo "rsem-calculate-expression -p $ths --output-genome-bam --sampling-for-bam --bowtie-e 60 --bowtie-m 30 --bowtie-chunkmbs 512 --fragment-length-mean 180 --fragment-length-sd 50 --estimate-rspd --ci-memory 2048 $lfq $gtf $rsemdir/$sname.$refpf.$rsempf >> $rsemdir/$sname.$refpf.$rsempf.log 2>&1" \
	| qsub -P prod -terse -N "RsemExp_$sname" -pe pthreads $ths -l h_vmem=5G -v HOLD_PIPELINE=1 $hjid -A $accounting \
	-o $sgelogsdir/RsemExp_$sname.log -e $sgelogsdir/RsemExp_$sname.err`
        echo "JID:$jid. Name:RsemExp_$sname. COMMAND:rsem-calculate-expression -p $ths --output-genome-bam --sampling-for-bam --bowtie-e 60 --bowtie-m 30 --bowtie-chunkmbs 512 --fragment-length-mean 180 --fragment-length-sd 50 --estimate-rspd --ci-memory 2048 $lfq $gtf $rsemdir/$sname.$refpf.$rsempf"
	hjid="-hold_jid $jid"
	
	jid=`echo "java -Xmx15G -jar $JAR_HOME/MarkDuplicates.jar I=$rsemdir/$sname.$refpf.$rsempf.genome.sorted.bam O=$rsemdir/$sname.$refpf.$rsempf.genome.sorted.markDup.bam M=$rsemdir/$sname.$refpf.$rsempf.genome.sorted.markDup.stats AS=true VALIDATION_STRINGENCY=LENIENT >> $rsemdir/$sname.$refpf.$rsempf.genome.sorted.markDup.log 2>&1" \
	| qsub -P prod -terse -N "MarkDup_$sname" -l h_vmem=20G -v HOLD_PIPELINE=1 $hjid -A $accounting \
	-o $sgelogsdir/MarkDup_$sname.log -e $sgelogsdir/MarkDup_$sname.err`
	echo "JID:$jid. Name:MarkDup_$sname. COMMAND:java -jar $JAR_HOME/MarkDuplicates.jar I=$rsemdir/$sname.$refpf.$rsempf.genome.sorted.bam O=$rsemdir/$sname.$refpf.$rsempf.genome.sorted.markDup.bam M=$rsemdir/$sname.$refpf.$rsempf.genome.sorted.markDup.stats AS=true VALIDATION_STRINGENCY=LENIENT"
	hjid="-hold_jid $jid"
	
	jid=`echo "samtools index $rsemdir/$sname.$refpf.$rsempf.genome.sorted.markDup.bam" \
	| qsub -P prod -terse -N "Index_$sname" -v HOLD_PIPELINE=1 $hjid -A $accounting \
	-o $sgelogsdir/Index_$sname.log -e $sgelogsdir/Index_$sname.err`
	echo "JID:$jid. Name:Index_$sname. COMMAND:samtools index $rsemdir/$sname.$refpf.$rsempf.genome.sorted.markDup.bam"
	hjid="-hold_jid $jid"
	
	jid=`echo "rm $rsemdir/$sname.$refpf.$rsempf.transcript.bam $rsemdir/$sname.$refpf.$rsempf.genome.bam $rsemdir/$sname.$refpf.$rsempf.genome.sorted.bam  $rsemdir/$sname.$refpf.$rsempf.genome.sorted.bam.bai $lfq" \
	| qsub -P prod -terse -N "Clean_$sname" -v HOLD_PIPELINE=1 $hjid -A $accounting \
	-o $sgelogsdir/Clean_$sname.log -e $sgelogsdir/Clean_$sname.err`
	echo "JID:$jid. Name:Clean_$sname. COMMAND:rm $rsemdir/$sname.$refpf.$rsempf.transcript.bam $rsemdir/$sname.$refpf.$rsempf.genome.bam $rsemdir/$sname.$refpf.$rsempf.genome.sorted.bam $rsemdir/$sname.$refpf.$rsempf.genome.sorted.bam.bai $lfq"
	hjid="-hold_jid $jid"
	
	holdalljid=$holdalljid","$jid
	echo ""
   else
	echo "EXIST: $rsemdir/$sname.$refpf.$rsempf.genes.results"
	echo ""
  fi
done

echo ""
holdalljid=${holdalljid#,}

hjid="-hold_jid $holdalljid"
echo $hjid

jid=`echo "rsem-generate-data-matrix $allsamplesgenes > $rsemdir/ALL2.$refpf.$rsempf.genes.results.txt; rsem-generate-data-matrix $allsamplestrans > $rsemdir/ALL.$refpf.$rsempf.isoforms.results.txt" \
| qsub -P prod -terse -N "RsMatrix_ALL" -v HOLD_PIPELINE=1 $hjid -A $accounting \
-o $sgelogsdir/RsMatrix_ALL.log -e $sgelogsdir/RsMatrix_ALL.err`
echo "JID:$jid. Name:RsMatrix_ALL. COMMAND:rsem-generate-data-matrix $allsamplesgenes > ALL.$refpf.$rsempf.genes.results.txt; rsem-generate-data-matrix $allsamplestrans > $rsemdir/ALL.$refpf.$rsempf.isoforms.results.txt"
echo ""

date +"%F %X"
echo "DONE"
