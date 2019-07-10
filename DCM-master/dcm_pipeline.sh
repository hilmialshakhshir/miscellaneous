#!/bin/sh

usage() {
  echo "Usage"
  echo " -id, --inputDir	    Input directory containing BAM files"
  echo " -rd, --referenceDir        Directory containing reference files"
  echo " -od, --outputDir	    Output directory containing ASE results"
  echo " -td, --tempDir             Temporary directory for intermediary results"
  echo " Note: these should all be absolute paths to the directories, and the directories should not contain links!"
  echo " -p   --numberOfProcessors  Number of processors to use."
}

until [ $# -eq 0 ]; do
  case "$1" in
    -id  | --inputDir ) 		ID=$2; shift 2 ;;
    -rd  | --referenceDir )             RD=$2; shift 2 ;;
    -od  | --outputDir )		OD=$2; shift 2 ;;
    -td  | --tempDir )                  TD=$2; shift 2 ;;
    -p   | --numberOfProcessors )        P=$2; shift 2 ;;
    * ) shift ;;
  esac
done

V_2BIT=hg38 # Fix for the fact that deeptools correctGCBias does not accept wildcard.
VCF=1000G_phase1.snps.high_confidence.hg38.vcf.gz
FASTA=GRCh38.fa

# TODO: Input dir is...
echo $ID
echo $RD
echo $OD
echo $TD
echo $P
echo '! '$V_2BIT
echo '! '$VCF
echo '! '$FASTA

# TODO: Add P
if [[ -z "$ID" || -z "$RD" || -z "$OD" || -z "$TD" ]]; then
  usage
  exit
fi

for f in $ID/*.bam ; do
  FILE=`basename ${f}`
  TMPFILE=`basename ${f%%.*}`

  
CMD1="docker run -i --rm -v ${ID}:/input -v ${RD}:/reference -v ${TD}:/tempdir gatk4-tool \
     gatk AddOrReplaceReadGroups -I=/input/${FILE} -O=/tempdir/${TMPFILE}.rg.bam \
       -RGID=${TMPFILE} -RGLB=DCM -RGPL=illumina -RGPU=dcm -RGSM=${TMPFILE} -CREATE_INDEX=true"

# TODO: Test incl. -p
CMD2="docker run -i --rm -v ${ID}:/input -v ${RD}:/reference -v ${TD}:/tempdir \
       deeptools computeGCBias -b /tempdir/${TMPFILE}.rg.bam --effectiveGenomeSize 3049315783 \
       --genome /reference/${V_2BIT}.2bit -l 100 --GCbiasFrequenciesFile /tempdir/${TMPFILE}.txt -p $P"

CMD3="docker run -i --rm -v ${ID}:/input -v ${OD}:/output -v ${RD}:/reference -v ${TD}:/tempdir \
       deeptools correctGCBias -b /tempdir/${TMPFILE}.rg.bam --effectiveGenomeSize 3049315783 \
       --genome /reference/${V_2BIT}.2bit --GCbiasFrequenciesFile /tempdir/${TMPFILE}.txt \
       -o /output/${TMPFILE}.gc.bam -p $P"

CMD4="docker run -i --rm -v ${ID}:/input -v ${OD}:/output -v ${RD}:/reference -v ${TD}:/tempdir \
       gatk4-tool gatk ASEReadCounter -I /output/${TMPFILE}.gc.bam -O /output/${TMPFILE}.csv \
        -V /reference/$VCF -R /reference/$FASTA -min-depth 10 -mmq 20 -mbq 5"



  
  echo -e "$CMD1\n"; eval $CMD1 
  echo 
  echo -e "$CMD2\n"; eval $CMD2
  echo 
  echo -e "$CMD3\n"; eval $CMD3
  echo 
  echo -e "$CMD4\n"; eval $CMD4
  echo 
  echo "removing all files in temp folder"
  
#  rm -f $TD/*
done
