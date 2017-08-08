#### create bigWigs to view in genome browser

module load python/2.7.11
for I in /well/got2d/rna-seq/data/STEMBANCC/BiPSC/merged/bam/*.bam; do
 base=`basename $I .bam`
 bamCoverage -b $I -o bigWigs/$base.bigWig -of bigwig -v --normalizeUsingRPKM
done

## now create one bigwig per group :

## DE - BiPSC 
bigWigMerge *DEBIPSC*.bigWig DE_BiPSC.bedGraph 

## DE - FiPSC
bigWigMerge *DEFIPSC*.bigWig DE_FiPSC.bedGraph 

## PP - BiPSC
bigWigMerge *PPBIPSC*.bigWig PP_BiPSC.bedGraph 

## PP - FiPSC
bigWigMerge *PPFIPSC*.bigWig PP_FiPSC.bedGraph 

## only take "normal" chromosomes 
grep ^chr DE_BiPSC.bedGraph > DE_BiPSC.chr.bedGraph 
grep ^chr DE_FiPSC.bedGraph > DE_FiPSC.chr.bedGraph 

## sort:

/well/got2d/agata/bin/ucsctools/bedSort DE_BiPSC.chr.bedGraph DE_BiPSC.chr.bedGraph
/well/got2d/agata/bin/ucsctools/bedSort DE_FiPSC.chr.bedGraph DE_FiPSC.chr.bedGraph

/well/got2d/agata/bin/ucsctools/bedGraphToBigWig DE_BiPSC.chr.bedGraph /well/got2d/agata/bin/NGseqBasic/conf/UCSCgenomeSizes/hg19.chrom.sizes DE_BiPSC.bw
/well/got2d/agata/bin/ucsctools/bedGraphToBigWig DE_FiPSC.chr.bedGraph /well/got2d/agata/bin/NGseqBasic/conf/UCSCgenomeSizes/hg19.chrom.sizes DE_FiPSC.bw

### cp to public html folder:
cp DE_BiPSC.bw /home/agata/public_html

### then these can be uploaded to the genome browser with:
track type=bigWig name="DE_BiPSC" description="BiPSC RNA-seq at DE stage" bigDataUrl=http://www.well.ox.ac.uk/~agata/DE_BiPSC.bw
track type=bigWig name="DE_FiPSC" description="FiPSC RNA-seq at DE stage" bigDataUrl=http://www.well.ox.ac.uk/~agata/DE_FiPSC3.bw


