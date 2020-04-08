#!/usr/bin/env python2
# coding=utf-8
import sys
import re
import os
import ConfigParser
import argparse

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../Bin/config.ini')

#### SOFT
PYTHON = config.get('SOFTWARE', 'python')
RSCRIPT = config.get('SOFTWARE', 'Rscript')
FASTP = config.get('SOFTWARE', 'fastp')
SNAKEMAKE = config.get('SOFTWARE', 'snakemake')

#### SCRIPT
STATDATA = config.get('SCRIPT', 'StatData')

#### DATABASE

class ReadList(object):
    def __init__(self, line_in):
        list_split = re.split('\t', line_in)
        self.Sample = list_split[0]
        self.Name = list_split[1]
        list_fq = re.split(',', list_split[2])
        self.fq1 = list_fq[0]
        self.fq2 = list_fq[1]
#        self.Genda = list_split[3]

class Snake(object):
    def __init__(self, process):
        self.process = process
        self.input = ''
        self.output = ''
        self.params = ''
        self.log = ''
        self.threads = ''
        self.shell = ''

    def UpdateInput(self, line_in):
        self.input = line_in

    def UpdateOutput(self, line_in):
        self.output = line_in

    def UpdateParams(self, line_in):
        self.params = line_in

    def UpdateLog(self, line_in):
        self.log = line_in

    def UpdateThreads(self, line_in):
        self.threads = line_in

    def UpdateShell(self, line_in):
        self.shell = line_in

    def WriteStr(self, fn):
        fn.write('rule '+self.process+':\n')
        fn.write('\tinput:\n\t\t'+self.input+'\n')
        if self.output:
            fn.write('\toutput:\n\t\t'+self.output+'\n')
        if self.params:
            fn.write('\tparams:\n\t\t'+self.params+'\n')
        if self.log:
            fn.write('\tlog:\n\t\t'+self.log+'\n')
        if self.threads:
            fn.write('\tthreads: '+self.threads+'\n')
        if self.shell:
            fn.write('\tshell:\n\t\t'+self.shell+'\n')
        fn.write('\n')

def MakeSnake(file_in, file_out, platform, cores):
    out = open(file_out, 'w')
    if platform == 'SGE':
        out.write(SNAKEMAKE+' --cluster "qsub -cwd -V -b y -S /bin/bash -o {log.o} -e {log.e}" -j '+cores+' -s '+file_in+' --printshellcmds --latency-wait 10')
    elif platform == 'local':
        out.write('{} -j {} -s {} --printshellcmds --latency-wait 10'\
                 .format(SNAKEMAKE, cores, file_in))
    out.close()

def main():
    parser = argparse.ArgumentParser(description="WES pipeline")
    parser.add_argument('-c', help='the input fasta list', required=True)
    parser.add_argument('-o', help='the output path', required=True)
    parser.add_argument('-p', help='the platform', choices=['SGE', 'local'], default='local')
    parser.add_argument('-j', help='the job parallel cores', default='2')
    parser.add_argument('-r', help='run now', action='store_true')
    argv=vars(parser.parse_args())
    outpath = argv['o']
    snakefile = open(os.path.join(outpath, 'snakefile.txt'), 'w')
    RawData = os.path.join(outpath, 'RawData')
    list_ob = []
    if not os.path.exists(RawData):
        os.mkdir(RawData)
    else:
        os.system('rm -rf '+RawData)
        os.mkdir(RawData)
    with open(argv['c'], 'r') as f:
        os.chdir(RawData)
        for line in f:
            if not line.startswith('#'):
                ob = ReadList(line.strip())
                list_ob.append(ob)
                os.system('ln -s {} {}'.format(ob.fq1, ob.Name+'_1.fq.gz'))
                os.system('ln -s {} {}'.format(ob.fq2, ob.Name+'_2.fq.gz'))
    os.chdir(outpath)
    ### config file
    snakefile.write('Samples = "{}".split()\n'.format(' '.join([i.Name for i in\
                                                              list_ob])))

    ###All
    All = Snake('All')
    All.UpdateInput('"6.Final/QCStat.xls"')
    All.WriteStr(snakefile)
    ### QC
    QC = Snake('QC')
    QC.UpdateInput('"RawData/{sample}_1.fq.gz"')
    QC.UpdateOutput('"1.QC/{sample}_1.clean.fq.gz"')
    QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
    QC.UpdateShell(r'"'+FASTP+r' -i RawData/{wildcards.sample}_1.fq.gz -o 1.QC/{wildcards.sample}_1.clean.fq.gz -I RawData/{wildcards.sample}_2.fq.gz -O 1.QC/{wildcards.sample}_2.clean.fq.gz -A -j 1.QC/{wildcards.sample}_QC_report.json -h 1.QC/{wildcards.sample}_QC_report.html"')
    BamIndex.UpdateOutput('"2.Align/{sample}.sort.dedup.bam.bai"')
    BamIndex.UpdateLog('e = "logs/{sample}.BamIndex.e", o = "logs/{sample}.BamIndex.o"')
    BamIndex.UpdateShell(r'"'+SAMTOOLS+r' index {input}"')
    BamIndex.WriteStr(snakefile)
    ### BaseRecalibrator
    BaseRecalibrator = Snake('BaseRecalibrator')
    BaseRecalibrator.UpdateInput('"2.Align/{sample}.sort.dedup.bam.bai"')
    BaseRecalibrator.UpdateOutput('"3.Call/{sample}.recalibration.report"')
    BaseRecalibrator.UpdateLog('e = "logs/{sample}.BaseRecalibrator.e", o = "logs/{sample}.BaseRecalibrator.o"')
    BaseRecalibrator.UpdateShell(r'"'+GATK+' BaseRecalibrator --input 2.Align/{wildcards.sample}.sort.dedup.bam --output {output} --reference {FASTA} --intervals {REGION} --known-sites {G1000} --known-sites {Mills}"')
    BaseRecalibrator.WriteStr(snakefile)
    ### GatherBQSRReports
    GatherBQSRReports = Snake('GatherBQSRReports')
    GatherBQSRReports.UpdateInput('expand("3.Call/{sample}.recalibration.report", sample=Samples)')
    insh = '-I '+' -I '.join(["3.Call/"+i.Name+".recalibration.report" for i in list_ob])
    GatherBQSRReports.UpdateOutput('"3.Call/all.recalibration.report"')
    GatherBQSRReports.UpdateLog('e = "logs/GatherBQSRReports.e", o = "logs/GatherBQSRReports.o"')
    GatherBQSRReports.UpdateShell(r'"'+GATK+' GatherBQSRReports '+insh+' -O {output}"')
    GatherBQSRReports.WriteStr(snakefile)
    ### ApplyBQSR
    ApplyBQSR = Snake('ApplyBQSR')
    ApplyBQSR.UpdateInput('report="3.Call/all.recalibration.report",bam="2.Align/{sample}.sort.dedup.bam"')
    ApplyBQSR.UpdateOutput('"3.Call/{sample}.sort.dedup.BQSR.bam"')
    ApplyBQSR.UpdateLog('e = "logs/{sample}.ApplyBQSR.e", o = "logs/{sample}.ApplyBQSR.o"')
    ApplyBQSR.UpdateShell(r'"'+GATK+' ApplyBQSR -R {FASTA} -I {input.bam} -O {output} -bqsr {input.report} --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --add-output-sam-program-record --use-original-qualities -L {REGION}"')
    ApplyBQSR.WriteStr(snakefile)
    ### HaplotypeCaller
    HaplotypeCaller = Snake('HaplotypeCaller')
    HaplotypeCaller.UpdateInput('"3.Call/{sample}.sort.dedup.BQSR.bam"')
    HaplotypeCaller.UpdateOutput('"3.Call/{sample}.gvcf"')
    HaplotypeCaller.UpdateLog('e = "logs/{sample}.HaplotypeCaller.e", o = "logs/{sample}.HaplotypeCaller.o"')
    HaplotypeCaller.UpdateShell(r'"'+GATK+' HaplotypeCaller -R {FASTA} --emit-ref-confidence GVCF -I {input} -D {DBSNP} -L {REGION} -O {output}"')
    HaplotypeCaller.WriteStr(snakefile)
    ### CombineGVCFs
    CombineGVCFs = Snake("CombineGVCFs")
    cmsh = '-V '+' -V '.join(["3.Call/"+i.Name+".gvcf" for i in list_ob])
    CombineGVCFs.UpdateInput('expand("3.Call/{sample}.gvcf", sample=Samples)')
    CombineGVCFs.UpdateOutput('"3.Call/All.gvcf"')
    CombineGVCFs.UpdateLog('e = "logs/CombineGVCFs.e", o = "logs/CombineGVCFs.o"')
    CombineGVCFs.UpdateShell(r'"'+GATK+' CombineGVCFs -R {FASTA} '+cmsh+' -O {output}"')
    CombineGVCFs.WriteStr(snakefile)
    ### GenotypeGVCFs
    GenotypeGVCFs = Snake("GenotypeGVCFs")
    GenotypeGVCFs.UpdateInput('"3.Call/All.gvcf"')
    GenotypeGVCFs.UpdateOutput('"3.Call/All.raw.gvcf"')
    GenotypeGVCFs.UpdateLog('e = "logs/GenotypeGVCFs.e", o = "logs/GenotypeGVCFs.o"')
    GenotypeGVCFs.UpdateShell(r'"'+GATK+' GenotypeGVCFs -R {FASTA} -V {input} -O {output} -L {REGION}"')
    GenotypeGVCFs.WriteStr(snakefile)
    ### SelectSNP
    SelectSNP = Snake("SelectSNP")
    SelectSNP.UpdateInput('"3.Call/All.raw.gvcf"')
    SelectSNP.UpdateOutput('"4.Mutation/All.raw.snp.vcf"')
    SelectSNP.UpdateLog('e = "logs/SelectSNP.e", o = "logs/SelectSNP.o"')
    SelectSNP.UpdateShell(r'"'+GATK+' SelectVariants -select-type SNP -V {input} -O {output}"')
    SelectSNP.WriteStr(snakefile)
    ### SelectIndel
    SelectIndel = Snake("SelectIndel")
    SelectIndel.UpdateInput('"3.Call/All.raw.gvcf"')
    SelectIndel.UpdateOutput('"4.Mutation/All.raw.indel.vcf"')
    SelectIndel.UpdateLog('e = "logs/SelectIndel.e", o = "logs/SelectIndel.o"')
    SelectIndel.UpdateShell(r'"'+GATK+' SelectVariants -select-type INDEL -V {input} -O {output}"')
    SelectIndel.WriteStr(snakefile)
    ### FilterSNP
    FilterSNP = Snake("FilterSNP")
    FilterSNP.UpdateInput('"4.Mutation/All.raw.snp.vcf"')
    FilterSNP.UpdateOutput('"4.Mutation/All.raw.snp.filter.vcf"')
    FilterSNP.UpdateLog('e = "logs/FilterSNP.e", o = "logs/FilterSNP.o"')
    FilterSNP.UpdateShell('\''+GATK+' VariantFiltration -V {input} --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O {output}'+'\'')
    FilterSNP.WriteStr(snakefile)
    ### FilterIndel
    FilterIndel = Snake("FilterIndel")
    FilterIndel.UpdateInput('"4.Mutation/All.raw.indel.vcf"')
    FilterIndel.UpdateOutput('"4.Mutation/All.raw.indel.filter.vcf"')
    FilterIndel.UpdateLog('e = "logs/FilterIndel.e", o = "logs/FilterIndel.o"')
    FilterIndel.UpdateShell('\''+GATK+' VariantFiltration -V {input} --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter" -O {output}'+'\'')
    FilterIndel.WriteStr(snakefile)
    ### MergeVcfs
    MergeVcfs = Snake("MergeVcfs")
    MergeVcfs.UpdateInput('indel="4.Mutation/All.raw.indel.filter.vcf",snp="4.Mutation/All.raw.snp.filter.vcf"')
    MergeVcfs.UpdateOutput('"4.Mutation/All.raw.filter.vcf"')
    MergeVcfs.UpdateLog('e = "logs/MergeVcfs.e", o = "logs/MergeVcfs.o"')
    MergeVcfs.UpdateShell(r'"'+GATK+' MergeVcfs -I {input.snp} -I {input.indel} -O {output}"')
    MergeVcfs.WriteStr(snakefile)
    ### Convert2Annovar
    Convert2Annovar = Snake("Convert2Annovar")
    Convert2Annovar.UpdateInput('"4.Mutation/All.raw.{tp}.filter.vcf"')
    Convert2Annovar.UpdateOutput('"4.Mutation/All.raw.{tp}.filter.Avinput"')
    Convert2Annovar.UpdateLog('e = "logs/{tp}.Convert2Annovar.e", o = "logs/{tp}.Convert2Annovar.o"')
    Convert2Annovar.UpdateShell(r'"'+ANNOVAR+'/convert2annovar.pl -format vcf4 -allsample -withfreq --filter PASS {input} > {output}"')
    Convert2Annovar.WriteStr(snakefile)
    ### TableAnnovar
    TableAnnovar = Snake("TableAnnovar")
    TableAnnovar.UpdateInput('"4.Mutation/All.raw.{tp}.filter.Avinput"')
    TableAnnovar.UpdateOutput('"5.Anno/1.Annovar/All.raw.{tp}.filter.anno.hg19_multianno.txt"')
    TableAnnovar.UpdateLog('e = "logs/{tp}.TableAnnovar.e", o = "logs/{tp}.TableAnnovar.o"')
    TableAnnovar.UpdateShell(r'"'+ANNOVAR+'/table_annovar.pl {input} '+ANNOVAR+'/humandb/ -buildver hg19 -out 5.Anno/1.Annovar/All.raw.{wildcards.tp}.filter.anno -remove -protocol refGene,phastConsElements46way,genomicSuperDups,exac03,snp138,esp6500siv2_all,1000g2015aug_all,clinvar_20190305,ljb26_all,cosmic68 -operation g,r,r,f,f,f,f,f,f,f -nastring . -polish"')
    TableAnnovar.WriteStr(snakefile)
    ### Intervar
    Intervar = Snake("Intervar")
    Intervar.UpdateInput('"4.Mutation/All.raw.{tp}.filter.Avinput"')
    Intervar.UpdateOutput('label="5.Anno/2.Intervar/All.raw.{tp}.filter.hg19_multianno.txt.intervar"')
    Intervar.UpdateLog('e = "logs/{tp}.Intervar.e", o = "logs/{tp}.Intervar.o"')
    Intervar.UpdateShell(r'"'+INTERVAR+'/Intervar.py -b hg19 -i {input} --input_type avinput -o 5.Anno/2.Intervar/All.raw.{wildcards.tp}.filter -t '+INTERVAR+'/intervardb --table_annovar '+ANNOVAR+'/table_annovar.pl --convert2annovar '+ANNOVAR+'/convert2annovar.pl --annotate_variation '+ANNOVAR+'/annotate_variation.pl -d '+ANNOVAR+'/humandb/"')
    Intervar.WriteStr(snakefile)
    ### MergeAnno
    MergeAnno = Snake("MergeAnno")
    MergeAnno.UpdateInput('annovar="5.Anno/1.Annovar/All.raw.{tp}.filter.anno.hg19_multianno.txt", interver="5.Anno/2.Intervar/All.raw.{tp}.filter.hg19_multianno.txt.intervar"')
    MergeAnno.UpdateOutput('"6.Final/All.{tp}.interver.annovar.txt"')
    MergeAnno.UpdateLog('e = "logs/{tp}.MergeAnno.e", o = "logs/{tp}.MergeAnno.o"')
    MergeAnno.UpdateShell(r'"'+PYTHON+' '+MERGEDATA+' -i {input.interver} -a {input.annovar} -o {output}"')
    MergeAnno.WriteStr(snakefile)
    ### AddAnno
    AddAnno = Snake('AddAnno')
    AddAnno.UpdateInput('"6.Final/All.{tp}.interver.annovar.txt"')
    AddAnno.UpdateOutput('"6.Final/All.{tp}.finalanno.txt"')
    AddAnno.UpdateLog('e = "logs/{tp}.AddAnno.e", o = "logs/{tp}.AddAnno.o"')
    AddAnno.UpdateShell(r'"'+PYTHON+' '+ADDANNO+' -i {input} -o {output}"')
    AddAnno.WriteStr(snakefile)
    ### FilterIntervar
    FilterIntervar = Snake('FilterIntervar')
    FilterIntervar.UpdateInput('"6.Final/All.{tp}.finalanno.txt"')
    FilterIntervar.UpdateOutput('"6.Final/Final{tp}Result.txt", "6.Final/Final{tp}ROI.txt"')
    FilterIntervar.UpdateLog('e = "logs/{tp}.FilterIntervar.e", o = "logs/{tp}.FilterIntervar.o"')
    FilterIntervar.UpdateShell(r'"'+PYTHON+' '+FILTERINTERVAR+' -i {input} -o 6.Final/"')
    FilterIntervar.WriteStr(snakefile)
    ### SepSample
    SepSample = Snake('SepSample')
    SepSample.UpdateInput('anno="6.Final/Final{tp}ROI.txt",vcf="4.Mutation/All.raw.{tp}.filter.vcf"')
    SepSample.UpdateOutput('"6.Final/FinalSample{tp}ROI.txt"')
    SepSample.UpdateLog('e = "logs/{tp}.SepSample.e", o = "logs/{tp}.SepSample.o"')
    SepSample.UpdateShell(r'"'+PYTHON+' '+SEPSAMPLEF+' -i {input.anno} -v {input.vcf} -o {output}"')
    SepSample.WriteStr(snakefile)
    ### SepSample1
    SepSample1 = Snake('SepSample1')
    SepSample1.UpdateInput('anno="6.Final/Final{tp}Result.txt",vcf="4.Mutation/All.raw.{tp}.filter.vcf"')
    SepSample1.UpdateOutput('"6.Final/FinalSample{tp}.txt"')
    SepSample1.UpdateLog('e = "logs/{tp}.SepSample1.e", o = "logs/{tp}.SepSample1.o"')
    SepSample1.UpdateShell(r'"'+PYTHON+' '+SEPSAMPLEF+' -i {input.anno} -v {input.vcf} -o {output}"')
    SepSample1.WriteStr(snakefile)
    ### SepSample2
    SepSample2 = Snake('SepSample2')
    SepSample2.UpdateInput('anno="6.Final/All.{tp}.finalanno.txt",vcf="4.Mutation/All.raw.{tp}.filter.vcf"')
    SepSample2.UpdateOutput('"6.Final/AllSample{tp}.txt"')
    SepSample2.UpdateLog('e = "logs/{tp}.SepSample2.e", o = "logs/{tp}.SepSample2.o"')
    SepSample2.UpdateShell(r'"'+PYTHON+' '+SEPSAMPLEF+' -i {input.anno} -v {input.vcf} -o {output}"')
    SepSample2.WriteStr(snakefile)
    ### IntergrateMutation
    IntergrateMutation = Snake('IntergrateMutation')
    IntergrateMutation.UpdateInput('snp="6.Final/{cls}Samplesnp.txt",indel="6.Final/{cls}Sampleindel.txt"')
    IntergrateMutation.UpdateOutput('"6.Final/{cls}IntergrateMutation.txt"')
    IntergrateMutation.UpdateLog('e = "logs/Final.IntergrateMutation.e", o = "logs/Final.IntergrateMutation.o"')
    IntergrateMutation.UpdateShell(r'"'+PYTHON+' '+INTERGRATEMUTATION+' -s {input.snp} -i {input.indel} -o {output}"')
    IntergrateMutation.WriteStr(snakefile)
    ### StatQC
    StatQC = Snake('StatQC')
    StatQC.UpdateInput('a="6.Final/AllIntergrateMutation.txt", b="6.Final/FinalIntergrateMutation.txt"')
    StatQC.UpdateOutput('"6.Final/QCStat.xls"')
    StatQC.UpdateLog('e = "logs/StatQC.e", o = "logs/StatQC.o"')
    StatQC.UpdateShell(r'"'+PYTHON+' '+STATDATA+' -i '+outpath+' -b {REGION} -o {output}"')
    StatQC.WriteStr(snakefile)

    ######RUN
    OutShell = os.path.join(outpath, 'work.sh')
    MakeSnake(os.path.join(outpath, 'snakefile.txt'), OutShell, argv['p'], argv['j'])
    if argv['r']:
        os.system('nohup sh {}&'.format(OutShell))

if __name__ == '__main__':
    main()
