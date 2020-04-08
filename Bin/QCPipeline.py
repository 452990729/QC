#!/usr/bin/env python2
# coding=utf-8
import sys
import re
import os
import ConfigParser
import argparse

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/config.ini')

def HandleConfig(path_in):
    if path_in.startswith('..'):
        return os.path.join(BasePath, path_in)
    else:
        return path_in

#### SOFT
PYTHON = HandleConfig(config.get('SOFTWARE', 'python'))
FASTP = HandleConfig(config.get('SOFTWARE', 'fastp'))
SNAKEMAKE = HandleConfig(config.get('SOFTWARE', 'snakemake'))

#### SCRIPT
STATDATA = HandleConfig(config.get('SCRIPT', 'StatData'))

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
    parser.add_argument('-c', help='the input fastq list', required=True)
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
    All.UpdateInput('"2.Final/QCStat.xlsx"')
    All.WriteStr(snakefile)
    ### QC
    QC = Snake('QC')
    QC.UpdateInput('a = "RawData/{sample}_1.fq.gz", b = "RawData/{sample}_2.fq.gz"')
    QC.UpdateOutput('a = "1.QC/{sample}_1.clean.fq.gz", b = "1.QC/{sample}_2.clean.fq.gz", c = "1.QC/{sample}_QC_report.json", d = "1.QC/{sample}_QC_report.html"')
    QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
    QC.UpdateShell(r'"'+FASTP+r' -i {input.a} -o {output.a} -I {input.b} -O {output.b} -A -j {output.c} -h {output.d}"')
    QC.WriteStr(snakefile)
    ### StatQC
    StatQC = Snake('StatQC')
    StatQC.UpdateInput('expand("1.QC/{sample}_QC_report.json", sample=Samples)')
    StatQC.UpdateOutput('"2.Final/QCStat.xlsx"')
    StatQC.UpdateLog('e = "logs/StatQC.e", o = "logs/StatQC.o"')
    StatQC.UpdateShell(r'"'+PYTHON+' '+STATDATA+' -i '+outpath+' -o {output}"')
    StatQC.WriteStr(snakefile)

    ######RUN
    OutShell = os.path.join(outpath, 'work.sh')
    MakeSnake(os.path.join(outpath, 'snakefile.txt'), OutShell, argv['p'], argv['j'])
    if argv['r']:
        os.system('nohup sh {}&'.format(OutShell))

if __name__ == '__main__':
    main()
