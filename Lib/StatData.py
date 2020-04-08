#!/usr/bin/env python2


import sys
import re
import os
import json
import argparse
import ConfigParser
import pandas as pd
from glob import glob


BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../../Bin/config.ini')

def GetSequenceQC(file_in):
    with open(file_in, 'r') as f:
        content = f.read()
    data = json.loads(content)
    raw = data['summary']['before_filtering']['total_reads']
    qc = data['summary']['after_filtering']['total_reads']
    q20 = data['summary']['after_filtering']['q20_rate']
    q30 = data['summary']['after_filtering']['q30_rate']
    gc = data['summary']['after_filtering']['gc_content']
    return raw,qc,q20,q30,gc

def StatQC(path_in):
    RawData = path_in+'/RawData'
    QC = path_in+'/1.QC'
    list_sample = []
    for fl in glob(RawData+'/*_1.fq.gz'):
        lb = re.split('_1\.fq\.gz', os.path.basename(fl))[0]
        list_sample.append(lb)
    pd_stat = pd.DataFrame(index=sorted(list_sample), columns=['#RawReads', '#ReadsAfterQC', 'QCRate(%)', 'Q20','Q30','GC'])
    for fl in list_sample:
        raw,qc,q20,q30,gc = GetSequenceQC(os.path.join(QC, fl+'_QC_report.json'))
        QCRate = round((float(qc)*100)/raw, 2)
        pd_stat.loc[fl, :] = [raw,qc,QCRate,q20,q30,gc]
    return pd_stat

def main():
    parser = argparse.ArgumentParser(description="QC STAT pipeline")
    parser.add_argument('-i', help='the qc analysis main path', required=True)
    parser.add_argument('-o', help='the output file', default='./qcstat.xlsx')
    argv=vars(parser.parse_args())
    pd_qc = StatQC(argv['i'])
    qc_out = pd.ExcelWriter(argv['o'])
    pd_qc.loc[:, ['#RawReads', '#ReadsAfterQC', 'QCRate(%)', 'Q20', 'Q30', 'GC']].to_excel(qc_out, sheet_name='Sequencing QC', header=True, index=True)
    qc_out.close()


if __name__ == '__main__':
    main()
