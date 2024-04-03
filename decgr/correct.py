#!/usr/bin/env python


## Required modules
import os, subprocess
from decgr.cnv.loadcnv import binCNV
from decgr.cnv.correctcnv import matrix_balance
from decgr.cnv import runcnv
import cooler
from decgr.cnv.segcnv import HMMsegment
import numpy as np
from rpy2.robjects import numpy2ri, Formula
from rpy2.robjects.packages import importr


def run(hic, genome, enzyme, resolution):

    cnv_file = "./result/CNV_profile.txt"
    cnv_seg_file = "./result/CNV_seg_profile.txt"
    cachefolder = "./result"

    # calculate CNV profile
    weblinks = {
        'hg38_mappability_100mer.1kb.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38_mappability_100mer.1kb.bw',
        'hg38.MboI.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.MboI.npz',
        'hg38.DpnII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.MboI.npz',
        'hg38.Arima.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.Arima.npz',
        'hg38.BglII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.BglII.npz',
        'hg38.uniform.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.uniform.npz',
        'hg38.HindIII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.HindIII.npz',
        'hg38_1kb_GC.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38_1kb_GC.bw',
        'hg19_mappability_100mer.1kb.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19_mappability_100mer.1kb.bw',
        'hg19.MboI.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.MboI.npz',
        'hg19.Arima.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.Arima.npz',
        'hg19.HindIII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.HindIII.npz',
        'hg19.DpnII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.MboI.npz',
        'hg19.BglII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.BglII.npz',
        'hg19.uniform.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.uniform.npz',
        'hg19_1kb_GC.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19_1kb_GC.bw',
        'mm10_mappability_100mer.1kb.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10_mappability_100mer.1kb.bw',
        'mm10_1kb_GC.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10_1kb_GC.bw',
        'mm10.Arima.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.Arima.npz',
        'mm10.BglII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.BglII.npz',
        'mm10.DpnII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.DpnII.npz',
        'mm10.HindIII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.HindIII.npz',
        'mm10.MboI.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.MboI.npz',
        'mm10.uniform.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.uniform.npz',
        'mm9_mappability_100mer.1kb.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9_mappability_100mer.1kb.bw',
        'mm9_1kb_GC.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9_1kb_GC.bw',
        'mm9.Arima.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.Arima.npz',
        'mm9.BglII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.BglII.npz',
        'mm9.DpnII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.DpnII.npz',
        'mm9.HindIII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.HindIII.npz',
        'mm9.MboI.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.MboI.npz',
        'mm9.uniform.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.uniform.npz'
    }

    cachefolder = os.path.abspath(os.path.expanduser(cachefolder))
    if not os.path.exists(cachefolder):
        os.makedirs(cachefolder)
    mapscore_fil = os.path.join(cachefolder, '{0}_mappability_100mer.1kb.bw'.format(genome))
    cutsites_fil = os.path.join(cachefolder, '{0}.{1}.npz'.format(genome, enzyme))
    gc_fil = os.path.join(cachefolder, '{0}_1kb_GC.bw'.format(genome))
    for fil in [mapscore_fil, cutsites_fil, gc_fil]:
        if not os.path.exists(fil):
            key = os.path.split(fil)[1]
            command = ['wget', '-O', fil, '-L', weblinks[key]]
            subprocess.check_call(' '.join(command), shell=True)
    mgcv = importr('mgcv')
    stats = importr('stats')
    table, res = runcnv.get_marginals(hic)
    table = runcnv.signal_from_bigwig(table, gc_fil, name='GC')
    table = runcnv.signal_from_bigwig(table, mapscore_fil, name='Mappability')
    table = runcnv.count_REsites(table, cutsites_fil, res)
    mask, filtered = runcnv.filterZeros(table)
    print("running cnv")
    fomula = Formula('Coverage ~ s(GC) + s(Mappability) + s(RE)')
    fomula.environment['Coverage'] = numpy2ri.numpy2rpy(filtered['Coverage'].values)
    fomula.environment['GC'] = numpy2ri.numpy2rpy(filtered['GC'].values)
    fomula.environment['Mappability'] = numpy2ri.numpy2rpy(filtered['Mappability'].values)
    fomula.environment['RE'] = numpy2ri.numpy2rpy(filtered['RE'].values)
    gam = mgcv.gam(fomula, family=stats.poisson(link='log'))
    rs = mgcv.residuals_gam(gam, type='working')
    residuals = numpy2ri.rpy2py(rs)
    residuals = residuals - residuals.min()
    idx = np.where(mask)[0]
    CNV = np.zeros(table.shape[0])
    CNV[idx] = residuals
    table['CNV'] = CNV
    bedgraph = table[['chrom', 'start', 'end', 'CNV']]
    bedgraph.to_csv(cnv_file, sep='\t', header=False, index=False)
    print("segmenting cnv")
    # segment CNV profile
    work = HMMsegment(cnv_file,
                      res=resolution,
                      nproc=1,
                      ploidy=2,
                      n_states=None)

    # logger.info('Perform segmentation ...')
    work.segment(
        min_seg=3, min_diff=0.4, max_dist=4, p=1e-5
    )
    work.output(cnv_seg_file)
    print("correcting cnv")
    # correct CNV profile
    hic_pool = cooler.Cooler(hic)
    bincnv = binCNV(cnv_seg_file, hic_pool.binsize)
    bincnv.assign_cnv(hic)
    matrix_balance(hic, nproc=12, mad_max=5, min_nnz=10, ignore_diags=1)
