#!/usr/bin/env python

## Required modules
import os, subprocess, sys, h5py
path = os.path.split(os.path.abspath(__file__))[0]
sys.path.append(path)
from cnv.loadcnv import binCNV
from cnv.correctcnv import matrix_balance
from cnv import runcnv
import cooler
from cnv.segcnv import HMMsegment
import numpy as np
from rpy2.robjects import numpy2ri, Formula
from rpy2.robjects.packages import importr
import wget

def run(hic, genome, enzyme, resolution, nproc, chromosomes, path):
    print(path)
    cachefolder = path + "/result"
    cnv_norm_data_folder = path + "/cnv_norm_data"

    # calculate CNV profile
    weblinks = {
        'hg38_mappability_100mer.1kb.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38_mappability_100mer.1kb.bw',
        'hg38.MboI.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.MboI.npz',
        'hg38.DpnII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.MboI.npz',
        'hg38.Arima.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.Arima.npz',
        'hg38.BglII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.BglII.npz',
        'hg38.uniform.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.uniform.npz',
        'hg38.HindIII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg38.HindIII.npz',
        'hg38_mappability_100mer.1kb.txt.gz': 'https://raw.githubusercontent.com/GaoLabXDU/DeCGR/main/cnv_norm_data/hg38_mappability_100mer.1kb.txt.gz',
        'hg38_1kb_GC.txt.gz': 'https://raw.githubusercontent.com/GaoLabXDU/DeCGR/main/cnv_norm_data/hg38_1kb_GC.txt.gz',
        'hg19.MboI.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.MboI.npz',
        'hg19.Arima.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.Arima.npz',
        'hg19.HindIII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.HindIII.npz',
        'hg19.DpnII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.MboI.npz',
        'hg19.BglII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.BglII.npz',
        'hg19.uniform.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/hg19.uniform.npz',
        'hg19_mappability_100mer.1kb.txt.gz': 'https://raw.githubusercontent.com/GaoLabXDU/DeCGR/main/cnv_norm_data/hg19_mappability_100mer.1kb.txt.gz',
        'hg19_1kb_GC.txt.gz': 'https://raw.githubusercontent.com/GaoLabXDU/DeCGR/main/cnv_norm_data/hg19_1kb_GC.txt.gz',
        'mm10.Arima.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.Arima.npz',
        'mm10.BglII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.BglII.npz',
        'mm10.DpnII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.DpnII.npz',
        'mm10.HindIII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.HindIII.npz',
        'mm10.MboI.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.MboI.npz',
        'mm10.uniform.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm10.uniform.npz',
        'mm10_mappability_100mer.1kb.txt.gz': 'https://raw.githubusercontent.com/GaoLabXDU/DeCGR/main/cnv_norm_data/mm10_mappability_100mer.1kb.txt.gz',
        'mm10_1kb_GC.txt.gz': 'https://raw.githubusercontent.com/GaoLabXDU/DeCGR/main/cnv_norm_data/mm10_1kb_GC.txt.gz',
        'mm9_mappability_100mer.1kb.bw': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9_mappability_100mer.1kb.bw',
        'mm9.Arima.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.Arima.npz',
        'mm9.BglII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.BglII.npz',
        'mm9.DpnII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.DpnII.npz',
        'mm9.HindIII.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.HindIII.npz',
        'mm9.MboI.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.MboI.npz',
        'mm9.uniform.npz': 'http://3dgenome.fsm.northwestern.edu/neoLoopFinder/mm9.uniform.npz',
        'mm9_mappability_100mer.1kb.txt.gz': 'https://raw.githubusercontent.com/GaoLabXDU/DeCGR/main/cnv_norm_data/mm9_mappability_100mer.1kb.txt.gz',
        'mm9_1kb_GC.txt.gz': 'https://raw.githubusercontent.com/GaoLabXDU/DeCGR/main/cnv_norm_data/mm9_1kb_GC.txt.gz'
    }
    cachefolder = os.path.abspath(os.path.expanduser(cachefolder))
    if not os.path.exists(cachefolder):
        os.makedirs(cachefolder)
    cnv_file = path + "/result/CNV_profile.txt"
    cnv_seg_file = path + "/result/CNV_seg_profile.txt"
    mapscore_fil = os.path.join(cnv_norm_data_folder, '{0}_mappability_100mer.1kb.txt.gz'.format(genome))
    cutsites_fil = os.path.join(cnv_norm_data_folder, '{0}.{1}.npz'.format(genome, enzyme))
    gc_fil = os.path.join(cnv_norm_data_folder, '{0}_1kb_GC.txt.gz'.format(genome))
    for fil in [mapscore_fil, cutsites_fil, gc_fil]:
        if not os.path.exists(fil):
            key = os.path.split(fil)[1]
            try:
                # Try downloading with wget
                wget.download(weblinks[key], fil)
            except Exception as e:
                return None
            

    mgcv = importr('mgcv')
    stats = importr('stats')
    table, res = runcnv.get_marginals(hic)
    table = runcnv.signal_from_txt(table, gc_fil, name='GC')
    table = runcnv.signal_from_txt(table, mapscore_fil, name='Mappability')
    table = runcnv.count_REsites(table, cutsites_fil, res)
    mask, filtered = runcnv.filterZeros(table)
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
    clr = cooler.Cooler(hic)
    bins = clr.bins()[:]
    selected_bins = bins[bins['chrom'].isin(chromosomes)].reset_index()
    bin_id_mapping = {old_id: new_id for new_id, old_id in zip(selected_bins.index, selected_bins['index'])}
    pixels = clr.pixels()[:]
    selected_pixels = pixels[pixels['bin1_id'].isin(bin_id_mapping.keys()) & pixels['bin2_id'].isin(bin_id_mapping.keys())]
    selected_pixels.loc[:, 'bin1_id'] = selected_pixels['bin1_id'].map(bin_id_mapping)
    selected_pixels.loc[:, 'bin2_id'] = selected_pixels['bin2_id'].map(bin_id_mapping)
    selected_pixels = selected_pixels.drop_duplicates(subset=['bin1_id', 'bin2_id'])
    new_cool_path = path + "/tmp/test.cool"
    cooler.create_cooler(new_cool_path, selected_bins, selected_pixels)
    hic = new_cool_path
    hic_pool = cooler.Cooler(hic)
    bincnv = binCNV(cnv_seg_file, hic_pool.binsize)
    bincnv.assign_cnv(hic)
    matrix_balance(hic, nproc=nproc, mad_max=5, min_nnz=10, ignore_diags=1)
    return new_cool_path