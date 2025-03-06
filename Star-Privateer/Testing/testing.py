import star_privateer as sp
import numpy as np
import os, pathos
import matplotlib.pyplot as plt
from tqdm import tqdm


if not os.path.exists ('rooster_training_features') :
    os.mkdir ('rooster_training_features')
if not os.path.exists ('rooster_instances') :
    os.mkdir ('rooster_instances')
    
list_kic = sp.get_list_targets ("kepler_dataset")

def analysis_wrapper (kic) :
    """
    Analysis wrapper to speed computation
    by parallelising process and control
    memory usage.
    """
    str_kic = str (kic).zfill (9)
    filename = sp.get_target_filename ("kepler_dataset", str_kic)
    fileout = 'rooster_training_features/{}.csv'.format(str_kic)
    fileplot = 'rooster_training_features/{}.png'.format(str_kic)
    if not os.path.exists (fileout) :
        t, s, dt = sp.load_resource (filename)
        (p_ps, p_acf,
         ps, acf,
         cs, features,
         feature_names,
         fig) = sp.analysis_pipeline (t, s, pmin=0.1, pmax=60,
                                      wavelet_analysis=False, plot=True,
                                      filename=fileplot, figsize=(10,16),
                                      lw=1, dpi=150, pfa_threshold=1e-6,
                                      ls_err_smooth=True)
        df = sp.save_features (fileout, kic, features, feature_names)
        plt.close ("all")


process_pool = pathos.pools._ProcessPool (processes=4,
                                          maxtasksperchild=10)
with process_pool as p :
    list (tqdm (p.imap (analysis_wrapper,
                        list_kic,
                        ),
                total=len (list_kic))
          )
    p.close ()

df = sp.build_catalog_features ('rooster_training_features')

df 