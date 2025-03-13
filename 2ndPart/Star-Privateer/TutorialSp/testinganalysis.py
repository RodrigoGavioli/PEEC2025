import star_privateer as sp
import os, pathos
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

if not os.path.exists ('stellar_analysis_features') :
    os.mkdir ('stellar_analysis_features')

filename = sp.get_target_filename (sp.timeseries, '003733735')
t, s, dt = sp.load_resource (filename)

(p_ps, p_acf, ps, acf,
 cs, features, feature_names, _) = sp.analysis_pipeline (t, s, pmin=0.1, pmax=60, figsize=(8,10),
                                                            wavelet_analysis=False, plot=True,
                                                            filename='stellar_analysis_features/003733735.png',
                                                        )



fileout = 'stellar_analysis_features/003733735.csv'
df = sp.save_features (fileout, 3733735, features, feature_names)
df = sp.build_catalog_features ('stellar_analysis_features')


chicken = sp.load_rooster_instance (filename='rooster_instances/rooster_tutorial')
(target_id, p_candidates,
 e_p_candidates, E_p_candidates,
 features, feature_names) = sp.create_rooster_feature_inputs (df, return_err=True)



rotation_score, prot, e_p, E_p = chicken.analyseSet (features, p_candidates, e_p_err=e_p_candidates,
                                                     E_p_err=E_p_candidates, feature_names=feature_names)


import plato_msap4_demonstrator_datasets.plato_sim_dataset as plato_sim_dataset

if not os.path.exists ('plato_sim_features') :
    os.mkdir ('plato_sim_features')

list_id = sp.get_list_targets (plato_sim_dataset)



def analysis_wrapper (star_id) :
    """
    Analysis wrapper to speed computation
    by parallelising process and control
    memory usage.
    """
    star_id = str (star_id).zfill (3)
    fileout = 'plato_sim_features/{}.csv'.format(star_id)
    fileplot = 'plato_sim_features/{}.png'.format(star_id)
    filename = sp.get_target_filename (plato_sim_dataset, star_id, filetype='csv')
    if not os.path.exists (fileout) :
        t, s, dt = sp.load_resource (filename)
        s = sp.preprocess (t, s, cut=60)
        (p_ps, p_acf,
         ps, acf,
         cs, features,
         feature_names, fig) = sp.analysis_pipeline (t, s, pmin=0.1, pmax=60,
                                                     wavelet_analysis=False, plot=True,
                                                     filename=fileplot, figsize=(10,16),
                                                     lw=1, dpi=150, pfa_threshold=1e-6,
                                                     ls_err_smooth=True)
        df = sp.save_features (fileout, star_id, features, feature_names)
        plt.close ("all")




process_pool = pathos.pools._ProcessPool (processes=4,
                                          maxtasksperchild=10)
with process_pool as p :
    list (tqdm (p.imap (analysis_wrapper,
                        list_id,
                        ),
                total=len (list_id))
          )
    p.close ()


df = sp.build_catalog_features ('plato_sim_features')

df


(target_id, p_candidates,
 e_p_candidates, E_p_candidates,
 features, feature_names) = sp.create_rooster_feature_inputs (df, return_err=True)
rotation_score, prot, e_p, E_p = chicken.analyseSet (features, p_candidates, e_p_err=e_p_candidates,
                                                     E_p_err=E_p_candidates, feature_names=feature_names)





prot_ref = sp.get_prot_ref (target_id, catalog='plato-sim')
cond_0 = (rotation_score>0.5)
cond_1 = (np.abs (prot - prot_ref) < 0.1 * prot_ref)
cond_2 = (np.abs (prot - prot_ref) < 0.1 * prot_ref) & (rotation_score>0.5)
score_0 = target_id[cond_0].size / target_id.size
score_1 = target_id[cond_1].size / target_id.size
score_2 = target_id[cond_2].size / target_id.size
score_0, score_1, score_2


fig, (ax1, ax2) = plt.subplots (1, 2, figsize=(10, 4))

bins = np.linspace (0, 1, 20, endpoint=False)
ax1.hist (rotation_score, bins=bins, color='darkorange')
ax1.axvline (0.5, ls='--', color='blue', lw=2)
bins = np.linspace (0, 80, 20, endpoint=False)
ax2.hist (prot, bins=bins, color='darkorange')
ax2.hist (prot_ref, bins=bins, facecolor='none',
         edgecolor='black', label='Ref')

ax1.set_ylabel (r'Number of stars')
ax1.set_xlabel (r'Rotation score')
ax2.set_xlabel (r'$P_\mathrm{rot}$ (day)')

ax1.set_xlim (0, 1)
ax2.set_xlim (0, 80)

plt.show()