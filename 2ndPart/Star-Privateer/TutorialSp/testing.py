import star_privateer as sp
import plato_msap4_demonstrator_datasets.kepler_dataset as kepler_dataset
import os, pathos
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

if not os.path.exists ('rooster_training_features') :
    os.mkdir ('rooster_training_features')
if not os.path.exists ('rooster_instances') :
    os.mkdir ('rooster_instances')

list_kic = sp.get_list_targets (kepler_dataset)


def analysis_wrapper (kic) :
    """
    Analysis wrapper to speed computation
    by parallelising process and control
    memory usage.
    """
    str_kic = str (kic).zfill (9)
    filename = sp.get_target_filename (kepler_dataset, str_kic)
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

    print(df)

    df.to_csv ("training_features.csv")
   


df_train = df.sample (n=df.index.size//2, random_state=49458493)
df_test = df.loc[np.setdiff1d (df.index, df_train.index)]


(training_id, training_p_candidates,
 training_features, feature_names) = sp.create_rooster_feature_inputs (df_train)
(test_id, test_p_candidates,
 test_features, test_feature_names) = sp.create_rooster_feature_inputs (df_test)

feature_names

seed = 104359357
chicken = sp.ROOSTER (n_estimators=100, random_state=np.random.RandomState (seed=seed))
chicken.RotClass, chicken.PeriodSel

chicken.train (training_id, training_p_candidates,
               training_features, feature_names=feature_names,
               catalog='santos-19-21', verbose=True)

results = chicken.test (test_id, test_p_candidates, test_features,
                        feature_names=test_feature_names,
                        catalog='santos-19-21', verbose=True)

chicken.getScore ()
chicken.getNumberEltTrain ()
chicken.getNumberEltTest ()


prot_rooster = results[3]
prot_ref = sp.get_prot_ref (results[2], catalog='santos-19-21')

fig, ax = plt.subplots (1, 1)

bins = np.linspace (0, 80, 20, endpoint=False)

ax.hist (prot_rooster, bins=bins, color='darkorange', label='ROOSTER')
ax.hist (prot_ref, bins=bins, facecolor='none',
        edgecolor='black', label='Ref')

ax.set_xlabel (r'$P_\mathrm{rot}$ (day)')
ax.set_ylabel (r'Number of stars')

ax.legend ()

fig, ax = plt.subplots (1, 1, figsize=(5, 5))

ax.scatter (prot_ref, prot_rooster,
            color='darkorange', s=3, marker="o")

ax.set_xlabel (r'$P_\mathrm{rot, true}$ (day)')
ax.set_ylabel (r"$P_\mathrm{rot, ROOSTER}$ (day)")

ax.plot ([0, 60], [0, 60], ls="--", color="grey")

fig.tight_layout ()


fig, (ax, ax0) = plt.subplots (1, 2, figsize=(6, 4),
                               width_ratios=[0.8, 0.2],
                               sharey=True)

ax.scatter (prot_ref, (prot_rooster - prot_ref) / prot_ref * 100,
            color='darkorange', s=3, marker="o")

ax0.hist ((prot_rooster - prot_ref) / prot_ref * 100,
          bins=np.linspace (-20, 20, 31), orientation="horizontal",
          color="darkorange")

ax.set_xlabel (r'$P_\mathrm{rot, true}$ (day)')
ax.set_ylabel (r"$\delta P_\mathrm{rot}$ (%)")

ax.axhline (0, ls="--", color="grey")

ax.set_ylim (-10, 10)

ax0.set_xlim (0, 200)
ax0.set_xlabel (r"$N_\mathrm{stars}$")

plt.show()


chicken.save ('rooster_instances/rooster_tutorial')

