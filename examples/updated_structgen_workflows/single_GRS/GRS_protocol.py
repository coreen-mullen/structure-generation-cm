import sys
from opt_tools import *
sys.path
import jax.numpy as np
import numpy as vnp
from jax import grad, jit, vmap
from jax import random
from functools import partial
import scipy as sp
# get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
import matplotlib.pyplot as plt
import lammps
import lammps.mliap
from lammps.mliap.loader import *
import pandas as pd
import copy
import pickle
import os
import shutil
import random
#plt.rcParams['figure.figsize']=[12,8]

#count descriptors from yace file
n_descs=get_desc_count('coupling_coefficients.yace')
elems=get_desc_count('coupling_coefficients.yace',return_elems=True)
nelements = len(elems)
dump_relax = True
data_path="./StructureDump"
cross_weight=0.0 #closeness to target as a joint distribution min/max just flip this sign (pos being minimize).
self_weight=1. #closeness to target within an atomic configuration. min/max just flip this sign (pos being minimize). 

soft_strength=1.0
ml_strength=1.0

# total number of configurations
n_totconfig=10
# max and minimum sizes for cells
maxcellsize=72
mincellsize=40
#  mask for your descriptors if you want (including all 11 of them here)
mask = list(range(n_descs))
target=vnp.zeros((len(mask)))
#composition in dictionary form: e.g. ( W:0.5, Be:0.5 )
target_comps = {'Ta':1.0}
# extra objective function variables to target exact matches of descriptors 
q1_wt=5.0
q1_cross_wt=5.0
q2_wt=0.0
q2_cross_wt=0.0
# flag to perform minimization with variable cell lengths
box_rlx = False
# flag to run at temperature first
run_temp = True
# flag to use fire style minimization (good for phase transitions & large minimizations)
fire_min = False
# flag to use linestyle quadratic minimization
line_min = True
# flag to use randomized compositions for elements in the dictionary: target_comps = {'Cr':1.0 }
randomize_comps=False
if not fire_min and box_rlx:
    min_typ_global = 'box'
else:
    min_typ_global = 'min'

#for now, read in target distribution (average and variance)
descriptors_target_1 = np.load('target_descriptors.npy')
descriptors_target_2 = np.load('target_var_descriptors.npy')

class GRSModel:
    def __init__(self, n_elements, n_descriptors_tot, mask):
        self.mask=mask
        self.n_descriptors=n_descriptors_tot
        self.n_descriptors_keep=len(mask)*n_elements
        self.n_elements=n_elements
        self.n_params=1
        self.sum_of_products=vnp.zeros((self.n_descriptors_keep,self.n_descriptors_keep))
        self.sum=vnp.zeros((self.n_descriptors_keep,))
        self.sumsq=vnp.zeros((self.n_descriptors_keep,))
        self.q_count=0
        self.first_moment_grad=grad(self.first_moment)
        self.V_grad=grad(self.V)
        self.K_self=self_weight
        self.K_cross=cross_weight
        self.first_mom_weight_cross = q1_cross_wt
        self.first_mom_weight = q1_wt
        self.second_mom_weight_cross = q2_cross_wt
        self.second_mom_weight = q2_wt
        self.whitening=np.identity(self.n_descriptors_keep)
        self.mode=  "update"
        self.data=[]

    def set_mode_update(self):
        self.mode="update"

    def set_mode_run(self):
        self.mode="run"

    def update(self,d):
        if self.n_elements > 1:
            dft=d.flatten()
        else:
            dft = d
        self.q_count += dft.shape[0]
        self.sum+=vnp.sum(dft,axis=0)
        self.sumsq+=vnp.sum(dft*dft,axis=0)
        self.set_mode_run()

    #match for mean descriptor value for the set of structures
    @partial(jit, static_argnums=(0))
    def first_moment_cross(self, descriptors):
        #n_atoms=descriptors.shape[0]
        #n_descriptors=descriptors.shape[1]
        avgs = self.sum/self.q_count
        abs_diffs_1 = np.abs(avgs - descriptors_target_1)
        abs_diffs_1 = np.array(abs_diffs_1)
        abs_diffs_1 = np.nan_to_num(abs_diffs_1)
        is_zero = np.isclose(abs_diffs_1,np.zeros(abs_diffs_1.shape))
        is_zero = np.array(is_zero,dtype=int)
        bonus=-np.sum(is_zero*self.first_mom_weight_cross)
        tst_residual_1 = np.sum(abs_diffs_1) +bonus
        return tst_residual_1

    #match for variance of descriptor value for the set of structures
    @partial(jit, static_argnums=(0))
    def second_moment_cross(self, descriptors):
        #n_atoms=descriptors.shape[0]
        #n_descriptors=descriptors.shape[1]
        vrs = self.sumsq/self.q_count 
        abs_diffs_2 = np.abs(vrs - descriptors_target_2)
        abs_diffs_2 = np.array(abs_diffs_2)
        abs_diffs_2 = np.nan_to_num(abs_diffs_2)
        is_zero = np.isclose(abs_diffs_2,np.zeros(abs_diffs_2.shape))
        is_zero = np.array(is_zero,dtype=int)
        bonus=-np.sum(is_zero*self.second_mom_weight_cross)
        tst_residual_2 = np.sum(abs_diffs_2) +bonus
        #print (tst_residual_2)
        return tst_residual_2

    # match of mean descriptor values within the current structure only
    @partial(jit, static_argnums=(0))
    def first_moment(self, descriptors):
        #n_atoms=descriptors.shape[0]
        #n_descriptors=descriptors.shape[1]
        avgs = np.average(descriptors,axis=0)
        abs_diffs_1 = np.abs(avgs - descriptors_target_1)
        abs_diffs_1 = np.array(abs_diffs_1)
        abs_diffs_1 = np.nan_to_num(abs_diffs_1)
        is_zero = np.isclose(abs_diffs_1,np.zeros(abs_diffs_1.shape))
        is_zero = np.array(is_zero,dtype=int)
        bonus=-np.sum(is_zero*self.first_mom_weight)
        tst_residual_1 = np.sum(abs_diffs_1) +bonus
        #print (tst_residual_1)
        return tst_residual_1

    # match of descriptor variance values within the current structure only
    @partial(jit, static_argnums=(0))
    def second_moment(self, descriptors):
        #n_atoms=descriptors.shape[0]
        #n_descriptors=descriptors.shape[1]
        vrs = np.var(descriptors,axis=0)
        abs_diffs_2 = np.abs(vrs - descriptors_target_2)
        abs_diffs_2 = np.array(abs_diffs_2)
        abs_diffs_2 = np.nan_to_num(abs_diffs_2)
        is_zero = np.isclose(abs_diffs_2,np.zeros(abs_diffs_2.shape))
        is_zero = np.array(is_zero,dtype=int)
        bonus=-np.sum(is_zero*self.second_mom_weight)
        tst_residual_2 = np.sum(abs_diffs_2) +bonus
        #print (tst_residual_2)
        return tst_residual_2

    #note that the current default weights for this "potential" turn off the variance contribution
    @partial(jit, static_argnums=(0))
    def V(self,descriptors,weights=[1.0,0.0]):
        if self.n_elements > 1:
            descriptors_flt = descriptors.flatten()
        else:
            descriptors_flt = descriptors
        vi = ((weights[0]*self.first_moment(descriptors_flt)) + (weights[1]*self.second_moment(descriptors_flt)))
        vj = ((weights[0]*self.first_moment_cross(descriptors_flt)) + (weights[1]*self.second_moment_cross(descriptors_flt)))
        return self.K_self*vi + self.K_cross*vj

    def __call__(self, elems, descriptors, beta, energy):
        self.last_descriptors=descriptors.copy()
        if self.mode=="run":
            b=descriptors[:,self.mask]
            ener=self.V(b)
            energy[:]=0
            energy[0]=ener
            b=self.V_grad(b)
            if not np.all(np.isfinite(b)):
                print("GRAD ERROR!")
                #print(b)

            beta[:,:]=0
            beta[:,self.mask]=b

        if self.mode=="update":
            b=descriptors[:,self.mask]
            self.update(b)

class GRSSampler:
    def __init__(self, model, before_loading_init):
        self.model=model
        self.lmp = lammps.lammps(cmdargs=['-screen','none'])
        lammps.mliap.activate_mliappy(self.lmp)
        self.lmp.commands_string(before_loading_init)
        lammps.mliap.load_model(em)

    def update_model(self):
        self.model.set_mode_update()
        self.lmp.commands_string("variable etot equal etotal")
        self.lmp.commands_string("variable ptot equal press")
        self.lmp.commands_string("variable pairp equal epair")
        self.lmp.commands_string("variable numat equal atoms")
        self.lmp.commands_string("run 0")
        self.lmp.commands_string("print \"${etot} ${pairp} ${ptot} ${numat} \" append Summary.dat screen no")
        self.model.set_mode_run()

    def run(self,cmd=None):
        if cmd==None:
            self.lmp.commands_string("run 0")
        else:
            self.lmp.commands_string(cmd)

#begin the actual optimization loop
try:
    shutil.rmtree(data_path)
except:
    pass
try:
    os.mkdir(data_path)
except:
    pass
lmpi = lammps.lammps(cmdargs=['-screen','none'])
activate_mliappy(lmpi)
em=GRSModel(nelements,n_descs,mask=mask)
# use internal structure generator to build candidates
atoms_basic= bulk('Ta',cubic=True)*(3,3,3)#read('ats_00057.cif')#internal_generate_cell(0,desired_size=vnp.random.choice(range(mincellsize,maxcellsize)),template=None,desired_comps=target_comps,use_template=None,min_typ=min_typ_global,soft_strength=soft_strength)
print('pos',atoms_basic.positions)
atoms_basic.rattle(stdev=0.5,seed=vnp.random.randint(40000))
print('pos2',atoms_basic.positions)
#g= internal_basic(atoms_basic,index=0,min_typ=None,soft_strength=1.0)
#g = internal_generate_cell(0,desired_size=vnp.random.choice(range(mincellsize,maxcellsize)),template=atoms_basic,desired_comps=target_comps,use_template=atoms_basic,min_typ=min_typ_global,soft_strength=soft_strength)
g = internal_generate_cell(0,desired_size=vnp.random.choice(range(mincellsize,maxcellsize)),template=None,desired_comps=target_comps,use_template=None,min_typ=min_typ_global,soft_strength=soft_strength)
print ('model',em)
sampler=GRSSampler(em,g)
print ('after sampler')
sampler.update_model()

i=1
while i <= n_totconfig:
    print(i,"/",n_totconfig,"Using indicies :",mask)
    atoms_basic.rattle(stdev=0.5,seed=vnp.random.randint(40000))
    #g = read('ats_00057.cif')
    #g= internal_basic(atoms_basic,index=0,min_typ=None,soft_strength=1.0)
    #g = internal_generate_cell(i,desired_size=vnp.random.choice(range(mincellsize,maxcellsize)),template=atoms_basic,desired_comps=target_comps,use_template=atoms_basic,min_typ=min_typ_global,soft_strength=soft_strength)
    g = internal_generate_cell(i,desired_size=vnp.random.choice(range(mincellsize,maxcellsize)),template=None,desired_comps=target_comps,use_template=None,min_typ=min_typ_global,soft_strength=soft_strength)
    sampler=GRSSampler(em,g)
    em.K_cross=cross_weight
    em.K_self=self_weight
    if run_temp:
        if dump_relax:
            sampler.run("dump    r%d all xyz 1 dyn_%d.xyz" %(i,i))
        sampler.run('velocity all create 800.0 4928459 loop geom')
        sampler.run("fix  a1  all nve")
        sampler.run("run 1000")
        sampler.run("unfix  a1")
    if fire_min:
        sampler.run("min_style  fire")
        sampler.run("""min_modify integrator eulerexplicit tmax 10.0 tmin 0.0 delaystep 5 dtgrow 1.1 dtshrink 0.5 alpha0 0.1 alphashrink 0.99 vdfmax 100000 halfstepback no initialdelay no""")
    if line_min:
        sampler.run("min_style  cg")
        sampler.run("min_modify  dmax 0.5 line quadratic")


    #sampler.run("minimize 1e-12 1e-12 10000 100000")
    #sampler.run("minimize 1e-6 1e-6 10000 100000")
    if dump_relax:
        sampler.run("dump    a%d all xyz 1 min_%d.xyz" %(i,i))
    sampler.run("minimize 1e-12 1e-12 10000 100000")
    sampler.run("write_data %s/sample.%i.dat " % (data_path,i) )
    sampler.update_model()

    i+=1


