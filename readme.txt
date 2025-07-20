MultiHollowPtycho

# MultiHollowPtycho

This source code is a modified version of the PtychoShelves package originally developed at Paul Scherrer Institut (PSI). It has been adapted for hollow multi-slice electron ptychographic reconstruction at the University of Warwick.

If this code, its datasets, or any derivatives are used in your research, please cite the following publication:

Yu Lei and Peng Wang, Hollow multi-slice electron ptychography for simultaneous 3D structural imaging and EELS in 4D-STEM, arXiv:2506.22352 (https://arxiv.org/abs/2506.22352)

You should also cite the original references for the base methods:

K. Wakonig et al., PtychoShelves, a versatile high-level framework for high-performance analysis of ptychographic data, J. Appl. Cryst. 53(2) (2020)
M. Odstrcil et al., Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics Express, 26(3):3108–3123 (2018)
Z. Chen et al., Electron ptychography achieves atomic-resolution limits set by lattice vibrations, arXiv:2101.00465

---2025/7/18 Yu Lei, Dr.Peng Wang


#######################

- The code is streamlined from the original PtychoShelves code for demonstration purpose only. 

- The complete package with many more algorithms and documentation can be found at https://www.psi.ch/en/sls/csaxs/software

- Copyright and license issues should follow the agreements in PSI's codes and/or refer to PtychoShelves website.  

#######################

Short notes for startup:

1. Large diffraction data needs to be loaded from PARADIM website: https://data.paradim.org/, and should be put in '../PtychoShelves_EM_Hollow/ptycho/Hollow_Data_exp_30nm/'.

2. For an initial test of the reconstruction, change the Matlab current Folder (work path) to ../PtychoShelves_EM/ptycho/, and run the main drive ptychographt_exp_30nm.m. 

3. For your own data, prepare a data file 'data_dp.mat' containing diffractions, 'data_position.hdf5' containing probe positions, and initial probe 'probe_initial.mat'. Please see the example script: ../ptycho/utils_EM/prepare_data_electron_exp.m. 

4. look into the script ptychographt_exp_30nm.m and modify parameters accordingly. More options can be refered to the full documentation for PtychoShelves. 

5. Only Matlab 2023b with the GPU engine has been tested. 

6. Further improvements may be acheived by using more probe modes, more layers, better sampling, or a better initial start of probe and probe positions. But keep in mind that using large dimensions is very compuational demanding in both memory and time.

#########################

Contact information:

Yu Lei (Yu.Lei.3@warwick.ac.uk) or Peng Wang (Peng.Wang.3@warwick.ac.uk)
