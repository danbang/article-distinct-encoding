# article-distinct-encoding

This repository contains anonymised data and code supporting:

<a href="http://www.pnas.org/content/early/2018/05/18/1800795115">Bang & Fleming (2018) “Distinct encoding of decision confidence in human medial prefrontal cortex” PNAS </a>

FigureX.m files will generate the specified plots from the paper by loading behavioural data, model predictions and/or fMRI ROI data.

The paths in these scripts require altering the directory variable ‘repoBase’ to point to your local version of the repository.

Unthresholded group-level statistical maps supporting whole-brain analyses are available on NeuroVault: https://neurovault.org/collections/3792/.

**Folders**

*Custom*

Contains custom scripts for data extraction and plotting.

*Data*

Contains anonymised behavioural data for the pre-scan (Data/Behaviour/Prescen) and the scan (Data/Behaviour/Prescen) sessions, ROI contrast estimates (Data/Scan/Contrast), ROI single-trial activity estimates (Data/Scan/Activity) and ROI activity time courses (Data/Scan/Timecourse).

*HDDM*

Contains scripts for hierarchical DDM fitting using the HDDM toolbox (http://ski.clps.brown.edu/hddm_docs/) and generation of posterior predictives using the DMAT toolbox (https://ppw.kuleuven.be/okp/software/dmat/). 

*Model*

Contains scripts for fitting confidence model to pre-scan data and generating out-of-sample predictions about subjective confidence for scan data.

*Permutation*

Contains scripts for permutation testing for analysis of ROI activity time courses.

*SPM*

Contains scripts for setting up GLMs of fMRI data. fMRI data were pre-processed using standard pipelines available here: https://github.com/metacoglab/MetaLabCore. All fMRI analyses were done using SPM 12 (www.fil.ion.ucl.ac.uk/spm). 

*Upsample*

Contains example script for creating ROI activity time courses.

*Preprint*

A preprint of this article can be found on <a href="https://www.biorxiv.org/content/early/2018/05/03/251330">bioRxiv</a>.

This code is being released with a permissive open-source license. You should feel free to use or adapt the  code as long as you follow the terms of the license. If you make use of the behavioural, modelling or neuroimaging analyses, we would appreciate that you cite the paper.
