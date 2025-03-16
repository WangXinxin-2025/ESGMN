![image](https://github.com/user-attachments/assets/f6473e00-057c-4ede-b24f-b355b7da5dcb)# ESGMN
An enhanced structure-guided molecular networking (E-SGMN) method was developed, which is specifically tailored for the Orbitrap Astral mass spectrometer (MS).

ESGMN User Manual
--
* *1. Introduction and Installation<br>

* *2. Input files<br>

	* *2.1 Background database<br>

	* *2.2 Files for structure annotation<br>

* *3. Software running instructions<br>

	* *3.1 Parameter settings for ESGMN<br>

	* *3.2 Metabolite annotation<br>
		* *3.2.1 Using the default background database
		* *3.2.2 Building a new background network


 
# 1. Introduction and Installation

* Structure-Guided Molecular Network Strategy (ESGMN) is a free network-based and spectral library-independent tool for deep annotation of untargeted ultra-performance liquid chromatography-high resolution mass spectrometry (UPLC-HRMS) metabolomics data. The source program is in the Python language, and the development environment is Spyder 3.7. ESGMN is published as an executable file and can be downloaded freely from the SourceForge via https://sourceforge.net/projects/ESGMN/files/ (All five RAR files need to be downloaded).

# 2. Input files

## 2.1 Background database
* If you want to use the default background network, you can skip this step.<br>
* If you want to build a new background network instead of using the default one for the subsequent metabolite annotation, you need to provide a background database form. This form should be stored in csv format file and be named ‘Background database_POS.csv’ for positive ion mode and ‘Background database_NEG.csv’ for negative ion mode (Figure 1). The background database contains the ID, SMILES, exact mass, retention time (RT/min), formula, and Name. 
