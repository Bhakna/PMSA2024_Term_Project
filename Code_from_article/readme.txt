At the very top of "main.m", use 

case_study = @ cases.CaseStudyOne;
case_study = @ cases.CaseStudyTwo;
case_study = @ cases.CaseStudyThree;

respectively to switch between the three case studies in the paper. The defined function is called at row 40 in "main.m" to fetch the case study parameters.

We have not tried all edge cases and cannot guarantee the simulation works for various extreme cases of particle size, physiology, etc. Feel free to contact us in case you run into trouble, we will do our best to help.

Notes:

* Currently particles are normally distributed, defined by mu/sigma fetched from the case study at row 40 in "main.m". The distribution is created in "+deposition\Deposition.m" at row 48. This is where you should insert your own custom distribution if you would like to work with something else than normal.

* The lung deposited dose "lung_dose_ug" and final time "T" are also fetched from the case study at row 40 in "main.m".

*Alternately Yeh/Schum physiology can be loaded by use of "anatomy = 'lung_yeh_schum.csv';" instead of "anatomy = 'lung_weibel.csv';" at the very top of "main.m".

* Outside the currently specified particle size range (0.00015*1e-5 to 20*1e-5 [dm], diameter) the deposition equations may be invalid. These bounds could change if another lung geometry is used.