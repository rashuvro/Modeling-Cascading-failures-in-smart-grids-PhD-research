1. First run, "FindingStateSpaceMultiV2" to generate the line failure data.
You can specify Fi, r, e, \theta and capacity values at the bengining of the program to generate different data.
Line Failure data will be saved with the file name given at the end of the program.

2. Then run "checkdata4" to generate Table 1 and Figure 1 of the paper.
Note that, at first in "checkdata4" you need to provide the file input (variable "filename") you got at step 1 and
aslo specify the capacities (variable "TrueCaps") you have used in step 1.
Failure probability data will be saved in a file named "failProb", which we will use for curve fitting.

3. To generate Table 2 run "numFailedLinesEachTypeRatio".
first you need to provide the file input (variable "filename") you got at step 1 and
aslo specify the capacities (variable "TrueCaps") you have used in step 1.

4. Then run "checkdata_numStepsToSS_v2" to generate figure 2(a) and Fig. 3 of the paper.
first, you need to provide the file input (variable "filename") you got at step 1.

5. To generate figure 2(b) run "numFailedLinesEachType". data will be saved in a excel file named "numFailedType".
Then run "power_grid_paper_data" in python to generate Fig. 2(b). 
first, you need to provide the file input (variable "filename") you got at step 1.

6. To generate figure 4 run "curveFit6_03092018". Provide filename input of failure probability you got 
from step 2.

7. Run "checkEquation_02112018" to generate Fig. 5 and Fig. 6 (a).

8. Run "newDiffCapFailedLinesEachTimeStep02252018" to generate Fig. 6 (b).

9. Run "checkEquation_02112018" to generate Fig. 7(b). However, in this case you need to uncomment for loop
loop that is used for kappa and also uncomment Fig. 3 in the "checkEquation_02112018" program.

