How to run the program:
1)from codes folder change path according to your system in the 'c' codes.
2)set the step size and current start and ending value in "par_diag.h" and "par_HH.h"
3)compile and run "HHneuron6Sep2014.c" 
4)compile and run "SpikeStatistics_18Nov2014.c"

Note: for plotting of the data refer to gnucodes in Results folder which are included in the .odt files

New Additions:
The updated version of the code generates an extra text file "len_freq" this file contains frequency of lengths in the range of current values selected. The first line has freqency of length 0 , The second line has frequency of length 1 and so on. By changing the Num_lengths variable in "SpikeStatistics_18Nov2014.c" the max length recorded can be changed currently the max length is set to 60.
