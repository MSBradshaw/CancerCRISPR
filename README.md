# CancerCRISPR

## How to use it
`sliding_window.py` requires the following parameters:
* --vcf VCF BGZIPED and TABIXED vcf file
* --plot_output PLOT_OUTPUT path for where to save plot, if not used plot is not generated
* --top_x TOP_X         number of highest z score windows to report on (default 10)
* --output OUTPUT       path for output file to be saved
* --dp DP               DP treshold for counting a SNP in the vcf as a real mutation (default 10)
* --window_size WINDOW  window size (default 1000)
* --shift_size SHIFT    amount to slide the window by (default 100), should be <= than the window size

### Example:
Example data can be download from this publicly accessible [google drive](https://drive.google.com/file/d/15dYuplueDpu7giw3ysCdMh2FR_uihCeY/view) https://drive.google.com/file/d/15dYuplueDpu7giw3ysCdMh2FR_uihCeY/view


`python sliding_window.py --vcf NCIH2170_LUNG.g.vcf.gz --plot_output sliding_window.png --top_x 10 --dp 10 --window_size 1000 --shift_size 100 --output output.txt`
