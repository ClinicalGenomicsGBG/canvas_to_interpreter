# canvas_to_interpreter
Converts Canvas VCF-format to format accepted by variant interpreter

## Requirements
Runs successfully on Python 3.7.3  
PyVCF==0.6.8  
biopython==1.73  

## Run on cluster
module load anaconda2  
source activate canvas  

python canvasvcf_to_interpreter_reftest.py -l -q <path to _CNV_germline.vcf> <path to output _CNV_germline_alissaformat.vcf> <path to reference genome>

