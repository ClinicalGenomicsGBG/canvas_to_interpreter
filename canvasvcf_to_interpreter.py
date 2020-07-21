#!/usr/bin/env python

import argparse
import csv
import vcf
from Bio import SeqIO
import logging
import datetime


# Reads the vcf file
def read_vcf_new(infile):
    alts = []
    
    vcf_reader = vcf.Reader(open(infile, 'r'))
    return vcf_reader

# This function actually modifies the vcf to be compliant with Alissa Interpret. 
# It has to do some awkward stuff but eh... it works
def mod_canvas_vcf(vcf):
    refposlist = []
    fasta = "/medstore/External_References/hs37d5/hs37d5.fa"
    chrom_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))  
    vcf_list = list(vcf)
    for rec in vcf_list:        
        ref_pos = (str(rec.CHROM), str(rec.POS))
        refposlist.append(ref_pos)

    results = []
    for tup in refposlist: 
        presults = chrom_dict[tup[0]][int(tup[1]) - 1 : int(tup[1])]
        results.append({f"{tup[0]}:{tup[1]}": str(presults.seq)})
    
    refbases = {}
    for d in results:
        refbases.update(d)
   
    # This for loop will modify the actual variant lines in the vcf file
    for record in vcf_list:
        chrom = record.CHROM
        pos = record.POS
        cn = record.samples[0].data[3]
        mcc = record.samples[0].data[4]
        record.REF = refbases[f"{chrom}:{pos}"]
        
        try:
            chrom = int(chrom)
        except Exception:            
            pass
        if cn == 0 and type(chrom) == int:
            record.ALT[0] = "<DEL:HOM>"            
        elif cn == 0 and "Y" in chrom:
            record.ALT[0] = "<DEL:HEMI>"
        elif cn == 0 and "X" in chrom:
            # borde stå hemi om kön:man på x också
            record.ALT[0] = "<DEL:HOM>"
        elif cn == 1:
            record.ALT[0] = "<DEL>"
        elif cn == 2 and mcc == 2:
            record.ALT[0] = "<LOH>"
        elif cn == 2 and "Y" in chrom:
            record.ALT[0] = "<DUP>"
            ## Det borde vara så att det aldrig callas en cn == 2 variant i en kvinna och pga det behöver jag inte kolla om det är en man.
        elif cn == 2 and "X" in chrom:
            record.ALT[0] = "<DUP>"
        elif cn >= 3:
            record.ALT[0] = "<DUP>"
        else:
            # It seems like Alissa doesnt like <NORMAL> even though thats the suggestion from their docs ..
            record.ALT[0] = "<NORMAL>"
        logging.info(f"Variant at chromosome {chrom} pos {record.POS} has reference base '{record.REF}' and is set to type {record.ALT}")
    return vcf_list


# Also not used. Was used for testing
def print_vcf(vcf, rectype):
    print(rectype)
    try:
        for record in vcf:
            exec(f"print(record.{rectype})")
    except:
        exec(f"({vcf}.{rectype})")


# Writes the modified vcf
def write_vcf(fixed_vcf, output, template, length = False, qual = False, normals = False):
    fixed_vcf = iter(fixed_vcf)
    
    # This will modify the header values for ALT
    template.alts["DUP"] = ["<DUP>", template.alts["DUP"][1]]
    CN0copy = "<DEL>"
    CN0desc = "Deletion has happened (copy)"
    template.alts["CN0"] = [CN0copy, CN0desc]
    CN2copy = "<NORMAL>"
    CN2desc = "Copy number 2. Reference"
    template.alts["CN2"] = [CN2copy, CN2desc]
    CN3copy = "<DEL:HEMI>"
    CN3desc = "Hemizygous deletion"
    template.alts["CN3"] = [CN3copy, CN3desc]
    CN4copy = "<LOH>"
    CN4desc = "Loss of Heterozygosity"
    template.alts["CN4"] = [CN4copy, CN4desc]
    CN5copy = "<DEL:HOM>"
    CN5desc = "Homozygous deletion"
    template.alts["CN5"] = [CN5copy, CN5desc]
    
    writer = vcf.Writer(open(output, "w"), template)
    for rec in fixed_vcf:
        if "<NORMAL>" in rec.ALT and not normals:
            continue
        elif "FailedFT" in rec.FILTER:
            if "q7" in rec.samples[0].data[7] and not qual:
                logging.info(f"Quality filter failed: {rec.samples[0].data[7]}. Omitting this variant")
                continue
            elif "L10kb" in rec.samples[0].data[7] and not length:
                logging.info(f"Length filter failed: {rec.samples[0].data[7]}. Omitting this variant")
                continue
        writer.write_record(rec)


# Checks the header of the vcf to determine its origin
def determine_input(indata):
    with open(indata, "r") as vcfin:
        top = [next(vcfin) for x in range(100)]
        for line in top:
            if line.startswith("##source="):
                if "Canvas" in line:
                    return "Canvas"
                elif "Manta" in line:
                    return "Manta"
                else:
                    try:
                        return str(line)
                    except:
                        return None



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("vcf_input", nargs="?", type=str, help=':Full path to your vcf input file')
    parser.add_argument("vcf_output", nargs="?", type=str, help=':Specify full output path and filename of the vcf file you want to create')
    parser.add_argument("-r", "--record", nargs="?", action='store', type=str, help=':Specify which record ya wanna see. Only used for testing')
    parser.add_argument("-n", '--save_normals', action='store_true', help="Set this flag if you do not want to remove the normal (reference) variants from the outptu")
    parser.add_argument('-l', '--length_filtered', action='store_true', help="Set this flag if you do not want to remove the 'L10kb' filter-fail variants")
    parser.add_argument('-q', '--qual_filtered', action='store_true', help="Set this flag if you do not want to remove the 'qual below 7' filter-fail variants")
    parser.add_argument('-s', '--sex', action='store_true', help="Set sex of patient")
        

    args = parser.parse_args()
    
    logging.getLogger().setLevel(logging.INFO)
    inp = args.vcf_input
    outp = args.vcf_output
    rec = args.record
    length_filter = args.length_filtered
    save_normal = args.save_normals
    qual_filter = args.qual_filtered

    source_program = determine_input(inp)
    print(source_program)

    try:
        if "M" in args.sex.upper():
            gender = "Male"
        elif "F" in args.sex.upper() or "K" in args.sex.upper():
            gender = "Female"
        else:
            gender = "Unknown"
    except:
        gender = "Unknown"

    if source_program == "Canvas":
        fixed_pyvcf = mod_canvas_vcf(read_vcf_new(inp))
        write_vcf(fixed_pyvcf, outp, read_vcf_new(inp), length_filter, qual_filter, save_normal)
    elif source_program == "Manta":
        print("Manta not yet implemented ... sorry")
    else:
        try:
            print(f"The source of the vcf file could not be determined from:\n {line}")
        except:
            print("Source line not found. Source could not be determined.")

if __name__ == "__main__":
    main()
