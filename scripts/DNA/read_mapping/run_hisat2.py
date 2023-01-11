
import os

#parameters retrieved from ../scripts/get_mapping_params.py
indv_params = (
    #(fastq_prefix, PU_value, SM_value, rg_id, sam_output_filename, [PAIRED/UNPAIRED]),
    ('Amp_177.0', 'C6WMHACXX:3:none', 'Ethan', 'Amp_0', 'Amp_177.0.sam', 'PAIRED'),
    ('Amp_177.1', 'C6WMHACXX:4:none', 'Ethan', 'Amp_1', 'Amp_177.1.sam', 'PAIRED'),
    ('Amp_259.0', 'C6WMHACXX:3:none', 'Eric', 'Amp_0', 'Amp_259.0.sam', 'PAIRED'),
    ('Amp_259.1', 'C6WMHACXX:4:none', 'Eric', 'Amp_1', 'Amp_259.1.sam', 'PAIRED'),
    ('Amp_261.0', 'C6WMHACXX:3:none', 'Elizabeth', 'Amp_0', 'Amp_261.0.sam', 'PAIRED'),
    ('Amp_261.1', 'C6WMHACXX:4:none', 'Elizabeth', 'Amp_1', 'Amp_261.1.sam', 'PAIRED'),
    ('Amp_272.0', 'C6WMHACXX:3:none', 'Eli', 'Amp_0', 'Amp_272.0.sam', 'PAIRED'),
    ('Amp_272.1', 'C6WMHACXX:4:none', 'Eli', 'Amp_1', 'Amp_272.1.sam', 'PAIRED'),
    ('LP6005672-DNA_A01.0', 'C32JPACXX_0:5:none', 'Alica', '0', 'LP6005672-DNA_A01.0.sam', 'PAIRED'),
    ('LP6005672-DNA_A01.1', 'C32JPACXX_0:6:none', 'Alica', '1', 'LP6005672-DNA_A01.1.sam', 'PAIRED'),
    ('LP6005672-DNA_A01.2', 'C32JPACXX_0:7:none', 'Alica', '2', 'LP6005672-DNA_A01.2.sam', 'PAIRED'),
    ('LP6005672-DNA_B01.0', 'C31HEACXX_0:5:none', 'Anne', '0', 'LP6005672-DNA_B01.0.sam', 'PAIRED'),
    ('LP6005672-DNA_B01.1', 'C31HEACXX_0:6:none', 'Anne', '1', 'LP6005672-DNA_B01.1.sam', 'PAIRED'),
    ('LP6005672-DNA_B01.2', 'C31HEACXX_0:7:none', 'Anne', '2', 'LP6005672-DNA_B01.2.sam', 'PAIRED'),
    ('LP6005672-DNA_C01.0', 'C31HEACXX_0:8:none', 'Adam', '0', 'LP6005672-DNA_C01.0.sam', 'PAIRED'),
    ('LP6005672-DNA_C01.1', 'C31YYACXX_1:7:none', 'Adam', '1', 'LP6005672-DNA_C01.1.sam', 'PAIRED'),
    ('LP6005672-DNA_C01.2', 'C31YYACXX_1:8:none', 'Adam', '2', 'LP6005672-DNA_C01.2.sam', 'PAIRED'),
    ('LP6005801-DNA_A01.0', 'C390HACXX_0:8:none', 'Brad', '0', 'LP6005801-DNA_A01.0.sam', 'PAIRED'),
    ('LP6005801-DNA_A01.1', 'C38YAACXX_1:1:none', 'Brad', '1', 'LP6005801-DNA_A01.1.sam', 'PAIRED'),
    ('LP6005801-DNA_A01.2', 'C38YAACXX_1:2:none', 'Brad', '2', 'LP6005801-DNA_A01.2.sam', 'PAIRED'),
    ('LP6005801-DNA_A02.0', 'C373KACXX_0:8:none', 'DeAnne', '0', 'LP6005801-DNA_A02.0.sam', 'PAIRED'),
    ('LP6005801-DNA_A02.1', 'C37GBACXX_1:1:none', 'DeAnne', '1', 'LP6005801-DNA_A02.1.sam', 'PAIRED'),
    ('LP6005801-DNA_A02.2', 'C37GBACXX_1:2:none', 'DeAnne', '2', 'LP6005801-DNA_A02.2.sam', 'PAIRED'),
    ('LP6005801-DNA_B01.0', 'C38YAACXX_0:3:none', 'Betty', '0', 'LP6005801-DNA_B01.0.sam', 'PAIRED'),
    ('LP6005801-DNA_B01.1', 'C38YAACXX_0:4:none', 'Betty', '1', 'LP6005801-DNA_B01.1.sam', 'PAIRED'),
    ('LP6005801-DNA_B01.2', 'C38YAACXX_0:5:none', 'Betty', '2', 'LP6005801-DNA_B01.2.sam', 'PAIRED'),
    ('LP6005801-DNA_B02.0', 'C37GBACXX_0:3:none', 'Elizabeth', '0', 'LP6005801-DNA_B02.0.sam', 'PAIRED'),
    ('LP6005801-DNA_B02.1', 'C37GBACXX_0:4:none', 'Elizabeth', '1', 'LP6005801-DNA_B02.1.sam', 'PAIRED'),
    ('LP6005801-DNA_B02.2', 'C37GBACXX_0:5:none', 'Elizabeth', '2', 'LP6005801-DNA_B02.2.sam', 'PAIRED'),
    ('LP6005801-DNA_C01.0', 'C38YAACXX_0:6:none', 'Chris', '0', 'LP6005801-DNA_C01.0.sam', 'PAIRED'),
    ('LP6005801-DNA_C01.1', 'C38YAACXX_0:7:none', 'Chris', '1', 'LP6005801-DNA_C01.1.sam', 'PAIRED'),
    ('LP6005801-DNA_C01.2', 'C38YAACXX_0:8:none', 'Chris', '2', 'LP6005801-DNA_C01.2.sam', 'PAIRED'),
    ('LP6005801-DNA_C02.0', 'C37GBACXX_0:6:none', 'Eli', '0', 'LP6005801-DNA_C02.0.sam', 'PAIRED'),
    ('LP6005801-DNA_C02.1', 'C37GBACXX_0:7:none', 'Eli', '1', 'LP6005801-DNA_C02.1.sam', 'PAIRED'),
    ('LP6005801-DNA_C02.2', 'C37GBACXX_0:8:none', 'Eli', '2', 'LP6005801-DNA_C02.2.sam', 'PAIRED'),
    ('LP6005801-DNA_D01.0', 'C38MJACXX_0:1:none', 'Catherine', '0', 'LP6005801-DNA_D01.0.sam', 'PAIRED'),
    ('LP6005801-DNA_D01.1', 'C38MJACXX_0:2:none', 'Catherine', '1', 'LP6005801-DNA_D01.1.sam', 'PAIRED'),
    ('LP6005801-DNA_D01.2', 'C38MJACXX_0:3:none', 'Catherine', '2', 'LP6005801-DNA_D01.2.sam', 'PAIRED'),
    ('LP6005801-DNA_D02.0', 'C39ALACXX_0:1:none', 'Bob', '0', 'LP6005801-DNA_D02.0.sam', 'PAIRED'),
    ('LP6005801-DNA_D02.1', 'C39ALACXX_0:2:none', 'Bob', '1', 'LP6005801-DNA_D02.1.sam', 'PAIRED'),
    ('LP6005801-DNA_D02.2', 'C39ALACXX_0:3:none', 'Bob', '2', 'LP6005801-DNA_D02.2.sam', 'PAIRED'),
    ('LP6005801-DNA_E01.0', 'C38MJACXX_0:4:none', 'Carl', '0', 'LP6005801-DNA_E01.0.sam', 'PAIRED'),
    ('LP6005801-DNA_E01.1', 'C38MJACXX_0:5:none', 'Carl', '1', 'LP6005801-DNA_E01.1.sam', 'PAIRED'),
    ('LP6005801-DNA_E01.2', 'C38MJACXX_0:6:none', 'Carl', '2', 'LP6005801-DNA_E01.2.sam', 'PAIRED'),
    ('LP6005801-DNA_F01.0', 'C373KACXX_0:1:none', 'Dave', '0', 'LP6005801-DNA_F01.0.sam', 'PAIRED'),
    ('LP6005801-DNA_F01.1', 'C38MJACXX_1:7:none', 'Dave', '1', 'LP6005801-DNA_F01.1.sam', 'PAIRED'),
    ('LP6005801-DNA_F01.2', 'C38MJACXX_1:8:none', 'Dave', '2', 'LP6005801-DNA_F01.2.sam', 'PAIRED'),
    ('LP6005801-DNA_G01.0', 'C36MNACXX_0:6:none', 'Douglas', '0', 'LP6005801-DNA_G01.0.sam', 'PAIRED'),
    ('LP6005801-DNA_G01.1', 'C36MNACXX_0:7:none', 'Douglas', '1', 'LP6005801-DNA_G01.1.sam', 'PAIRED'),
    ('LP6005801-DNA_G01.2', 'C36MNACXX_0:8:none', 'Douglas', '2', 'LP6005801-DNA_G01.2.sam', 'PAIRED'),
    ('LP6005801-DNA_H01.0', 'C373KACXX_0:5:none', 'Ethan', '0', 'LP6005801-DNA_H01.0.sam', 'PAIRED'),
    ('LP6005801-DNA_H01.1', 'C373KACXX_0:6:none', 'Ethan', '1', 'LP6005801-DNA_H01.1.sam', 'PAIRED'),
    ('LP6005801-DNA_H01.2', 'C373KACXX_0:7:none', 'Ethan', '2', 'LP6005801-DNA_H01.2.sam', 'PAIRED'),
    ('NA12878.0', 'C0L54ACXX_0:4:none', 'Wendy', '0', 'NA12878.0.sam', 'PAIRED'),
    ('NA12878.1', 'C0L54ACXX_0:5:none', 'Wendy', '1', 'NA12878.1.sam', 'PAIRED'),
    ('NA12878.2', 'C0L54ACXX_0:6:none', 'Wendy', '2', 'NA12878.2.sam', 'PAIRED'),
    ('NA12891.0', 'C0L0AACXX_0:1:none', 'Wanda', '0', 'NA12891.0.sam', 'PAIRED'),
    ('NA12891.1', 'C0L54ACXX_1:7:none', 'Wanda', '1', 'NA12891.1.sam', 'PAIRED'),
    ('NA12891.2', 'C0L54ACXX_1:8:none', 'Wanda', '2', 'NA12891.2.sam', 'PAIRED'),
    ('NA12892.0', 'C0L0AACXX_0:2:none', 'Wilber', '0', 'NA12892.0.sam', 'PAIRED'),
    ('NA12892.1', 'C0L0AACXX_0:3:none', 'Wilber', '1', 'NA12892.1.sam', 'PAIRED'),
    ('NA12892.2', 'C0L0AACXX_0:4:none', 'Wilber', '2', 'NA12892.2.sam', 'PAIRED'),

    #beads
    ('Om_177.0', 'C6WMHACXX:3:none', 'Ethan', 'Om_177_0', 'Om_177.0.sam', 'PAIRED'),
    ('Om_177.1', 'C6WMHACXX:4:none', 'Ethan', 'Om_177_1', 'Om_177.1.sam', 'PAIRED'),
    ('Om_259.0', 'C6WMHACXX:3:none', 'Eric', 'Om_259_0', 'Om_259.0.sam', 'PAIRED'),
    ('Om_259.1', 'C6WMHACXX:4:none', 'Eric', 'Om_259_1', 'Om_259.1.sam', 'PAIRED'),
    ('Om_261.0', 'C6WMHACXX:3:none', 'Elizabeth', 'Om_261_0', 'Om_261.0.sam', 'PAIRED'),
    ('Om_261.1', 'C6WMHACXX:4:none', 'Elizabeth', 'Om_261_1', 'Om_261.1.sam', 'PAIRED'),
    ('Om_272.0', 'C6WMHACXX:3:none', 'Eli', 'Om_272_0', 'Om_272.0.sam', 'PAIRED'),
    ('Om_272.1', 'C6WMHACXX:4:none', 'Eli', 'Om_272_1', 'Om_272.1.sam', 'PAIRED'),

    #new reads for bad dads
    ('bob.i5_ACAGTG_L008_R1_001', 'C86JUANXX:8:none', 'Bob', '4', 'bob.05-2016-reads.sam', 'UNPAIRED'),
    ('doug.i3_TTAGGC_L008_R1_001', 'C86JUANXX:8:none', 'Douglas', '4', 'doug.05-2016-reads.sam', 'UNPAIRED'),
)

#indv_params = (
    #(fastq_prefix, PU_value, SM_value, rg_id, sam_output_filename, [PAIRED/UNPAIRED]),
#    ('LP6005672-DNA_B01.0', 'C31HEACXX_0:5:none', 'LP6005672-DNA_B01', '0', 'LP6005672-DNA_B01.0.sam', 'PAIRED'),
#    ('LP6005801-DNA_D01.1', 'C38MJACXX_0:2:none', 'LP6005801-DNA_D01', '1', 'LP6005801-DNA_D01.1.sam', 'PAIRED'),
#    ('LP6005801-DNA_E01.2', 'C38MJACXX_0:6:none', 'LP6005801-DNA_E01', '2', 'LP6005801-DNA_E01.2.sam', 'PAIRED'),
#    ('LP6005801-DNA_G01.0', 'C36MNACXX_0:6:none', 'LP6005801-DNA_G01', '0', 'LP6005801-DNA_G01.0.sam', 'PAIRED'),
#)


if __name__=="__main__":
    for fastq_pre, pu, sm, rg_id, sam, read_type in indv_params:
        os.system('sbatch --export=FASTQPRE="' + fastq_pre + '",SAM="' + sam + '",PU="' + pu + '",SM="' + sm + '",RGID="' + rg_id + 
              '",TYPE="' + read_type + '" mappinghisat2.sh')


