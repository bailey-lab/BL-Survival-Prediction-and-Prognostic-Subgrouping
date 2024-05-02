import pandas as pd
import sys
import math

tsv_file = open(sys.argv[1], 'r')
df = pd.read_table(tsv_file)

sample_file = open(sys.argv[2], 'r')
sample_file_lines = sample_file.readlines()

raw_samples = set()
for line in sample_file_lines:
    raw_samples.add(line.rstrip('\n'))

known_samples_set = set()
for raw_sample in raw_samples:
    string = str(raw_sample) + ('_GT')
    known_samples_set.add(string)


variant_dict = {}

num_var = set()
all_variants_set = set()

#EBV variants
for i, j in df.iterrows():
    possible_key = '_(' + str(j['POS']) + ',' + str(j['REF']) + ',' + str(j['ALT']) + ')'
    num_var.add(possible_key)
    if possible_key not in variant_dict:
        variant_dict[possible_key] = []
        all_variants_set.add(possible_key)
    for sample in known_samples_set:
        if possible_key == "_(1432,G,A)":
        if j[sample] == "0" or j[sample] == 0:
            val = [str(sample.rstrip('_GT')), 0]
            variant_dict[possible_key].append(val)
        elif j[sample] == "1" or j[sample] == 1:
            val = [str(sample.rstrip('_GT')), 1]
            variant_dict[possible_key].append(val)

survival_data_file = open(sys.argv[3], 'r')
df_survival_data = pd.read_table(survival_data_file)

survival_dict = {}
sex_dict = {}
type_dict = {}
age_dict = {}
geography_dict = {}
tumor_resection_site_dict = {}
tumor_presentation_site_dict = {}
location_dict = {}
for i, j in df_survival_data.iterrows():
    if str(j['Study ID']) in raw_samples:
        survival_dict[j['Study ID']] = (j['Dead Status'],j['Survival'])
        if j['Sex'] == "M": #convert sex to binary (0/1) instead of M/F
            sex_dict[j['Study ID']] = 0
        elif j['Sex'] == "F":
            sex_dict[j['Study ID']] = 1
        type_dict[j['Study ID']] = j['EBV_Type']
        age_dict[j['Study ID']] = j['Age']
        tumor_resection_site_dict[j['Study ID']] = j['Tumor_Resection_Site']
        tumor_presentation_site_dict[j['Study ID']] = j['Tumor_Presentation_Site']
        geography_dict[j['Study ID']] = j['Geography']
        location_dict[j['Study ID']] = j['Location']

output_df_columns = []
output_df_columns.append('Study ID')
for variant in all_variants_set:
    output_df_columns.append(variant)

output_df_columns.append('Sex')
output_df_columns.append('EBV_Type')
output_df_columns.append('Age')
output_df_columns.append('Tumor_Resection_Site')
output_df_columns.append('Tumor_Presentation_Site')
output_df_columns.append('Geography')
output_df_columns.append('Location')
output_df_columns.append('Dead Status')
output_df_columns.append('Survival')

output_df = pd.DataFrame(columns = output_df_columns)

#for the BLGSP samples, mapping the file name (wgs_realn) to study ID (BLGSP-)
dataset1_fileid_to_snpid_file = open(sys.argv[4], 'r')
fileid_to_studyid = {}
studyid_to_fileid = {}
dataset1_fileid_to_snpid_file_lines = dataset1_fileid_to_snpid_file.readlines()
for line in dataset1_fileid_to_snpid_file_lines:
    array = line.rsplit('\t')
    studyid = array[0]
    fileid = array[1].rstrip('\n')
    fileid_to_studyid[fileid] = studyid
    studyid_to_fileid[studyid] = fileid

#filling in the dataframe output
for sample in raw_samples:
    value = {}
    for column in output_df_columns:
        substring = ""
        if column == 'Study ID':
            if sample in fileid_to_studyid:
                value['Study ID'] = fileid_to_studyid[sample]
            else:
                value['Study ID'] = sample
        elif column == 'Sex':
            value[column] = sex_dict[sample]
        elif column == 'EBV_Type':
            value[column] = type_dict[sample]
        elif column == 'Age':
            value[column] = age_dict[sample]
        #elif column in all_variants_set and column in variant_dict:
        elif column in all_variants_set:
            for possible_sample in variant_dict[column]:
                if sample == possible_sample[0]:
                    value[column] = possible_sample[1]
        elif column == "Tumor_Resection_Site":
            value[column] = tumor_resection_site_dict[sample]
        elif column == "Tumor_Presentation_Site":
            value[column] = tumor_presentation_site_dict[sample]
        elif column == "Geography":
            value[column] = geography_dict[sample]
        elif column == "Location":
            value[column] = location_dict[sample]
        elif column == "Dead Status" and sample in survival_dict:
            value[column] = survival_dict[sample][0]
        elif column == "Survival" and sample in survival_dict:
            value[column] = survival_dict[sample][1]

    #print(value)
    output_df = output_df.append(value, ignore_index=True)
output_df.to_csv('df_ebvvariants_231107.csv')
