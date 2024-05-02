import pandas as pd
import sys

blgsp_coding_ssms_file = open(sys.argv[1], 'r')
df_blgsp_coding_ssms = pd.read_csv(blgsp_coding_ssms_file)

sample_file = open(sys.argv[2], 'r')
sample_file_lines = sample_file.readlines()

dataset1_studyids_file = open("/Users/isaacekimjr/Desktop/dataset1_studyids.txt", 'r')
dataset1_studyids_file_lines = dataset1_studyids_file.readlines()
raw_samples = set()
for line in dataset1_studyids_file_lines:
    raw_samples.add(line.rstrip('\n'))

sample_to_gene_dict = {}
genes_dict = {}

dataset1_samples = set()

samples_with_human_tumor_mutated_genes = set()

for i, j in df_blgsp_coding_ssms.iterrows():

    if "IGH" in str(j['Hugo_Symbol']):
        gene_name = "IGH"
    elif "SIGL" not in str(j['Hugo_Symbol']) and "IGLON" not in str(j['Hugo_Symbol']) and "IGL" in str(j['Hugo_Symbol']):
        gene_name = "IGL"
    elif "PIGK" not in str(j['Hugo_Symbol']) and "IGK" in str(j['Hugo_Symbol']):
        gene_name = "IGK"
    else:
        gene_name = j['Hugo_Symbol']

    if gene_name not in genes_dict:
        genes_dict[gene_name] = set()

    specimen = str(j['tumor_biospecimen_id'])
    sample_name = ""
    isSample = False
    for sample in raw_samples:
        if str(sample) in specimen: #if given sample has survival data + EBV(+)
            samples_with_human_tumor_mutated_genes.add(sample)
            sample_name = str(sample)
            isSample = True
            dataset1_samples.add(sample)
            break
    if isSample:
        if sample_name in sample_to_gene_dict:
            if gene_name in sample_to_gene_dict[sample_name]:
                if j['Variant_Classification'] == 'Missense_Mutation':
                    sample_to_gene_dict[sample_name][gene_name][0] += 1
                    genes_dict[gene_name].add(sample_name)
                elif j['Variant_Classification'] == 'Frame_Shift_Del' or j['Variant_Classification'] == 'Frame_Shift_Ins':
                    sample_to_gene_dict[sample_name][gene_name][1] += 1
                    genes_dict[gene_name].add(sample_name)
                elif j['Variant_Classification'] == 'Nonsense_Mutation' or j['Variant_Classification'] == 'Nonstop_Mutation':
                    sample_to_gene_dict[sample_name][gene_name][2] += 1
                    genes_dict[gene_name].add(sample_name)
                elif j['Variant_Classification'] == 'Splice_Region' or j['Variant_Classification'] == 'Splice_Site':
                    sample_to_gene_dict[sample_name][gene_name][3] += 1
                    genes_dict[gene_name].add(sample_name)
                elif j['Variant_Classification'] == 'In_Frame_Ins' or j['Variant_Classification'] == 'In_Frame_Del'or j['Variant_Classification'] == 'Targeted_Region':
                    sample_to_gene_dict[sample_name][gene_name][4] += 1
                    genes_dict[gene_name].add(sample_name)
            else:
                sample_to_gene_dict[sample_name][gene_name] = [0,0,0,0,0]
                if j['Variant_Classification'] == 'Missense_Mutation':
                    sample_to_gene_dict[sample_name][gene_name][0] += 1
                    genes_dict[gene_name].add(sample_name)
                elif j['Variant_Classification'] == 'Frame_Shift_Del' or j['Variant_Classification'] == 'Frame_Shift_Ins':
                    sample_to_gene_dict[sample_name][gene_name][1] += 1
                    genes_dict[gene_name].add(sample_name)
                elif j['Variant_Classification'] == 'Nonsense_Mutation' or j['Variant_Classification'] == 'Nonstop_Mutation':
                    sample_to_gene_dict[sample_name][gene_name][2] += 1
                    genes_dict[gene_name].add(sample_name)
                elif j['Variant_Classification'] == 'Splice_Region' or j['Variant_Classification'] == 'Splice_Site':
                    sample_to_gene_dict[sample_name][gene_name][3] += 1
                    genes_dict[gene_name].add(sample_name)
                elif j['Variant_Classification'] == 'In_Frame_Ins' or j['Variant_Classification'] == 'In_Frame_Del'or j['Variant_Classification'] == 'Targeted_Region':
                    sample_to_gene_dict[sample_name][gene_name][4] += 1
                    genes_dict[gene_name].add(sample_name)
        else:
            sample_to_gene_dict[sample_name] = dict()
            sample_to_gene_dict[sample_name][gene_name] = [0,0,0,0,0]
            if j['Variant_Classification'] == 'Missense_Mutation':
                sample_to_gene_dict[sample_name][gene_name][0] += 1
                genes_dict[gene_name].add(sample_name)
            elif j['Variant_Classification'] == 'Frame_Shift_Del' or j['Variant_Classification'] == 'Frame_Shift_Ins':
                sample_to_gene_dict[sample_name][gene_name][1] += 1
                genes_dict[gene_name].add(sample_name)
            elif j['Variant_Classification'] == 'Nonsense_Mutation' or j['Variant_Classification'] == 'Nonstop_Mutation':
                sample_to_gene_dict[sample_name][gene_name][2] += 1
                genes_dict[gene_name].add(sample_name)
            elif j['Variant_Classification'] == 'Splice_Region' or j['Variant_Classification'] == 'Splice_Site':
                sample_to_gene_dict[sample_name][gene_name][3] += 1
                genes_dict[gene_name].add(sample_name)
            elif j['Variant_Classification'] == 'In_Frame_Ins' or j['Variant_Classification'] == 'In_Frame_Del'or j['Variant_Classification'] == 'Targeted_Region':
                sample_to_gene_dict[sample_name][gene_name][4] += 1
                genes_dict[gene_name].add(sample_name)

#Sandeep's variants
sandeep_somatic_variants_file = open(sys.argv[3], 'r')
df_sandeep_somatic_variants = pd.read_csv(sandeep_somatic_variants_file)
for i,j in df_sandeep_somatic_variants.iterrows():
    if "IGH" in str(j['Gene.refGene']):
        gene_name = "IGH"
    elif "SIGL" not in str(j['Gene.refGene']) and "IGLON" not in str(j['Gene.refGene']) and "IGL" in str(j['Gene.refGene']):
        gene_name = "IGL"
    elif "PIGK" not in str(j['Gene.refGene']) and "IGK" in str(j['Gene.refGene']):
        gene_name = "IGK"
    else:
        gene_name = j['Gene.refGene']

    if gene_name not in genes_dict:
        genes_dict[gene_name] = set()
    for sample_name in raw_samples:
        if sample_name in j:
            dataset1_samples.add(sample_name)
            if int(j[sample_name]) == 1:
                if sample_name in sample_to_gene_dict:
                    if gene_name in sample_to_gene_dict[sample_name]:
                        if j['ExonicFunc.refGene'] == 'nonsynonymous_SNV':
                            sample_to_gene_dict[sample_name][gene_name][0] += 1
                            genes_dict[gene_name].add(sample_name)
                        elif j['ExonicFunc.refGene'] == 'frameshift_deletion' or j['ExonicFunc.refGene'] == 'frameshift_insertion':
                            sample_to_gene_dict[sample_name][gene_name][1] += 1
                            genes_dict[gene_name].add(sample_name)
                        elif j['ExonicFunc.refGene'] == 'stopgain' or j['ExonicFunc.refGene'] == 'stoploss':
                            sample_to_gene_dict[sample_name][gene_name][2] += 1
                            genes_dict[gene_name].add(sample_name)
                        elif j['Func.refGene'] == 'splicing':
                            sample_to_gene_dict[sample_name][gene_name][3] += 1
                            genes_dict[gene_name].add(sample_name)
                    else:
                        sample_to_gene_dict[sample_name][gene_name] = [0,0,0,0,0]
                        if j['ExonicFunc.refGene'] == 'nonsynonymous_SNV':
                            sample_to_gene_dict[sample_name][gene_name][0] += 1
                            genes_dict[gene_name].add(sample_name)
                        elif j['ExonicFunc.refGene'] == 'frameshift_deletion' or j['ExonicFunc.refGene'] == 'frameshift_insertion':
                            sample_to_gene_dict[sample_name][gene_name][1] += 1
                            genes_dict[gene_name].add(sample_name)
                        elif j['ExonicFunc.refGene'] == 'stopgain' or j['ExonicFunc.refGene'] == 'stoploss':
                            sample_to_gene_dict[sample_name][gene_name][2] += 1
                            genes_dict[gene_name].add(sample_name)
                        elif j['Func.refGene'] == 'splicing':
                            sample_to_gene_dict[sample_name][gene_name][3] += 1
                            genes_dict[gene_name].add(sample_name)
                else:
                    sample_to_gene_dict[sample_name] = dict()
                    sample_to_gene_dict[sample_name][gene_name] = [0,0,0,0,0]
                    if j['ExonicFunc.refGene'] == 'nonsynonymous_SNV':
                        sample_to_gene_dict[sample_name][gene_name][0] += 1
                        genes_dict[gene_name].add(sample_name)
                    elif j['ExonicFunc.refGene'] == 'frameshift_deletion' or j['ExonicFunc.refGene'] == 'frameshift_insertion':
                        sample_to_gene_dict[sample_name][gene_name][1] += 1
                        genes_dict[gene_name].add(sample_name)
                    elif j['ExonicFunc.refGene'] == 'stopgain' or j['ExonicFunc.refGene'] == 'stoploss':
                        sample_to_gene_dict[sample_name][gene_name][2] += 1
                        genes_dict[gene_name].add(sample_name)
                    elif j['Func.refGene'] == 'splicing':
                        sample_to_gene_dict[sample_name][gene_name][3] += 1
                        genes_dict[gene_name].add(sample_name)

output_df_columns = []
output_df_columns.append('Study ID')

gene_mutation_count_dict = {}
for gene in genes_dict:
    gene_mutation_count_dict[gene] = len(genes_dict[gene])

sorted_genes_dict = dict(sorted(gene_mutation_count_dict.items(), key=lambda x:x[1], reverse=True))

output_df = pd.DataFrame(columns = output_df_columns)

gene_count_dict = {}

num_samples = 0
for sample in sample_to_gene_dict:
    value = {}
    value['Study ID'] = sample
    for gene in sorted_genes_dict:
        if gene not in gene_count_dict:
            gene_count_dict[gene] = 0
        if gene in sample_to_gene_dict[sample]:
            nonsynonymous_n = sample_to_gene_dict[sample][gene][0]
            frameshift_n = sample_to_gene_dict[sample][gene][1]
            stopgain_n = sample_to_gene_dict[sample][gene][2]
            splicing_n = sample_to_gene_dict[sample][gene][3]
            nonframeshift_n = sample_to_gene_dict[sample][gene][4]
            if nonsynonymous_n + frameshift_n + stopgain_n + splicing_n + nonframeshift_n > 0:
                value[gene] = 1
                gene_count_dict[gene] += 1
            else:
                value[gene] = 0
        else:
            value[gene] = 0

    num_samples += 1
    # print(yo)
    # print("hello")
    # print("hello")
    # print("hello")
    # print(value)

for gene in sorted_genes_dict:
    if gene_count_dict[gene] > 1:
        output_df_columns.append(gene)

#print(gene_count_dict)

igh_list = []
igl_list = []
igk_list = []
for gene in gene_count_dict:
    if "IGH" in gene:
        igh_list.append(gene)
    if "SIGL" not in gene and "IGLON" not in gene and "IGL" in gene:
        igl_list.append(gene)
    if "PIGK" not in gene and "IGK" in gene:
        igk_list.append(gene)

for sample in sample_to_gene_dict:
    value = {}
    value['Study ID'] = sample
    for gene in sorted_genes_dict:
        if gene_count_dict[gene] > 1: #only have columns/variants/human tumor genes that are mutated in more than 6 sample samples (aka > 5%)
            if gene in sample_to_gene_dict[sample]:
                nonsynonymous_n = sample_to_gene_dict[sample][gene][0]
                frameshift_n = sample_to_gene_dict[sample][gene][1]
                stopgain_n = sample_to_gene_dict[sample][gene][2]
                splicing_n = sample_to_gene_dict[sample][gene][3]
                nonframeshift_n = sample_to_gene_dict[sample][gene][4]
                if nonsynonymous_n + frameshift_n + stopgain_n + splicing_n + nonframeshift_n > 0:
                    value[gene] = 1
                else:
                    value[gene] = 0
            else:
                value[gene] = 0
    output_df = output_df.append(value, ignore_index=True)

# print(sorted_genes_dict)

for item in igk_list:
    print(item)

print(gene_count_dict['IGH'])
print(gene_count_dict['IGL'])
print(gene_count_dict['IGK'])
# print(num_samples)

print(gene_count_dict)
unique_genes = set()
for key in gene_count_dict:
    unique_genes.add(key)

print(len(unique_genes))

output_df.to_csv('df_human_tumor_mutated_genes_231115_240225trial.csv')
