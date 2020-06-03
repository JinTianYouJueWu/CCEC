import sys
import os
import pandas as pd


out_Cyton_file = sys.argv[1]

#Targetp部分
os.system(f"cat ./{out_Cyton_file}/1-Targetp/1-Targetp_process/*_summary.targetp2 | grep -v '#'|sort |uniq > ./{out_Cyton_file}/1-Targetp/all_targetp_result_1.txt")
os.system(f"cat ./{out_Cyton_file}/1-Targetp/all_targetp_result_1.txt| grep -v 'noTP' > ./{out_Cyton_file}/1-Targetp/all_targetp_result_2.txt")
file_out = open(f'./{out_Cyton_file}/1-Targetp/all_targetp_result.txt','w')
with open(f'./{out_Cyton_file}/1-Targetp/all_targetp_result_2.txt','r') as in_file :
    for line in in_file :
        line = line.strip().split()
        file_out.write(f"{line[0]}\t{line[1]}\n")
file_out.close()

os.system(f"rm ./{out_Cyton_file}/1-Targetp/all_targetp_result_1.txt")
os.system(f"rm ./{out_Cyton_file}/1-Targetp/all_targetp_result_2.txt")

#LOCALIZER 预测核质协同基因部分
os.system(f'cat ./{out_Cyton_file}/2-LOCALIZER/1-Localizer_process/*/chloroplast_predicted.fasta > ./{out_Cyton_file}/2-LOCALIZER/1-LOCALIZER_c.fasta')
os.system(f'cat ./{out_Cyton_file}/2-LOCALIZER/1-Localizer_process/*/mitochondria_predicted.fasta > ./{out_Cyton_file}/2-LOCALIZER/1-LOCALIZER_m.fasta')
os.system(f'cat ./{out_Cyton_file}/2-LOCALIZER/1-Localizer_process/*/chloroplast_mitochondria_predicted.fasta > ./{out_Cyton_file}/2-LOCALIZER/1-LOCALIZER_Dual.fasta')
os.system(f'cat ./{out_Cyton_file}/2-LOCALIZER/1-Localizer_process/*/mitochondria_chloroplast_predicted.fasta >> ./{out_Cyton_file}/2-LOCALIZER/1-LOCALIZER_Dual.fasta')

file_out = open(f'./{out_Cyton_file}/2-LOCALIZER/Localizer.result','w')

with open(f'./{out_Cyton_file}/2-LOCALIZER/1-LOCALIZER_c.fasta','r')  as in_file :
    for line in in_file :
        if line[0] == ">" :
           line = line.strip().split()
           gene = line[0][1:]
           file_out.write(f'{gene}\tPlastid\n')
with open(f'./{out_Cyton_file}/2-LOCALIZER/1-LOCALIZER_m.fasta','r')  as in_file :
    for line in in_file :
        if line[0] == ">" :
            line = line.strip().split()
            gene = line[0][1:]
            file_out.write(f'{gene}\tMitochondria\n') 
with open(f'./{out_Cyton_file}/2-LOCALIZER/1-LOCALIZER_Dual.fasta','r')  as in_file :
    for line in in_file :
        if line[0] == ">" :
            line = line.strip().split()
            gene = line[0][1:]
            file_out.write(f'{gene}\tDual\n') 

file_out.close()
os.system(f'cat ./{out_Cyton_file}/2-LOCALIZER/Localizer.result |sort |uniq > ./{out_Cyton_file}/2-LOCALIZER/Localizer_1.result')
os.system(f'rm ./{out_Cyton_file}/2-LOCALIZER/Localizer.result')

#CyMIRA
headers = ["gene", "Arabidopsis_thaliana", "length_rate", "pident"]
CyMIRA_data = pd.read_excel("CyMIRA.xlsx")
out_list = ['gene','CyMIRA targeting',
            'CyMIRA Interaction','CyMIRA Interaction Category','CyMIRA Interaction Subcategory']


Gene_data = pd.read_table(f'./{out_Cyton_file}/3-CyMIRA/1-blast/2-CyMIRA_blast_filter_1.txt' ,names=headers)
Gene_data['gene'] = Gene_data['gene'].astype(str)
Gene_data['Arabidopsis_thaliana_gene'] = Gene_data['Arabidopsis_thaliana'].apply(lambda x: x.split('.')[0])
Gene_data['Result'] = Gene_data['length_rate']*Gene_data['pident']
Gene_data_1 = Gene_data.groupby(by=['gene'])['Result'].max().to_frame()
Gene_data_1['Gene_1'] = Gene_data_1.index
Gene_data_result = pd.merge(Gene_data_1, Gene_data, how='left', left_on='Result', right_on='Result')
Gene_data_result = pd.merge(Gene_data_result, CyMIRA_data, how='left', left_on='Arabidopsis_thaliana_gene', right_on='AGI Identifier')
Gene_data_result = Gene_data_result[out_list].drop_duplicates()
Gene_data_result.to_excel(f"./{out_Cyton_file}/3-CyMIRA/Protein_CyMIRA_Result.xlsx", index=False)

#结果整合
data_Targetp = pd.read_table(f"./{out_Cyton_file}/1-Targetp/all_targetp_result.txt", names=['gene','Targetp targeting'])
data_LOCALIZER = pd.read_table(f"./{out_Cyton_file}/2-LOCALIZER/Localizer_1.result", names=['gene','LOCALIZER targeting'])
data_CyMRIA = pd.read_excel(f"./{out_Cyton_file}/3-CyMIRA/Protein_CyMIRA_Result.xlsx")

data_CyMRIA_Targetp = pd.merge(data_CyMRIA, data_Targetp, how='outer', left_on='gene', right_on = 'gene')
data_CyMRIA_Targetp_LOCALIZER = pd.merge(data_CyMRIA_Targetp, data_LOCALIZER, how='outer', left_on='gene', right_on = 'gene')
data_CyMRIA_Targetp_LOCALIZER.fillna("None",inplace=True)
data_CyMRIA_Targetp_LOCALIZER.to_excel(f"./{out_Cyton_file}/{out_Cyton_file}-Cytonuclear_Protein.xlsx",index=False)
