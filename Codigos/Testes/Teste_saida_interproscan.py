import json

file_path = '/home/mantovani/Documents/trabalhos/Em_produção/IC/HMC_TE-s/Temp_Data/Dados_Dominos_Consevados/sequences.json'

with open(file_path, 'r') as file:

    dados = json.load(file)

    gene_1 = dados['results'][0]['openReadingFrames']

    for i in gene_1:

        print(i['protein']['matches'], '\n')

    file.close