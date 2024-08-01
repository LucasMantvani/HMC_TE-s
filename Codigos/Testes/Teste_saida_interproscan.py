import json
import pandas as pd

file_path = '/home/mantovani/Documents/trabalhos/Em_produção/IC/HMC_TE-s/Temp_Data/Dados_Dominos_Consevados/sequences.json'

with open(file_path, 'r') as file:

    dados = json.load(file)

    for i in dados['results']:

        gene_1 = i['openReadingFrames']

        df = pd.DataFrame([{
        'id_DC': item['id'],
        'start': item['start'],
        'end': item['end'],
        'strand': item['strand'],
        'protein_sequence': item['protein']['sequence'],
        } for item in gene_1])

            # print(j['protein']['matches'], '\n')
        print(df, '\n')

    file.close