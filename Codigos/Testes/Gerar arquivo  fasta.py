import pandas as pd 

sequences = {
    "seq1": "CAGTGGTTCTCAAACTATGCTCCCACGGAACACTGGTGTTCCGCGAGCTGGACCGAGGTGTTCCGCCGCTAAAATCGAATAATGGCAGTCTTTTCCCAATTCGCAGAAAAAATTAAGTTAATGGATAGAATTTTCTCAATTTTAAGTTTTTCTCCCATAAATTTTTCTTAAACCTTTGGCGCCTGCTACCGCTTGTCTGCTAGCGAGCGCTAGCAGATGACAAAACGAATGAACAAGTAGCAAATGAATATTAACGGAAAAACTCGGAAACAGCCGATCATGCTTAGTTCGACAGCTTGTTTTTTTCATGGACACGTGTGCTCGTGGTTGTGATTTTATATTGCATACATATAATATAGTCATATTTTCTGCTAGTAAAATTGTTCCTAGTTCGTATGATGGAATTAAATATATCAGTATTTCATGGTGTTCCGCAAAGAGCACCATGACTTCCAGGTGCTCTGCCACCTGAACAAGTTTGAGAACCACTG"
}



def file_fasta(sequencias: dict, file_path: str) -> None:


    with open(file_path, "w") as fasta_file:

        for seq_id, sequence in sequencias.items():

            fasta_file.write(f">{seq_id}\n")
            fasta_file.write(f"{sequence}\n")

    print(f"Arquivo {file_path} criado com sucesso!")

    return None



df = pd.read_csv('/home/mantovani/Documents/trabalhos/Em_produção/IC/HMC_TE-s/Temp_Data/Dados_Genes/dado.csv')

df = df[['accession', 'consensus_sequence']]

df        = df.iloc[range(0,99)]
dict_list = df.set_index('accession').to_dict()['consensus_sequence']

print(dict_list)

diretorio = "/home/mantovani/Documents/trabalhos/Em_produção/IC/HMC_TE-s/Temp_Data/Arquivos_Fasta/"

file_fasta(dict_list, diretorio + "sequences.fasta")
