import requests 
import subprocess
from Bio.Seq       import Seq
from Bio.SeqRecord import SeqRecord
from Bio           import SeqIO
import pandas as pd
import numpy  as np

#LEMBRE-SE DE ORDENAR A LISTA DE FORMA QIUE A CLASSIFICAÇÃO CORRETA SEJA VIAVEL DESSE JEITO, VAI FUNCIONAR PARA TESTES MAS NA 
#HORA DE TREINAR O MODELO, ELE SIMPLESMENTE NÃO FUNCIONARÁ!

todas_hierarquias:tuple = ('Transib', 'Mariner', 'Zator', 'L1', 'RTE-X', '5S-RNA_Promoter', 'Group-I', 'Hero', 'scRNA', 'Poseidon', 
    'MuDR', 'Mutator-like', 'ISL2EU', 'Tad1_End', 'PiggyBac-X', 'DNA_Polymerase', 'Hydra-specific_Branch', 
    'Academ-2', 'Pararetroviridae', 'CR1-end', 'Odin', 'P_Element', 'tRNA_and_7SL_RNA', 'Stowaway', 'L1-group', 
    'Fot1', 'Restless', 'Retrotransposon', 'Casposon', 'CRE-1', 'Crypton-I', 'Orthoretrovirinae', 'Tc4', 
    'Unknown_LINE-dependent', 'Hydra', 'Charlie', 'ERV2-group', 'Viper-group', 'R4-like', 'Tc1-Mariner', 
    'RTE-end', 'R1-end', 'R1-subgroup', 'DIRS', 'Kolobok', 'R1', 'ISRm11', 'Gypsy-ERV', 'Crypton-R', 'Zorro', 
    'RTE-group', 'Chapaev_group', 'Ty1-Copia', 'Academ-1', 'Sola', 'Tyrosine_Recombinase', 'Q', 'Rex-end', 
    'Tyrosine_Recombinase_Elements', 'HarbS', 'Long_Terminal_Repeat_Element', 'Crypton-S', 'R1-group', 
    'Zisupton', 'LINE', 'hobo', 'BovB-end', 'Pseudogene', 'Mermaid', 'Penelope-like_Elements', 'snRNA', 
    'Tad1', 'Ginger-2', 'ERV2', 'IS3EU', 'CMC', 'SINE-like', 'ERV1', 'T2', 'Novosib', 'Transposable_Element', 
    'Crypton-V', 'Dong-R4', 'Deceiver', 'ERV3', 'R2', 'Alu', 'L2-derived', 'L1-derived', 'SINE', 'ERV4', 'Ant1', 
    'V_and_MIR-core', 'Group-II', 'Rex-Babar', 'Ceph-core', 'Sola-3', 'Class_II_DNA_Transposition', 'Crypton', 
    'Neptune', 'Helicase', 'CACTA', 'CRE-2', 'CR1', 'Coprina', 'Athena', 'Blackjack', 'TRIM', 'm44', 'Mirage', 
    'Caulimoviridae', 'Academ-H', 'LOA', 'NOF', 'Tx1', 'CR1-group', 'I-end', 'Bel-Pao', 'Maverick', 'Nematis', 
    'Transposase', 'Crypton-A', 'ORTE', 'Ginger', 'Tag1', 'R2-like', 'Merlin', 'PiggyBac', 'Kolobok-H', 
    'I-group', 'Dada', 'Tc2', 'Tc2-group', 'Pegasus', 'Gypsy', 'tRNA', 'hATx', 'Class_I_Retrotransposition', 
    'Ambal', 'Ricksha', 'Jockey', 'RTE', 'BovB', 'Chlamys', 'Meta-core', 'Proto-1', 'hATw', 'RNA', 'hAT1', 
    'L1-dependent', 'Genie', 'Lenti', 'I', 'L2-group', 'Sauria-core', 'EnSpm', 'Harbinger', 'No-core', 
    'Lacking_Small_RNA_pol_III_Promoter', 'SVA', 'I-derived', 'CRE', 'Tigger', 'Proto-2', 'Cweed', 'DRE', 
    'Sola-1', 'R2-end', 'hATm', 'R4-derived', 'Helitron-1', 'Crypton-C', 'Unknown_Promoter', 
    'Interspersed_Repeat', 'Crypton-F', 'Crypton-H', 'RTE-derived', 'No_or_Unknown_Core', 'L2', 
    'Centromeric', 'L2-end', 'Group-2', 'RTE-like', 'PiggyBac_Group', 'Chapaev', 'Fungi-specific_Branch', 
    'MaLR', 'IS885', 'V-core', 'Gizmo', 'Unknown', 'Tip100', 'Spy', 'Academ', 'Activator', 'Jockey-end', 
    'tRNA_and_5S_RNA', 'hAT5', 'TATE', 'Tc1', 'Zenon', 'F', 'Mogwai', '7SL-RNA_Promoter', 'PIF-Harbinger', 
    'snoRNA', 'Chapaev-3', 'L1-like', 'Helitron-2', 'Ginger-1', 'Kolobok-E', 'hAT19', 'B2', 'rRNA', 
    'U-RNA_Promoter', 'Dualen', 'Maverick-Mavirus', 'Crypton-X', 'NeSL', 'Ngaro', 'Sagan', 'MIR-core', 
    'Sola-2', 'Retroviridae', 'Naiad', 'hAT', 'PiggyBac-A', 'LINE-dependent_Retroposon', 'Pogo', 'hAT6', 
    'tRNA_Promoter', 'R1-like', 'Viper', 'Group-1', 'Spumaretrovirinae', 'Deu-core')

def lista_binaria(classificacao_H: str) -> tuple:
    
    classificacao_H = classificacao_H.split(';')

    return (1 if i in classificacao_H else 0 for i in todas_hierarquias)

def reverso_lista_binaria(lista_binaria: tuple) -> list:

    classificacao:list = [j for i, j in zip(lista_binaria, todas_hierarquias) if i]

    return classificacao

class Data:

    def __init__(self, relative_path:str, taxa_dowload: int = 1000) -> None: #colocar um pedaço de código para calcular o caminho relativo sozinho
        
        self.relative_path = relative_path
        self.url:          str  = 'https://dfam.org/api/families'
        self.taxa_dowload: int  = taxa_dowload
        self.params:       dict = {"format": "full", 
                                   "start" : 0, 
                                   "limit" : self.taxa_dowload}
        return None
    
    def __consulta_Dfan(self) -> pd.DataFrame:

        resposta: requests = requests.get(self.url, self.params)

        try:

            resultados: dict = resposta.json()['results']
            
        except:

            raise ValueError("O formato da entrada não é compativel com jason, veja se eles não mudaram a API.")
        
        saida: pd.DataFrame = pd.DataFrame(resultados)[['accession', 'classification', 'consensus_sequence']]
        saida['classification']     = saida['classification'].apply(lista_binaria)
        saida['consensus_sequence'] = saida['consensus_sequence'].apply(Seq)

        self.params["start"] += self.taxa_dowload
        
        return saida
    
    def __fasta_file(self, data_frame: pd.DataFrame) -> None:

        file_path = self.relative_path + 'Temp_Data/Arquivos_Fasta/sequences'
        
        a = [SeqRecord(row.consensus_sequence, id=row.accession) for row in data_frame.itertuples(index=False)]

        SeqIO.write(a, file_path, "fasta")

        return None

    def __conserved_domain(self) -> None:
        #PRINTS, Pfam, CDD, ProSitePatterns

        comando_0 = './interproscan.sh ' + ' -appl PRINTS' +' -t n '
        comando_1 = ' -i ' + self.relative_path + 'Temp_Data/Arquivos_Fasta/sequences' + ' -b ' + self.relative_path + 'Temp_Data/Dados_Dominios_Conservados/'

        comando_f = comando_0 + comando_1

        subprocess.run(comando_f, cwd='/home/mantovani/interproscan-5.69-101.0', shell=True)

        return None
    
    def get_data(self) -> None:

        def __uma_iteracao() -> None:

            dados = self.__consulta_Dfan()

            self.__fasta_file(dados)
            self.__conserved_domain()

            return None
        
        __uma_iteracao()

        return None
    
def main() -> None:

    data = Data('/home/mantovani/Documents/IC/HMC_TE-s/', 10)
    
    
    a = data.get_data()

    return None

 
if __name__ == '__main__':

    main()