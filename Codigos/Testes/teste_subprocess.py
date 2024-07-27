import subprocess

relative_path = '/home/mantovani/Documents/trabalhos/Em_produção/IC/HMC_TE-s/'

comando_0 = './interproscan.sh ' + ' -appl PRINTS, PANTHER, Pfam ' +' -t n ' 
comando_1 = ' -i ' + relative_path + 'Temp_Data/Arquivos_Fasta/sequences' + ' -b ' + relative_path + 'Temp_Data/Dados_Dominos_Consevados/'

comando_f = comando_0 + comando_1
print(comando_f)

subprocess.run(comando_f, cwd='/home/mantovani/interproscan-5.68-100.0', shell=True)