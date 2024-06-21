#!/opt/anaconda3/envs/genomic-tools/bin/python3

### --- SET IMPORTS --- ###
from bioservices.kegg import KEGG
from UCSC import hg38
import gene-lookup
import json, wget, argparse, os, csv, shutil
import pandas as pd 

### --- SET ARGUMENTS --- ###
parser = argparse.ArgumentParser(description='case control analysis based on KEGG human disease pathways')
parser.add_argument('-r', '--reset', action='store_true', help='redownload kegg disease database json')
parser.add_argument('-c', '--csv', action='store_true', help='begin from gtf.csv file')
args = parser.parse_args()

### --- SET RESET DATA DOWNLOAD SWITCH --- ###
if args.reset:
    dir = os.path.expanduser('~/db')
    print('downloading kegg human disease database...')
    # download kegg human disease json database
    os.chdir(dir)
    url = 'https://www.genome.jp/kegg-bin/download_htext?htext=br08402.keg&format=json&filedir=' # may need to update this link in the future
    wget.download(url, 'kegg-disease.json')
    hg38.get_gtf() # see if this downloads gtf file 
else:
    pass

### --- PARSE KEGG JSON INTO PANDAS DF --- ###
filter_dict = {}
def parse_kegg():
    dir = os.path.expanduser('~/db')
    f = open(dir + '/kegg-disease.json', 'r')
    j = json.load(f)
    data = {} # dict formatted for pandas
    database = j["children"]  
    iteration = 1
    key = 1
    print('DISEASE CLASSIFICATION OPTIONS: ')
    for classification in database:
            # disease classifications
            classif = classification['name']

            for subclass in classification['children']:
                # disease subclassification
                subclassif = subclass['name']
                # print out subclassifications
                print(str(key) + ' - ' + subclassif)
                filter_dict[key] = subclassif
                key = key + 1

                for disease in subclass['children']:
                    # disease name and kegg identifyer
                    kegg_id = disease['name'].split('  ')[0]
                    disease_name = disease['name'].split('  ')[1]
                    
                    # parse out disease name from pathway id
                    if '[' in disease_name:
                        disease_name = disease['name'].split('  ')[1].split('[')[0]
                        pathway_id = disease['name'].split('  ')[1].split('[')[1][:-1].split(':')[1]

                        # parse double pathway ids
                        if len(pathway_id) > 8:
                            path_list = disease['name'].split('  ')[1].split('[')[1][:-1].split(':')[1].split(' ')
                            path_ids = []
                            for item in path_list:
                                path_ids.append(item)
                            pathway_id = path_ids
                        else:
                            pathway_id = [pathway_id]
                    elif '[' not in disease_name:
                        pathway_id = ''
                    # now create dictionary for pandas df
                    data[iteration] = [disease_name, kegg_id, classif, subclassif, pathway_id]
                    iteration = iteration + 1
    df = pd.DataFrame.from_dict(data, orient='index', columns = ['Name', 'Kegg_ID', 'Classification', 'Subclassification', 'Pathway_ID'])
    return df

### --- DEFINITE DISEASE FILTER AND QUERY KEGG BIOSERVICES API --- ###
id_dict = {} # dictionary for disease names
def set_filter(dataframe): # eventually this will need to work on keywords
    # list of genes to return
    s = KEGG()
    genes = []
    
    # user input
    value = input("Enter Number of Disease Classification:\n")
    print(f'Pulling {filter_dict[int(value)]} associated genes from KEGG...')
    
    # filter dataframe 
    is_disease = dataframe['Subclassification'] == filter_dict[int(value)]
    df_disease =  dataframe[is_disease]
    
    # loop through ids
    for id in df_disease['Kegg_ID'].to_list(): 
        entry = df_disease.loc[df_disease['Kegg_ID'] == id].squeeze()
        ailment = entry['Name']
        id_dict[id] = ailment # add to id dictionary 

        # Begin storing genes 
        data = s.get(id)
        dict_data = s.parse(data)
        # print(id)
        # only take entries with genes
        if 'GENE' in dict_data:
            # loop through each gene:
            for symbol in dict_data['GENE']:
                # only take all caps gene symbols
                if symbol.strip('(').strip(')').isupper() == True:
                    genes.append(symbol.strip('(').strip(')').strip(';'))
                else:
                    pass         
        else:
            pass
    return list(set(genes))


### --- CREATE HG38 BED FILES FOR EACH GENE --- ###
def extract_bed(genelist):
    
    dir = os.path.expanduser('~/db')
    if os.path.exists('tmp/'):
        shutil.rmtree('tmp/')
        os.makedirs('tmp/')
    else:
        os.makedirs('tmp/')
    out_object = open('tmp/tmp.bed', 'w')
    out_writer = csv.writer(out_object, delimiter='\t')


    ### --- open HGNC and OMIM as dataframe --- ###
    hgnc = gene_lookup.genetable.hgnc()
    omim = gene_lookup.genetable.omim()
    
    ### --- Set switch for reading gtf-dataframe from previously written csv --- ###
    if args.reset == False:
        print('Loading dataframe from gtf.csv...') 
        data = pd.read_csv(dir + '/hg38.ncbiRefSeq.csv') # load from pre-written csv 
    else:
        data = hg38.parse_gtf(dir + '/hg38.ncbiRefSeq.gtf') # parse from re-downloaded csv
        out = open(dir + '/hg38.ncbiRefSeq.csv', 'w')
        data.to_csv(out, index=False)

    ### --- lookup gene names in gtf df --- ###
    print('Creating interval bed file...')
    for symbol in genelist:     
        entry = data.loc[(data['feature'] == 'transcript') & ((data['gene_id'] == symbol) | (data['gene_name'] == symbol))].squeeze()

        ### --- SEARCH HGNC AND OMIM FOR MISSING SYMBOLS--- ###
        if entry.empty:
            print('NOTE: ' + symbol + ' is absent from gtf... searching HGNC')
            matches = []
            
            ### --- SEARCH HGNC FOR EQUIVALENT SYMBOL --- ###
            df = gene_lookup.lookup(hgnc)
            result1 = df.match_hgnc(symbol)
            
            if result1 == False:
                ### --- SEARCH OMIM FOR EQUIVALENT SYMBOL --- ###
                print('Searching OMIM database for ' + symbol + ' synonyms...')
                df2 = gene_lookup.lookup(omim)
                result2 = df2.match_omim(symbol)

                ### --- REPORT NO MATCH PRESENT --- ###
                if result2 == False:
                    print('NOTE: no equivelent symbols found for ' + symbol)
                    pass

                ### --- CREATE BED INTERVAL FOR LONGEST TRANSCRIPT --- ###
                else:
                    for syn in result2:
                        out = hg38.make_bed(data, syn)
                        if out == False:
                            pass
                        else:
                            print('Located ' + syn + ' in hg38 GTF: ' + str(out))
                            out_writer.writerow(out)
            
            ### --- CREATE BED INTERVAL FOR LONGEST TRANSCRIPT  --- ###
            else:
                for syn in result1:
                    out = hg38.make_bed(data, syn)
                    if out == False:
                        pass
                    else:
                        print('Located ' + syn + ' in hg38 GTF: ' + str(out))
                        out_writer.writerow(out)

        ### --- CREATE BED FROM MATCHES IN GTF ---###
        else: 
            out = hg38.make_bed(data,symbol)
            if out == False:
                pass
            else:
                print('Located ' + symbol + ' in hg38 GTF: ' + str(out))
                out_writer.writerow(out)
            

        
        #print(entry)



###--- BEGIN MAIN EXECUTION ---###
df = parse_kegg()
genes = set_filter(df)
extract_bed(genes)
