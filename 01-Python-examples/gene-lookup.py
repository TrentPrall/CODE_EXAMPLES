#!/opt/anaconda3/envs/genomic-tools/bin/python3
import pandas as pd 
import os
import wget


class genetable:
    
    # IMPORTANT: need to add download and parser functions to OMIM download
    def omim():
        print('Downloading OMIM database...')
        omim = pd.read_csv('~/db/full_omim_table.txt', dtype='str', engine='c', delimiter='\t')
        return omim

    def hgnc():
        base = os.getcwd()
        if os.path.exists('~/db/hgnc_complete_set.txt'):
            os.remove('~/db/hgnc_complete_set.txt')
        else:
            pass
        print('Downloading HGNC database...')
        dir = os.path.expanduser('~/db')
        os.chdir(dir)
        wget.download('ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt')
        hgnc = pd.read_csv('hgnc_complete_set.txt', dtype='str', engine='c', delimiter= '\t')
        os.chdir(base)
        return hgnc


class lookup:   
    def __init__(self, dataframe):
        self.df = dataframe
        
    def match_omim(self, symbol):
        df2 = self.df[self.df.apply(lambda row: row.astype(str).str.contains(symbol).any(), axis=1)][['genes', 'hgnc_synonyms']]
        df2 = df2.dropna()
        

        ###--- return false for no matches in dataframe ---###
        check = df2.squeeze()
        if check.empty:
            return False
        else:
            match = []


            ###--- generate list of symbols per row ---###
            for row in range(df2.shape[0]):
                index = df2.index[row]
                potentials = []
                for col in range(df2.shape[1]):
                    # print(df2.iat[row,col])
                    if '|' in str(df2.iat[row,col]):
                        potentials.extend(str(df2.iat[row,col]).split('|'))
                    else:
                        potentials.append(df2.iat[row,col])
                potentials = list(set(potentials))


                ###--- iterate through list of symbols per row ---###
                for sym in potentials:    
                    if sym == symbol:
                        matches = self.df.at[index,'hgnc_genes'] # pull from original dataframe
                        match.extend(str(matches).split(','))
                        break
                    else:
                        pass
        

        ###--- Check if there are any perfect matches ---###
        if len(match) == 0:
            print('No equivelent hgnc symbols in OMIM database for ' + symbol)
            return False
        else:
            print('Found ' + symbol + ' in OMIM database, converted to HGNC symbols: ' + str(list(set(match))))
            match = list(set(match))
            return match





        # now iterate through each row and column to find a perfect match

    def match_hgnc(self, symbol):
        df2 = self.df[self.df.apply(lambda row: row.astype(str).str.contains(symbol).any(), axis=1)]['alias_symbol']
        df2 = df2.dropna() # drop nan from df2
        
        if df2.empty:
            print('No synonyms found for ' + symbol + ' in HGNC')
            return False
        else:
            match = []
            
            ###--- create list of synonyms for each row with match ---###
            for row in df2.iteritems():
                potentials = []
                if '|' in row[1]:
                    # print(row[1].split('|'))
                    potentials.extend(row[1].split('|'))
                else:
                    # print(row[1])
                    potentials.append(row[1])
                potentials = list(set(potentials))

                
                ###--- iterate through list of synonyms ---###
                for syn in potentials:
                    if syn == symbol:
                        matches = self.df.at[row[0],'symbol']
                        print('Found ' + symbol + ' synonym in HGNC database: ' + matches)
                        match.append(matches)
                        break
                    else:
                        pass
                    
            ###--- Check if any direct matches occured ---###        
            if len(match) == 0:
                print('No synonyms found for ' + symbol + ' in HGNC')
                return False
            else:
                match = list(set(match))
                return match


###--- TEST METHODS ---###
# omim = genetable.omim()
# genetable = lookup(omim)
# match = genetable.match_omim('VAMP7')
# print(match)

# hgnc = genetable.hgnc()
# df = lookup(hgnc)
# match = df.match_hgnc('C')
# print(match)