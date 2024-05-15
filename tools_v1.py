
#project 3: miRNA-Disease Database & Network
#author: fangfei lin 
#ID Number: 11466444


#importing the pandas library
import pandas as pd


#printing a decorative header for the console output
print('#'*30)
print('## LOADING DATABASE')
print('#'*30)
print('\n')

#source data
source_data = pd.read_csv('clean_data.tsv',sep='\t')

print('#'*30)
print('## DATABASE LOADED')
print('#'*30)
print('\n')

#print(source_data.head(10).to_markdown())
'''
	miRNA	Locus	Score	Organism	Gene	geneSymbol	diseaseName	score
0	hsa-let-7a-2-3p	NM_004281	90.98065	Homo sapiens	BAG3	BAG3	CARDIOMYOPATHY, DILATED, 1HH	0.9
1	hsa-let-7g-3p	NM_004281	90.98065	Homo sapiens	BAG3	BAG3	CARDIOMYOPATHY, DILATED, 1HH	0.9
2	hsa-miR-129-1-3p	NM_004281	82.90620	Homo sapiens	BAG3	BAG3	CARDIOMYOPATHY, DILATED, 1HH	0.9
'''


#import the argparse module
import argparse

#initialize the argument parser with a description of the tool's functionality
parser=argparse.ArgumentParser(description="a simple local lookup database to look up diseases a given miRNA is associated to, or look up miRNAs a disease is associated to")

#add command line arguments for various functionalities of the tool
#the arguments allow the user to specify their query and filter preferences
parser.add_argument('--mirna', help='Input mirna.')
parser.add_argument('--disease', help='Input gene')
parser.add_argument('--browse', help='Given part of a miRNA name, it should report back all the miRNAs that match it', action="store_true")
parser.add_argument('--score', choices=['True', 'False'], help='Set "False" for without score.', default='True')
parser.add_argument('--gene_conf', type=int, help='If no confidence score is provided, all miRNA to Gene associations above 80 will be reported.', default=80)
parser.add_argument('--disease_conf', type=float, help='If no confidence score is provided, all Gene to Disease associations above 0.8 should be reported.', default=0.8)

#parse the arguments provided by the user
args = parser.parse_args()

#assign parsed arguments to variables for easier access
mirna = args.mirna
disease = args.disease
browse = args.browse
score = args.score
gene_conf = args.gene_conf
disease_conf = args.disease_conf

#print the current query parameters for confirmation
print('Current query parameters:')
print(f'mirna: {mirna}')
print(f'disease: {disease}')
print(f'browse: {browse}')
print(f'score: {score}')
print(f'gene_conf: {gene_conf}')
print(f'disease_conf: {disease_conf}')
print('\n')


#function to look up diseases a given miRNA is associated to, whilst reporting the associated gene(s).

def miRNA_related_diseases(mirna:str = None, gene_conf:int = 80, disease_conf:float = 0.8, score:str = 'True') -> str:
    '''
    Func: 
        diseases a given miRNA is associated to, whilst reporting the associated gene(s).
     
     Parameters:
    - mirna (str): The miRNA to query.
    - gene_conf (int): Minimum confidence score for miRNA to Gene associations.
    - disease_conf (float): Minimum confidence score for Gene to Disease associations.
    - score (str): Determines whether to include confidence scores in the output.

    input:
        mirna
    output:
        related diseases and gene(s)
    '''
#select columns based on whether score information should be included
    if score == 'True':
        cols = ['Gene','Score','diseaseName','score']
    elif score == 'False':
        cols = ['Gene','diseaseName']        
#query the dataset for the specified miRNA, gene confidence, and disease confidence        
    mirna_related_diseases = source_data.loc[(source_data['miRNA']==mirna) & (source_data['Score']>gene_conf) & (source_data['score']>disease_conf),cols]
    if mirna_related_diseases.empty:
        print('No relevant information was found.')
        
    else:
        print(mirna_related_diseases.reset_index().drop('index',axis=1).to_markdown())
    print('\n')


#function to look up miRNAs a disease is associated to, whilst reporting the associated gene(s).

def diseases_related_miRNA_Gene(disease:str = None, gene_conf:int = 80, disease_conf:float = 0.8, score:str = 'True') -> str:
    '''
    Func: 
        look up miRNAs a disease is associated to, whilst reporting the associated gene(s).
    
    Parameters:
    - disease (str): The disease to query.
    - gene_conf (int): Minimum confidence score for miRNA to Gene associations.
    - disease_conf (float): Minimum confidence score for Gene to Disease associations.
    - score (str): Determines whether to include confidence scores in the output.

    input:
        disease
    output:
        related mirnas and gene(s)
    '''
    print('diseases related miRNA.\n')
# Select columns based on whether score information should be included
    if score == 'True':
        cols = ['miRNA','Score','Gene','score']
    elif score == 'False':
        cols = ['miRNA','Gene']
#query the dataset for the specified disease, gene confidence, and disease confidence    
    diseases_related_mirna_gene = source_data.loc[(source_data['diseaseName']==disease) & (source_data['Score']>gene_conf) & (source_data['score']>disease_conf),cols]
    if diseases_related_mirna_gene.empty:
        print('No relevant information was found.')
    else:
        print(diseases_related_mirna_gene.reset_index().drop('index',axis=1).to_markdown())
    print('\n')


# browse what miRNA are in the database

def browse_miRNA(term: str = None) -> str:
    '''
    Func: 
        Have an option to browse what miRNA are in the database.
        - E.g. given part of a miRNA name, it should report back all the miRNAs that match it, e.g. searching for “cfa” should list everything that has “cfa” in the miRNA name such as cfa-miR-1185, cfa-miR-544, cfa-miR-8808, etc.
    input:
        miRNA_term
    output:
    - Prints the list of diseases that match the search term. If 'ALL' is used, prints all disease names.
    - Returns a string containing the list of matching disease names, each separated by a newline.
    '''

#check if the term is set to 'ALL' to list all miRNA
    if term == 'ALL':
        print(list(source_data['miRNA'].unique()))
        print('\n')
        return '\n'.join(source_data['miRNA'].unique())
        
    else:
        if source_data.loc[source_data['miRNA'].str.contains(term),'miRNA'].empty:
            print('No relevant information was found.')
        
        else:
            print(list(source_data.loc[source_data['miRNA'].str.contains(term),'miRNA'].unique()))
        print('\n')
        return '\n'.join(source_data.loc[source_data['miRNA'].str.contains(term),'miRNA'].unique())
        

# browse what diseases are in the database

def browse_diseases(term: str = None) -> str:
    '''
    Func:
        Allows browsing of diseases in the database based on a search term.
    input:
        disease_term
    output:
    - Prints the list of diseases that match the search term. If 'ALL' is used, prints all disease names.
    - Returns a string containing the list of matching disease names, each separated by a newline.
    '''
    
#check if the term is set to 'ALL' to list all diseases
    if term == 'ALL':
        print(list(source_data['diseaseName'].unique()))
        print('\n')
        return '\n'.join(source_data['diseaseName'].unique())
    
#search for diseases that contain the given term, if term is not None       
    else:
        if source_data.loc[source_data['diseaseName'].str.contains(term),'diseaseName'].empty:
            print('No relevant information was found.')
        else:
            print(list(source_data.loc[source_data['diseaseName'].str.contains(term),'diseaseName'].unique()))
        print('\n')
        return '\n'.join(source_data.loc[source_data['diseaseName'].str.contains(term),'diseaseName'].unique())


#check if both disease and mirna arguments are provided, which is not allowed

if disease and mirna:
    #print error message if both miRNA and disease arguments are present
    print('miRNA and disease cannot exist at the same time')
    print('\n')

elif disease:
    #if only the disease argument is provided
    if browse:
        #if the browse flag is set, perform a browsing search for diseases
        browse_diseases(term = disease)
    
    else:
        #if the browse flag is not set, find miRNAs related to the provided disease
        diseases_related_miRNA_Gene(disease = disease, gene_conf = gene_conf, disease_conf = disease_conf, score=score)
    
elif mirna:
    #if only the mirna argument is provided
    if browse:
        #if the browse flag is set, perform a browsing search for miRNAs
        browse_miRNA(term = mirna)
    
    else:
        #if the browse flag is not set, find diseases related to the provided miRNA
        miRNA_related_diseases(mirna = mirna, gene_conf = gene_conf, disease_conf = disease_conf, score=score)
    
else:
    #print an error message if neither miRNA nor disease arguments are provided
    print('miRNA and disease cannot exist at the same time')
    print('\n')