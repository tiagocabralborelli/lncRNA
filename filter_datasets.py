import json
import requests
import pandas as pd
from Bio import Entrez
import re
import ftplib
import ppx
import datetime

'''
from bs4 import BeautifulSoup
soup = BeautifulSoup(requests.get('https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10090').content, 'html.parser')
info = soup.find_all("dd")
for dl in info:
    print(dl.string)
'''

#Get a list of all datasets
def getAllDts(base_url = "massive.ucsd.edu/"):
    """
    Search on all datasets from massIVE
    """
    datasets_url = 'https://' + base_url + '/ProteoSAFe/datasets_json.jsp'
    json_obj = json.loads(requests.get(datasets_url).text)
    return json_obj["datasets"]

def isPlant(taxonid, taxa, filt=True):
    handle = Entrez.efetch(db="Taxonomy", id=taxonid, retmode="xml")
    records = Entrez.read(handle)
    if filt:
        return bool(re.search(taxa, records[0]['Lineage']))
    else:
        return records

def filterMassIVE(df, ntoexclude=5, taxa='plantae'):
    exclude = df.species.value_counts().index[:ntoexclude].tolist()
    search_list = df.species[~df.species.isin(exclude)].value_counts().index.tolist()
    plants = []
    missing = []
    for s in search_list:
        id_list = [x for x in s.split(';') if x not in exclude]
        if len(id_list):
            for idx in id_list:
                taxonid = re.sub('.+NCBITaxon:(\d+).+', '\\1', idx)
                print("taxonid: %s" % taxonid)
                try:
                    if isPlant(re.sub('\D', '', idx), taxa):
                        print('Found a plant: %s' % idx)
                        plants.append(df[df.species==s])
                    else:
                        print('Not a plant: %s' % idx)
                except:
                    print('cant process id: %s' % idx)
                    missing.append(idx)
                    continue
    pdf = pd.concat(plants)
    return pdf[~pdf.dataset.duplicated()].reset_index(drop=True)

def findSearchMassIVE(datasets):
    url="massive.ucsd.edu"
    search = []
    for dr in datasets:
        ftp = ftplib.FTP(url)
        ftp.login("anonymous", "password")
        ftp.cwd(dr)
        if 'search' in ftp.nlst():
            ftp.cwd('search')
            ftp.dir()
            search.append(1)
        else:
            search.append(0)
        ftp.close()
    return search

def downloadPublicFTP(url="massive.ucsd.edu", dr="MSV000088167"):
    ftp = ftplib.FTP(url)
    ftp.login("anonymous", "password")
    ftp.cwd(dr)
    if 'search' in ftp.nlst():
        ftp.cwd('search')
    ftp.cwd('Prep1_searched')
    fls = ftp.nlst()
    [ftp.retrbinary("RETR " + x, open(x, 'wb').write) for x in fls]
    ftp.close()

#https://ppx.readthedocs.io/en/latest/api/functions.html
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489246/
# https://github.com/PRIDE-Archive?language=python

def filterPride(taxa='plantae'):
    plant_projects = {}
    projects = ppx.pride.list_projects()
    for p in projects:
        project = ppx.PrideProject(p)
        acession = [(x['accession'], x['name']) for x in project.metadata['organisms']]
        for a in acession:
            try:
                if isPlant(a, taxa):
                    print('Found plant! Acession: %s' % a[1])
                    plant_projects[p] = {}
                    plant_projects[p]['project'] = project
                    plant_projects[p]['accession'] = a[0]
                    plant_projects[p]['name'] = a[1]
                    break
                else:
                    print('Not a plant! Acession: %s' % a[1])
            except:
                    print('cant process acession: %s' % a[1])

    plant_projects = {}
    for k,v in plant_projects.items():
        tmp = {}
        tmp['name'] = v['name']
        tmp['accession'] = v['accession']
        plant_projects[k] = tmp

    with open('pride_plant.json', 'w+') as f:
        json.dump(plant_projects, f, indent=4)

def findmztab(identifier):
    proj = ppx.find_project(identifier)
    remote_files = proj.remote_files()
    mztab = [x for x in remote_files if 'mztab' in x]
    if len(mztab):
        return {identifier: mztab}

def get_mztabs(pride_dataset:dict) -> dict:
    #Função para baixar e abrir todos os arquivos mztabe de um projeto
    project_dict_map = {}
    for projects,values in pride_dataset.items():
        try:
            mztabs_from_projects = findmztab(projects)
            if mztabs_from_projects != None:
                project_dict_map[projects] = mztabs_from_projects[projects]
        except:
            print(f"ERROR:{datetime.datetime.now().strftime('%d/%m/%Y %H:%M')} - {projects} \n")
    return project_dict_map
def pride_summary():
    #Organiza os dados do pride
    pass
