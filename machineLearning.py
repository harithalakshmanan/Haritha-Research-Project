from bravado.client import SwaggerClient
from collections import Counter
from lifelines.plotting import plot_lifetimes
from lifelines import KaplanMeierFitter
from sklearn.naive_bayes import BernoulliNB
from collections import defaultdict
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import RFECV
from sklearn.svm import SVR

import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import mygene

cbioportal = SwaggerClient.from_url('https://www.cbioportal.org/api/api-docs',
                                config={"validate_requests":False,"validate_responses":False})

def featureSelection(matrix, survival):
    #train test split
    matrix_train, matrix_test, survival_train, survival_test = train_test_split(matrix, survival)
    clf = BernoulliNB()
    clf.fit(matrix_train, survival_train)
    print("bernoulli classification accuracy")
    classificationAccuracy(matrix_test, survival_test)
    
    estimator = BernoulliNB()
    selector = RFECV(estimator, step=50, verbose=1)
    selector = selector.fit(matrix_train, survival_train)
    
    print(selector.ranking_)
    print(selector.predict(matrix_train))
    print(selector.predict(matrix_test))
    print("train data classification accuracy")
    classificationAccuracy(selector.predict(matrix_train), survival_train)
    print("test data classification accuracy")
    classificationAccuracy(selector.predict(matrix_test), survival_test)
    
def classificationAccuracy(data, survival):
    data = np.array(data)
    survival = np.array(survival)
    data = data.flatten()
    x = 0
    i = 0
    while i<len(survival):
        if survival[i] == data[i]:
            x += 1
        i += 1
    percent = (x/len(survival))*100
    print("The Classification Algorithm was {}% accurate ".format(percent))
    print()

def dataSubset(survival_status, patient_matrix):
    #astype(bool) converts 1 to True and 0 to False
    #in survival_status, 1 represents deceased
    #a is a list of the indices where survival_status is 0
    a = list(np.where(~survival_status.astype(bool))[0])
    #survival is a random list of 150 indices
    random.seed(a=20)
    survival = random.sample(a,150)

    #deceased is a list of the indices where survival_status is 1
    deceased = list(np.where(survival_status.astype(bool))[0])

    #patients concatenates or joins the arrays together
    patients = survival + deceased
    patients = np.array(patients)

    patient_matrix = patient_matrix.to_numpy()
    living_matrix = patient_matrix[patients] #np.ndarray
    status = survival_status[patients] #np.ndarray
    return living_matrix, status
    
def getSurvivalData(patientIds):    
    living = [cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_STATUS', patientId=i, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value'] for i in patientIds]
    survival_status = np.array(living) == '1:DECEASED'
    print(survival_status)
    return survival_status  
    
def anomolies(Ids):
    #this method returns a lists of patient Ids that are not repeated and that do not belong to the anomolies list
    #TCGA-OL-A66H does not have the attribute OS_MONTHS
    #TCGA-BH-A0B2 does not have OS_MONTHS, only has AGE, AJCC_PATHOLOGIC_TUMOR_STAGE, AJCC_STAGING_EDITION, CANCER_TYPE_ACRONYM, CENTER, DAYS_LAST_FOLLOWUP, DAYS_TO_BIRTH, DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS, ETHNICITY, FORM_COMPLETION_DATE, HISTORY_NEOADJUVANT_TRTYN, ICD_10, ICD_O_3_HISTOLOGY, ICD_O_3_SITE, INFORMED_CONSENT_VERIFIED, "IN_PANCANPATHWAYS_FREEZE, OTHER_PATIENT_ID, PATH_M_STAGE, PATH_N_STAGE, PATH_T_STAGE, PERSON_NEOPLASM_CANCER_STATUS, PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT, PRIOR_DX, RACE, SAMPLE_COUNT, SEX
    anomolies = (['TCGA-BH-A0B2', 'TCGA-OL-A66H'])
    mask = np.isin(Ids, anomolies)
    Ids = Ids[~mask] #error message here #TypeError: only integer scalar arrays can be converted to a scalar index
    Ids = np.unique(Ids)  
    return Ids

def get_sample_matrix(studyId):
    # using defaultdict as suggested on stackoverflow https://stackoverflow.com/a/41165807
    # this is a modified version of a dictionary, which means we don't have to check if a key is already in the dict when we add a new value
    print("querying mutations...")
    all_mutations = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(molecularProfileId= studyId + '_mutations', sampleListId= studyId + '_all').result()

    # exclude invalid entrez ids; negative id is invalid, 101243544 has been removed from database
    all_mutations = [m for m in all_mutations if m.entrezGeneId not in [-66, 101243544]]
    
    updated_entrez_lookup = {9142: 84631, 23285: 284697, 26148: 84458, 83935: 143872, 114299: 445815, 117153: 4253, \
                             284083: 5414, 348738: 6241, 401388: 7979, 645840: 114112, 100127889: 387707}
    
    print("building patient matrix...")
    # use defaultdict to keep track of which patients have which genes mutated
    # edit updated entrez ids where necessary
    genes_by_patient = defaultdict(list)
    for mutation in all_mutations:
        if mutation.entrezGeneId in updated_entrez_lookup:
            mutation.entrezGeneId = updated_entrez_lookup[mutation.entrezGeneId]
        genes_by_patient[mutation.patientId].append(mutation.entrezGeneId)

    # not all patients are present in this dict (1066 total) because some don't have mutations in the selected mollecular profile
    # remove the pateients with missing clinical data
    genes_by_patient.pop('TCGA-BH-A0B2')
    genes_by_patient.pop('TCGA-OL-A66H')
    

    # remove duplicate genes (patient that have >1 mutation in these genes)
    genes_by_patient = {p:np.unique(genes_by_patient[p]) for p in genes_by_patient.keys()}

    # use multilabelBinarizer to binarize this dict, as suggested on stackoverflow https://stackoverflow.com/a/47209945
    s = pd.Series(genes_by_patient)
    mlb = MultiLabelBinarizer()
    d = mlb.fit_transform(s)

    patient_matrix = pd.DataFrame(d, s.index, mlb.classes_)
    
    # use mygene to translate entrez ids to gene symbol https://pypi.org/project/mygene/
    mg = mygene.MyGeneInfo()
    # use querymany to make a lookup dictionary of entrez ids to gene symbol
    # similar to entrez id's with updated versions.
    print('querying gene names...')
    gene_names_lookup = mg.querymany(patient_matrix.columns, fields = 'symbol')
    gene_names_lookup = {int(g['query']): g['symbol'] for g in gene_names_lookup}
    #patient_matrix.columns = patient_matrix.columns.astype(str)
    patient_matrix = patient_matrix.rename(gene_names_lookup, axis = 'columns')

    return s.index, patient_matrix

def main():
    # extended documentation available here https://www.cbioportal.org/api/swagger-ui.html
    # select patients in the cohort of interest (TCGA pan cancer project)
    patients = cbioportal.Patients.getAllPatientsInStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
    patientIds = np.array([x.patientId for x in patients])
    print("The brca_tcga_pan_can_atlas_2018 study spans {} patients".format(len(patients)))

    patient = anomolies(patientIds)
    print("Patients used to graph {} ".format(len(patient)))

    patientIds, patient_matrix = get_sample_matrix(studyId='brca_tcga_pan_can_atlas_2018')
    survival_status = getSurvivalData(patientIds)
    matrix, status = dataSubset(survival_status, patient_matrix)
    featureSelection(matrix, status)

if __name__ == '__main__':
    main()
