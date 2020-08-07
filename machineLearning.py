from bravado.client import SwaggerClient
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from collections import Counter
from lifelines.plotting import plot_lifetimes
from lifelines import KaplanMeierFitter
from sklearn.naive_bayes import BernoulliNB

cbioportal = SwaggerClient.from_url('https://www.cbioportal.org/api/api-docs',
                                config={"validate_requests":False,"validate_responses":False})

def bernoulliNaiveBayes(data, survival_status):
    a=data # gets for one mutation
    b=survival_status
    clf=BernoulliNB()
    clf.fit(a,b)
    print(clf.predict(a))
    return (clf.predict(a))

def classificationAccuracy(survival_status, naive_bayes):
    x=0
    i=0
    while i<len(survival_status):
        if survival_status[i]==naive_bayes[i]:
            x+=1
        i+=1
    percent=(i/len(survival_status))*100
    print("The Bernoulli Naive Bayes Classification Algorithm was {}% accurate ".format(percent))

def getSurvivalData(patientIds):    
    living=[cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_STATUS', patientId=i, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value'] for i in patientIds]
    survival_status=np.array(living)=='1:DECEASED'
    
    return survival_status  
    
def overallMutations(Ids):
    data=pd.DataFrame({'GATA3': mutation(genes('GATA3'),Ids),
                       'TP53': mutation(genes('TP53'),Ids),
                       'EP300': mutation(genes('EP300'),Ids),
                       'PIK3CA': mutation(genes('PIK3CA'),Ids),
                       'CDH1': mutation(genes('CDH1'),Ids),
                       })
    print(data)
    return data

def mutation(entrezId,patient):
    mutatedIds=[]
    #create a matrix where every row is full patient (1082 = size)
    mutation = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(entrezGeneId=entrezId, molecularProfileId='brca_tcga_pan_can_atlas_2018_mutations', sampleListId='brca_tcga_pan_can_atlas_2018_all').result()
    mutatedIds.append(np.array([x.patientId for x in mutation]))
    mutatedIds=np.unique(mutatedIds)
    overall_mutations=np.isin(patient,mutatedIds)
    #often times, the same patient ID has multiple mutations of the gene, so unique prevents the same patient ID from being in the list of mutated Ids
    return overall_mutations
    
def anomolies(Ids):
    #this method returns a lists of patient Ids that are not repeated and that do not belong to the anomolies list
    #TCGA-OL-A66H does not have the attribute OS_MONTHS
    #TCGA-BH-A0B2 does not have OS_MONTHS, only has AGE, AJCC_PATHOLOGIC_TUMOR_STAGE, AJCC_STAGING_EDITION, CANCER_TYPE_ACRONYM, CENTER, DAYS_LAST_FOLLOWUP, DAYS_TO_BIRTH, DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS, ETHNICITY, FORM_COMPLETION_DATE, HISTORY_NEOADJUVANT_TRTYN, ICD_10, ICD_O_3_HISTOLOGY, ICD_O_3_SITE, INFORMED_CONSENT_VERIFIED, "IN_PANCANPATHWAYS_FREEZE, OTHER_PATIENT_ID, PATH_M_STAGE, PATH_N_STAGE, PATH_T_STAGE, PERSON_NEOPLASM_CANCER_STATUS, PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT, PRIOR_DX, RACE, SAMPLE_COUNT, SEX
    anomolies=(['TCGA-BH-A0B2', 'TCGA-OL-A66H'])
    mask=np.isin(Ids, anomolies)
    Ids=Ids[~mask] #error message here #TypeError: only integer scalar arrays can be converted to a scalar index
    Ids=np.unique(Ids)  
    return Ids
    
def genes(name):
    # select genes in the cohort of interest
    #TP53 = 7157, EP300 = 2033, PIK3CA=5290, CDH1=999, GATA3=2625, MAP3K1=4214
    genes = cbioportal.Genes.getGeneUsingGET(geneId=name).result()
    return genes.entrezGeneId
	
def main():
    # extended documentation available here https://www.cbioportal.org/api/swagger-ui.html
    # select patients in the cohort of interest (TCGA pan cancer project)
    patients = cbioportal.Patients.getAllPatientsInStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
    patientIds = np.array([x.patientId for x in patients])
    print("The brca_tcga_pan_can_atlas_2018 study spans {} patients".format(len(patients)))

    patient = anomolies(patientIds)
    print("Patients used to graph {} ".format(len(patient)))

    overall_mutations = overallMutations(patient)
    survival_status = getSurvivalData(patient)
    naive_bayes = bernoulliNaiveBayes(overall_mutations, survival_status)
    classificationAccuracy(survival_status, naive_bayes)

if __name__ == '__main__':
	main()
