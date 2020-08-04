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

def bernoulliNaiveBayes(has_mutation, survival_status):
    rng = np.random.RandomState(1)
    a = has_mutation.reshape(len(has_mutation),1)
    b=survival_status
    clf=BernoulliNB()
    clf.fit(a,b)
    print(clf.predict(a))
          
def statisticalSignificance(survival_status, has_mutation):
    i=0
    aliveNo=0     #if patient is alive and does not have mutation
    aliveYes=0    #if patient is alive and has mutation
    deceasedNo=0  #if patient died and does not have mutation
    deceasedYes=0 #if patient died and has mutation
    print(len(has_mutation))
    while i<len(survival_status):
        if has_mutation[i] == 1:
            if survival_status[i] == 0:
                aliveYes = aliveYes + 1
            else:
                deceasedYes = deceasedYes + 1
        else:
            if survival_status[i] == 0:
                aliveNo = aliveNo + 1
            else:
                deceasedNo = deceasedNo + 1
        i = i + 1
        #increments i
    print("aliveNo {} ".format(aliveNo))
    print("aliveYes {} ".format(aliveYes))
    print("deceasedNo {} ".format(deceasedNo))
    print("deceasedYes {} ".format(deceasedYes))
    oddsratio, pvalue = stats.fisher_exact([[deceasedYes, aliveYes], [deceasedNo, aliveNo]])
    print("pvalue {} ".format(pvalue))
    
def graph(months, survival_status, has_mutation, name):

    survival_data=pd.DataFrame({'OS_MONTHS': months,
                                'OS_STATUS': survival_status # 0 if living, 1 if dead
                                })
    #0 if don't have mutation, 1 if do have mutation in has_mutation

    ## create an kmf object
    kmf=KaplanMeierFitter()

    ## fit the data into a model for each group
    kmf.fit(survival_data.OS_MONTHS[has_mutation], survival_data.OS_STATUS[has_mutation], label="have mutation")
    layer1=kmf.plot(ci_show=True)

    kmf.fit(survival_data.OS_MONTHS[~has_mutation], survival_data.OS_STATUS[~has_mutation], label="no mutation")
    layer2=kmf.plot(ax=layer1, ci_show=True)

    plt.title('{} survival plot'.format(name))

    ## view plot
    plt.show()

def getSurvivalData(patientIds, mutatedIds):
    overall_mutations=np.isin(patientIds, mutatedIds)
    
    months=[cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_MONTHS', patientId=j, studyId='brca_tcga_pan_can_atlas_2018').result()[0] for j in patientIds]
    months=[float(x.value) for x in months]
    
    living=[cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_STATUS', patientId=i, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value'] for i in patientIds]
    survival_status=np.array(living)=='1:DECEASED'
    
    return months, survival_status, overall_mutations  
    

def anomolies(patientIds, mutatedIds):
    #TCGA-OL-A66H does not have the attribute OS_MONTHS
    #TCGA-BH-A0B2 does not have OS_MONTHS, only has AGE, AJCC_PATHOLOGIC_TUMOR_STAGE, AJCC_STAGING_EDITION, CANCER_TYPE_ACRONYM, CENTER, DAYS_LAST_FOLLOWUP, DAYS_TO_BIRTH, DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS, ETHNICITY, FORM_COMPLETION_DATE, HISTORY_NEOADJUVANT_TRTYN, ICD_10, ICD_O_3_HISTOLOGY, ICD_O_3_SITE, INFORMED_CONSENT_VERIFIED, "IN_PANCANPATHWAYS_FREEZE, OTHER_PATIENT_ID, PATH_M_STAGE, PATH_N_STAGE, PATH_T_STAGE, PERSON_NEOPLASM_CANCER_STATUS, PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT, PRIOR_DX, RACE, SAMPLE_COUNT, SEX
    anomolies=(['TCGA-BH-A0B2', 'TCGA-OL-A66H'])
    
    mask=np.isin(patientIds, anomolies)
    patientIds=patientIds[~mask]
    
    mutate=np.isin(mutatedIds, anomolies)
    mutatedIds=mutatedIds[~mutate]
    mutatedIds=np.unique(mutatedIds) #often times, the same patient ID has multiple mutations of the gene, so unique prevents the same patient ID from being in the list of mutated Ids
    return(patientIds, mutatedIds)
    
def genes(name):
    # select genes in the cohort of interest
    #TP53 = 7157, EP300 = 2033, PIK3CA=5290, CDH1=999, GATA3=2625, MAP3K1=4214
    genes = cbioportal.Genes.getGeneUsingGET(geneId=name).result()
    print("The Entrez Gene ID for gene {} is {} ".format(name, genes.entrezGeneId))
    return genes.entrezGeneId
	
def main():
    name='GATA3'
    geneId=genes(name)

    # extended documentation available here https://www.cbioportal.org/api/swagger-ui.html
    # select patients in the cohort of interest (TCGA pan cancer project)
    patients = cbioportal.Patients.getAllPatientsInStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
    patientIds = np.array([x.patientId for x in patients])
    print("The brca_tcga_pan_can_atlas_2018 study spans {} patients".format(len(patients)))

    # what kind of mutations do the patients in this cohort have?
    mutation = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(entrezGeneId=geneId, molecularProfileId='brca_tcga_pan_can_atlas_2018_mutations', sampleListId='brca_tcga_pan_can_atlas_2018_all').result()
    mutatedIds=np.array([x.patientId for x in mutation])
    print("The number of mutations of the {} gene is {} ".format(name, len(mutation))) #this outputs the total number of mutations of a particular gene which does not need to be the total number of people with the mutation, as the same person could have multiple mutations for the same gene

    patient, mutated = anomolies(patientIds, mutatedIds)
    print("Patients used to graph {} ".format(len(patient)))
    print("total mutated {} ".format(len(mutated)))
    months, survival_status, overall_mutations = getSurvivalData(patient, mutated)
    graph(months, survival_status, overall_mutations, name)
    statisticalSignificance(survival_status, overall_mutations)
    bernoulliNaiveBayes(overall_mutations, survival_status)

if __name__ == '__main__':
	main()
