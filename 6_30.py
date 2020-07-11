from bravado.client import SwaggerClient
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from lifelines.plotting import plot_lifetimes
from lifelines import KaplanMeierFitter

cbioportal = SwaggerClient.from_url('https://www.cbioportal.org/api/api-docs',
                                config={"validate_requests":False,"validate_responses":False})

def graph(full, mutation):
    months, survival_status, overall_mutations = no_mutation(full, mutation)

    survival_data=pd.DataFrame({'OS_MONTHS': months,
                                'OS_STATUS': survival_status # 0 if living, 1 if dead
                                })

    has_mutation=np.array(overall_mutations) == 1 #0 if don't have mutation, 1 if do have mutation

    ## create an kmf object
    kmf=KaplanMeierFitter()

    ## fit the data into a model for each group
    kmf.fit(survival_data.OS_MONTHS[has_mutation], survival_data.OS_STATUS[has_mutation], label="have mutation")
    layer1=kmf.plot(ci_show=True)

    kmf.fit(survival_data.OS_MONTHS[~has_mutation], survival_data.OS_STATUS[~has_mutation], label="no mutation")
    layer2=kmf.plot(ax=layer1, ci_show=True)

    ## view plot
    plt.show()
    
def no_mutation (n,m):
        #list containing patient IDs with mutation of gene
        no=[]

        #list containing patient IDs without mutation of gene
        use=[]

        overall_mutations=[]

        everything=[]
        #x is declared as 0 and will be incremented inside the while loop until x is the number of patient IDs in the original list with all patient IDs
        x=0
        while x<len(n):
                for i in m:
                        #TCGA-OL-A66H does not have the attribute OS_MONTHS
                        #TCGA-BH-A0B2 does not have OS_MONTHS, only has AGE, AJCC_PATHOLOGIC_TUMOR_STAGE, AJCC_STAGING_EDITION, CANCER_TYPE_ACRONYM, CENTER, DAYS_LAST_FOLLOWUP, DAYS_TO_BIRTH, DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS, ETHNICITY, FORM_COMPLETION_DATE, HISTORY_NEOADJUVANT_TRTYN, ICD_10, ICD_O_3_HISTOLOGY, ICD_O_3_SITE, INFORMED_CONSENT_VERIFIED, "IN_PANCANPATHWAYS_FREEZE, OTHER_PATIENT_ID, PATH_M_STAGE, PATH_N_STAGE, PATH_T_STAGE, PERSON_NEOPLASM_CANCER_STATUS, PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT, PRIOR_DX, RACE, SAMPLE_COUNT, SEX
                        if n[x]==i or n[x]=='TCGA-BH-A0B2' or n[x]=='TCGA-OL-A66H':
                                #if the patient ID is in the list of those with a mutation then the patient ID will be added to the list that will not be used
                                no.append(i)
                                everything.append(i)
                                overall_mutations.append(1)                                
                                x=x+1

                #if the patient ID is not in the list of those with a mutation then the patient ID will be added to the list that will be outputted by this method
                use.append(n[x])
                everything.append(n[x])
                overall_mutations.append(0)

                #increments x
                x=x+1

        print("all the patient IDs used in the study {} ".format(everything))
        print("overall mutation list {} ".format(overall_mutations))
        print(" ")
        survival_status=[]
        months_living=[]
        for j in everything:
                #outputs how many months that person has been alive after diagnosis
                months_living.append(cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_MONTHS', patientId=j, studyId='brca_tcga_pan_can_atlas_2018').result()[0])
    
                #outputs current age of person without mutation
                living = cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_STATUS', patientId=j, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value']
                if living == '0:LIVING':
                    survival_status.append(0)
                elif living == '1:DECEASED':
                    survival_status.append(1)

        months = [x.value for x in months_living]
        return months, survival_status, overall_mutations

def main():

	# some examples of information we can query
	studies = cbioportal.Studies.getStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
	cancer_types = cbioportal.Cancer_Types.getAllCancerTypesUsingGET().result()
	

	# extended documentation available here https://www.cbioportal.org/api/swagger-ui.html

	# select patients in the cohort of interest (TCGA pan cancer project)
	patients = cbioportal.Patients.getAllPatientsInStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
	print("The brca_tcga_pan_can_atlas_2018 study spans {} patients".format(len(patients)))

	# select genes in the cohort of interest
	genes = cbioportal.Genes.getGeneUsingGET(geneId='MAP3K1').result() #TP53 = 7157, EP300 = 2033, PIK3CA=5290, CDH1=999, GATA3=2625, MAP3K1=4214
	print("The Entrez Gene ID for gene MAP3K1 is {} ".format(genes.entrezGeneId))

	
	# what kind of mutations do the patients in this cohort have? 
	mutations = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
    	entrezGeneId=4214,
        molecularProfileId='brca_tcga_pan_can_atlas_2018_mutations',
    	sampleListId='brca_tcga_pan_can_atlas_2018_all',
    	projection='DETAILED'
	).result()
	print("The number of patients with a mutation of the EP300 gene is {} ".format(len(mutations)))
	patients = cbioportal.Patients.getAllPatientsInStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
	patientIds = [x.patientId for x in patients]

	mutation_EP300 = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(entrezGeneId=4214, molecularProfileId='brca_tcga_pan_can_atlas_2018_mutations', sampleListId='brca_tcga_pan_can_atlas_2018_all').result()
	patient_EP300=[x.patientId for x in mutation_EP300]
	
	graph(patientIds,patient_EP300)

if __name__ == '__main__':
	main()
