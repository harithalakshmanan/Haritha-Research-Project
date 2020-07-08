from bravado.client import SwaggerClient
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

cbioportal = SwaggerClient.from_url('https://www.cbioportal.org/api/api-docs',
                                config={"validate_requests":False,"validate_responses":False})

def plot_mutations(m):

	# TODO: make a barplot of mutation counts

	mdf = pd.DataFrame.from_dict(m, orient='index')
	mdf = mdf.sort_values(by=[0], ascending = False)
	mdf.head(10).plot(
    	kind='bar',
    	color='green',
    	legend=None
	)
	plt.xlabel('')
	plt.xticks(rotation=300)
	plt.ylabel('Number of samples', labelpad=10)
	plt.title('Number of mutations in 10 genes in\n Breast Invasive Carcinoma (TCGA, PanCancer Atlas)',pad=10)
	plt.subplots_adjust(bottom=0.15)
	plt.show()
	plt.savefig("mutations_top10.png")

def no_mutation (n,m):
        #list containing patient IDs with mutation of gene
        no=[]

        #list containing patient IDs without mutation of gene
        use=[]

        x=0
        while x<len(n):
                for i in m:
                        #TCGA-BH-A0B2 does not have OS_MONTHS, only has AGE, AJCC_PATHOLOGIC_TUMOR_STAGE, AJCC_STAGING_EDITION, CANCER_TYPE_ACRONYM, CENTER, DAYS_LAST_FOLLOWUP, DAYS_TO_BIRTH, DAYS_TO_INITIAL_PATHOLOGIC_DIAGNOSIS, ETHNICITY, FORM_COMPLETION_DATE, HISTORY_NEOADJUVANT_TRTYN, ICD_10, ICD_O_3_HISTOLOGY, ICD_O_3_SITE, INFORMED_CONSENT_VERIFIED, "IN_PANCANPATHWAYS_FREEZE, OTHER_PATIENT_ID, PATH_M_STAGE, PATH_N_STAGE, PATH_T_STAGE, PERSON_NEOPLASM_CANCER_STATUS, PRIMARY_LYMPH_NODE_PRESENTATION_ASSESSMENT, PRIOR_DX, RACE, SAMPLE_COUNT, SEX
                        if n[x]==i or n[x]=='TCGA-BH-A0B2' or n[x]=='TCGA-OL-A66H':
                                no.append(i)
                                x=x+1
                use.append(n[x])
                x=x+1
        print("List of those without mutation {} ".format(use))
        print(" ")
        print("This list will output the information regarding those without a mutation for the selected gene")
        for j in use:
                #outputs patient Id of someone without mutation in gene
                print("Patient Id: {} ".format(j))

                #outputs how many months that person has been alive after diagnosis
                living = cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_MONTHS', patientId=j, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value']
                print("Months Living: {} ".format(living))

                #outputs current age of person without mutation
                age = cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='AGE', patientId=j, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value']
                print("Age: {} ".format(age))

                print(" ")
	
def mutation (z):
        # method in charge of outputting data concerning those with a mutation in gene EP300
        print(" ")
        print("This list will output the information regarding those with a mutation for the selected gene")
        for y in z:
                #outputs patient Id of someone with a mutation in gene EP300
                print("Patient Id: {} ".format(y))

                #outputs how many months that person has been alive after diagnosis
                living = cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_MONTHS', patientId=y, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value']
                print("Months Living: {} ".format(living))

                #outputs current age of person with mutation
                age = cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='AGE', patientId=y, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value']
                print("Age: {} ".format(age))

                print(" ")
def main():

	# some examples of information we can query
	studies = cbioportal.Studies.getStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
	cancer_types = cbioportal.Cancer_Types.getAllCancerTypesUsingGET().result()
	

	# extended documentation available here https://www.cbioportal.org/api/swagger-ui.html

	# select patients in the cohort of interest (TCGA pan cancer project)
	patients = cbioportal.Patients.getAllPatientsInStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
	print("The brca_tcga_pan_can_atlas_2018 study spans {} patients".format(len(patients)))

	# select genes in the cohort of interest
	genes = cbioportal.Genes.getGeneUsingGET(geneId='EP300').result()
	print("The Entrez Gene ID for gene EP300 is {} ".format(genes.entrezGeneId))

	
	# what kind of mutations do the patients in this cohort have? 
	mutations = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
    	entrezGeneId=2033,
        molecularProfileId='brca_tcga_pan_can_atlas_2018_mutations',
    	sampleListId='brca_tcga_pan_can_atlas_2018_all',
    	projection='DETAILED'
	).result()
	print("The number of patients with a mutation of the EP300 gene is {} ".format(len(mutations)))
	patients = cbioportal.Patients.getAllPatientsInStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
	patientIds = [x.patientId for x in patients]
	
                
	# Id stores the patientId
	#Id=['TCGA-A2-A0YF']
	print("patient ID {} ".format(patientIds))

	mutation_EP300 = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(entrezGeneId=2033, molecularProfileId='brca_tcga_pan_can_atlas_2018_mutations', sampleListId='brca_tcga_pan_can_atlas_2018_all').result()
	patient_EP300=[x.patientId for x in mutation_EP300]
	print("patient Ids with a mutation in gene EP300 {} ".format(patient_EP300))

	# which genes have mutations? create plot to visualize
	mutated_genes = Counter([m.gene.hugoGeneSymbol for m in mutations])

	print("The brca_tcga_pan_can_atlas_2018 study spans {} mutations, in {} genes".format(
		  len(mutations),
		  len(mutated_genes)
	))

	mutation(patient_EP300)
	no_mutation(patientIds,patient_EP300)
	plot_mutations(mutated_genes)

if __name__ == '__main__':
	main()
