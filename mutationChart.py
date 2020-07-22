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
	print("The number of patients with a mutation of the 7157 gene is {} ".format(len(mutations)))
	
	# Id stores the patientId
	Id='TCGA-A2-A0YF'

	# how many months has this patient been alive (overall survival in months)?
	living = cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_MONTHS', patientId=Id, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value']
	print("months living {} ".format(living))

	# how old is the patient?3
	status = cbioportal.Clinical_Data.getAllClinicalDataOfPatientInStudyUsingGET(attributeId='OS_STATUS', patientId=Id, studyId='brca_tcga_pan_can_atlas_2018').result()[0]['value']
	print("living {} ".format(status))

	patients = cbioportal.Patients.getAllPatientsInStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
	patientIds = [x.patientId for x in patients]	

	# which genes have mutations? create plot to visualize
	mutated_genes = Counter([m.gene.hugoGeneSymbol for m in mutations])

	print("The brca_tcga_pan_can_atlas_2018 study spans {} mutations, in {} genes".format(
		  len(mutations),
		  len(mutated_genes)
	))

	mutation = cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(entrezGeneId=4214, molecularProfileId='brca_tcga_pan_can_atlas_2018_mutations', sampleListId='brca_tcga_pan_can_atlas_2018_all').result()
	patient=[x.patientId for x in mutation]
	print(patient)
	plot_mutations(mutated_genes)

if __name__ == '__main__':
	main()
