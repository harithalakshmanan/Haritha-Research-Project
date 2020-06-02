# Example script for querying TCGA data using the cBio portal API
# Adapted from https://github.com/mskcc/cbsp-hackathon/blob/master/0-introduction/cbsp_hackathon.ipynb
# Version: 1.0

# import required packages
from bravado.client import SwaggerClient
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter


# set the database that is to be queried
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
	plt.savefig("images/mutations_top10.png")


def main():

	# some examples of information we can query
	studies = cbioportal.B_Studies.getAllStudiesUsingGET().result()
	cancer_types = cbioportal.A_Cancer_Types.getAllCancerTypesUsingGET().result()
	
	print("In total there are {} studies in cBioPortal, spanning {} different types of cancer.".format(
	      len(studies),
	      len(cancer_types)
	))

	# extended documentaiton available here https://www.cbioportal.org/api/swagger-ui.html

	# select patients in the cohort of interest (TCGA pan cancer project)
	patients = cbioportal.C_Patients.getAllPatientsInStudyUsingGET(studyId='brca_tcga_pan_can_atlas_2018').result()
	print("The brca_tcga_pan_can_atlas_2018 study spans {} patients".format(len(patients)))

	# what kind of mutations do the patients in this cohort have? 
	mutations = cbioportal.K_Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
    	molecularProfileId='brca_tcga_pan_can_atlas_2018_mutations',
    	sampleListId='brca_tcga_pan_can_atlas_2018_all',
    	projection='DETAILED'
	).result()

	# which genes have mutations? create plot to visualize
	mutated_genes = Counter([m.gene.hugoGeneSymbol for m in mutations])

	print("The brca_tcga_pan_can_atlas_2018 study spans {} mutations, in {} genes".format(
		  len(mutations),
		  len(mutated_genes)
	))

	plot_mutations(mutated_genes)

if __name__ == '__main__':
	main()