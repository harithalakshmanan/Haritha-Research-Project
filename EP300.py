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

	# which genes have mutations? create plot to visualize
	mutated_genes = Counter([m.gene.hugoGeneSymbol for m in mutations])

	print("The brca_tcga_pan_can_atlas_2018 study spans {} mutations, in {} genes".format(
		  len(mutations),
		  len(mutated_genes)
	))

	plot_mutations(mutated_genes)

if __name__ == '__main__':
	main()
