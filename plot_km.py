# full tutorial at https://github.com/anurag-code/Survival-Analysis-Intuition-Implementation-in-Python/blob/master/Survival%20Analysis%20-%20Quick%20Implementation.ipynb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from lifelines.plotting import plot_lifetimes   
from lifelines import KaplanMeierFitter

## Example Data 
survival_data = pd.DataFrame({'OS_MONTHS': [50,60,60,25,10,22],
                              'OS_STATUS': [1,0,0,1,1,1], # 0 if living, 1 if died
                             })

has_mutation = np.array([1,0,0,0,1,1]) == 1 # 0 if don't have mutaiton, 1 if do

## create an kmf object
kmf = KaplanMeierFitter() 

## Fit the data into the model for each group
kmf.fit(survival_data.OS_MONTHS[has_mutation], survival_data.OS_STATUS[has_mutation], label="have mutation")
layer1 = kmf.plot(ci_show=True)

kmf.fit(survival_data.OS_MONTHS[~has_mutation], survival_data.OS_STATUS[~has_mutation], label="don't have mutation")
layer2 = kmf.plot(ax=layer1, ci_show=True)

## view plot
plt.show()