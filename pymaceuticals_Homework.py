#!/usr/bin/env python
# coding: utf-8

# ## Observations and Insights 

# 1- There seems to be a strong correlation (.84)between mouse weight and tumor volume that could be a contributing factor to the effectiveness of any drug regimen.
# 2- There appears to be one potential outlier within the Infubinol regimen (see box plot below). While most mice showed tumor volume increase, there was one mouse that had a reduction in tumor growth in the study.
# 3-Capomulin had the most number of mice complete the study (See bar graph below).

# In[225]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
mouse_study_df = pd.merge(study_results, mouse_metadata,  how = "left", on = ["Mouse ID","Mouse ID"]) 

# Display the data table for preview
mouse_study_df.head()


# In[227]:


# Checking the number of mice.
mouse_study_df['Mouse ID'].nunique()


# In[233]:


#Crosscheck
mice=len(mouse_study_df['Mouse ID'].unique())
mice


# In[234]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_mouse=mouse_study_df.loc[mouse_study_df.duplicated(["Mouse ID", "Timepoint"]), "Mouse ID"].unique()
duplicate_mouse


# In[231]:


# Optional: Get all the data for the duplicate mouse ID. 
mouse_study_df=mouse_study_df.drop_duplicates(subset=['Mouse ID','Timepoint','Tumor Volume (mm3)', 'Age_months', 'Weight (g)'])
mouse_study_df[mouse_study_df['Mouse ID']=='g989']


# In[238]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_study_data_complete=mouse_study_df[mouse_study_df['Mouse ID'].isin(duplicate_mouse)==False]
clean_study_data_complete.head()
                            


# In[243]:


# Checking the number of mice in the clean DataFrame.
mice=len(clean_study_data['Mouse ID'].unique())
mice


# ## Summary Statistics

# In[245]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
means=clean_study_data_complete.groupby('Drug Regimen').mean()['Tumor Volume (mm3)']
medians=clean_study_data_complete.groupby('Drug Regimen').median()['Tumor Volume (mm3)']
variances=clean_study_data_complete.groupby('Drug Regimen').var()['Tumor Volume (mm3)']
sds=clean_study_data_complete.groupby('Drug Regimen').std()['Tumor Volume (mm3)']
sem=clean_study_data_complete.groupby('Drug Regimen').sem()['Tumor Volume (mm3)']
summary_table=pd.DataFrame({"Mean Tumor Volume":means,
                           "Median Tumor Volume": medians,
                           "Tumor Volume Variance": variances,
                           "Tumor Volume Standard Deviation":sds,
                           "Tumor Volume Standard Error": sem})
summary_table.head()


# In[292]:


study_table=clean_study_data_complete.groupby(['Drug Regimen']).agg({'Tumor Volume (mm3)':["mean","median","var","std","sem"]})
study_table


# ## Bar and Pie Charts

# In[417]:


# Generate a bar plot showing the total number of mice for each treatment throughout the course of the study using pandas. 
counts = clean_study_data_complete['Drug Regimen'].value_counts()
counts.plot(kind="bar")
plt.title("Number of Mice per Treatment",fontsize = 12)
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Number of Data Points")
plt.show()

plt.savefig("../Images/mat_mice_per_treat.png", bbox_inches = "tight")


# In[615]:


# Generate a bar plot showing the total number of mice for each treatment throughout the course of the study using pyplot.
counts = clean_study_data_complete['Drug Regimen'].value_counts()
x_axis = np.arange(len(count))
fig1, ax1 = plt.subplots(figsize=(7, 5))
plt.bar(x_axis, count, color='blue', alpha=0.75, align='center')

tick_locations = [value for value in x_axis]

plt.xticks(tick_locations, ['Capomulin', 'Ceftamin', 'Infubinol', 'Ketapril', 'Naftisol', 'Placebo', 'Propriva', 'Ramicane', 'Stelasyn', 'Zoniferol'],  rotation='vertical')

plt.xlim(-0.55, len(x_axis)-0.20)

plt.ylim(0, max(mice_list)+10)

plt.title("Number of Mice per Treatment",fontsize = 12)
plt.xlabel("Drug Regimen",fontsize = 12)
plt.ylabel("Number of Mice",fontsize = 12)

plt.savefig("../Images/mat_mice_per_treat.png", bbox_inches = "tight")


# In[407]:


# Generate a pie plot showing the distribution of female versus male mice using pandas
labels = ["Male", "Female"]
sizes = ["51","49"]
colors = ["blue", "orange"]
explode = (0, 0)
plt.ylabel("Sex")
plt.title('Male vs Female Mouse Population',fontsize = 15)
plt.pie(sizes, explode=explode, labels=labels, colors=colors,
        autopct="%1.1f%%", shadow=True, startangle=120)

plt.savefig("../Images/pi_plot.png", bbox_inches = "tight")
plt.show()


# In[409]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot

labels = ["Male", "Female"]
sizes = ["51","49"]
colors = ["blue", "orange"]
explode = (0, 0)
fig1, ax1 = plt.subplots(figsize=(3.5, 5))
plt.pie(sizes, explode=explode,labels=labels, colors=colors, autopct="%1.1f%%", shadow=True, startangle=140,)
plt.title('Male vs Female Mouse Population',fontsize = 15)
plt.ylabel('Sex',fontsize = 14)
#Set equal axis
plt.axis("equal",fontsize = 14)

plt.savefig("../Images/pi_plot.png", bbox_inches = "tight")
plt.show()


# ## Quartiles, Outliers and Boxplots

# In[414]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin
Capomulin_df = clean_study_data_complete.loc[clean_study_data_complete["Drug Regimen"] == "Capomulin",:]
Ramicane_df = clean_study_data_complete.loc[clean_study_data_complete["Drug Regimen"] == "Ramicane", :]
Infubinol_df = clean_study_data_complete.loc[clean_study_data_complete["Drug Regimen"] == "Infubinol", :]
Ceftamin_df = clean_study_data_complete.loc[clean_study_data_complete["Drug Regimen"] == "Ceftamin", :]
Capomulin_last = Capomulin_df.groupby('Mouse ID').max()['Timepoint']
Capomulin_vol = pd.DataFrame(Capomulin_last)
Capomulin_merge = pd.merge(Capomulin_vol, clean_study_data_complete, on=("Mouse ID","Timepoint"),how="left")
Capomulin_merge.head()


# In[429]:


# Put treatments into a list for for loop (and later for plot labels)
Capomulin_tumors = Capomulin_merge["Tumor Volume (mm3)"]

quartiles =Capomulin_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


print(f"The lower quartile of Capomulin tumors: {lowerq}")
print(f"The upper quartile of Capomulin tumors: {upperq}")
print(f"The interquartile range of Capomulin tumors: {iqr}")
print(f"The median of Capomulin tumors: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)

print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[433]:


Ramicane_last = Ramicane_df.groupby('Mouse ID').max()['Timepoint']
Ramicane_vol = pd.DataFrame(Ramicane_last)
Ramicane_merge = pd.merge(Ramicane_vol, clean_study_data_complete, on=("Mouse ID","Timepoint"),how="left")
Ramicane_merge.head()
Ramicane_merge.to_csv("output.csv")


# In[434]:


Ramicane_tumors = Ramicane_merge["Tumor Volume (mm3)"]

quartiles =Ramicane_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


print(f"The lower quartile of Ramicane tumors is: {lowerq}")
print(f"The upper quartile of Ramicane tumors is: {upperq}")
print(f"The interquartile range of Ramicane tumors is: {iqr}")
print(f"The median of Ramicane tumors is: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)

print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[435]:


Infubinol_last = Infubinol_df.groupby('Mouse ID').max()['Timepoint']
Infubinol_vol = pd.DataFrame(Infubinol_last)
Infubinol_merge = pd.merge(Infubinol_vol, clean_study_data_complete, on=("Mouse ID","Timepoint"),how="left")
Infubinol_merge.head()


# In[436]:


Infubinol_tumors = Infubinol_merge["Tumor Volume (mm3)"]

quartiles =Infubinol_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq


print(f"The lower quartile of Infubinol tumors is: {lowerq}")
print(f"The upper quartile of Infubinol tumors is: {upperq}")
print(f"The interquartile range of Infubinol tumors is: {iqr}")
print(f"The median of Infubinol tumors is: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)


print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")
Infubinol_merge.to_csv("output.csv")


# In[439]:


Ceftamin_last = Ceftamin_df.groupby('Mouse ID').max()['Timepoint']
Ceftamin_vol = pd.DataFrame(Ceftamin_last)
Ceftamin_merge = pd.merge(Ceftamin_vol, clean_study_data_complete, on=("Mouse ID","Timepoint"),how="left")
Ceftamin_merge.head()


# In[440]:


Ceftamin_tumors = Ceftamin_merge["Tumor Volume (mm3)"]

quartiles = Ceftamin_tumors.quantile([.25,.5,.75])
lowerq = quartiles[0.25]
upperq = quartiles[0.75]
iqr = upperq-lowerq

print(f"The lower quartile of treatment is: {lowerq}")
print(f"The upper quartile of temperatures is: {upperq}")
print(f"The interquartile range of temperatures is: {iqr}")
print(f"The the median of temperatures is: {quartiles[0.5]} ")

lower_bound = lowerq - (1.5*iqr)
upper_bound = upperq + (1.5*iqr)

print(f"Values below {lower_bound} could be outliers.")
print(f"Values above {upper_bound} could be outliers.")


# In[521]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
data_to_plot = [Capomulin_tumors, Ramicane_tumors, Infubinol_tumors, Ceftamin_tumors]
Regimen= ['Capomulin', 'Ramicane', 'Infubinol','Ceftamin']
fig1, ax1 = plt.subplots(figsize=(10, 5))
ax1.set_title('Tumor Volume at Selected Mouse',fontsize =20)
ax1.set_ylabel('Final Tumor Volume (mm3)',fontsize = 12)
ax1.set_xlabel('Drug Regimen',fontsize = 12)
ax1.boxplot(data_to_plot, labels=Regimen, widths = 0.6, patch_artist=True,vert=True)
plt.ylim(20, 70)

plt.savefig("../Images/box_plot.png", bbox_inches = "tight")
plt.show()


# ## Line and Scatter Plots

# In[477]:


# Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin
forline_df = Capomulin_df.loc[Capomulin_df["Mouse ID"] == "b128",:]
forline_df.head()
x_axis = forline_df["Timepoint"]
tumsiz = forline_df["Tumor Volume (mm3)"]
fig1, ax1 = plt.subplots(figsize=(7, 5))
plt.title('Capomulin treatmeant of mouse b128',fontsize =25)
plt.plot(x_axis, tumsiz,linewidth=2, markersize=15,marker="o",color="blue", label="Fahreneit")
plt.xlabel('Timepoint (Days)',fontsize =12)
plt.ylabel('Tumor Volume (mm3)',fontsize =12)


plt.savefig("../Images/line_graph.png", bbox_inches = "tight")
plt.show()


# In[593]:


# Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin regimen
fig1, ax1 = plt.subplots(figsize=(7, 5))
avg_capm_volume =Capomulin_df.groupby(['Mouse ID']).mean().round(3)
marker_size=10
plt.scatter(avg_capm_volume['Weight (g)'],avg_capm_volume['Tumor Volume (mm3)'],s=170, color="blue")
plt.title('Mouse Weight Vs Average Tumor Volume',fontsize =20)
plt.xlabel('Weight (g)',fontsize =12)
plt.ylabel('Averag Tumor Volume (mm3)',fontsize =12)

plt.savefig("../Images/line_graph.png", bbox_inches = "tight")
plt.show()


# ## Correlation and Regression

# In[594]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen
corr=round(st.pearsonr(avg_capm_volume['Weight (g)'],avg_capm_volume['Tumor Volume (mm3)'])[0],2)
print(f"The correlation between mouse weight and average tumor volume is {corr}")


# In[655]:


x_values = avg_capm_volume['Weight (g)']
y_values = avg_capm_volume['Tumor Volume (mm3)']
slope, int, r, p, std_err = st.linregress(mouse_weight, tumor_volume)
fit = slope * mouse_weight + int

print(f"slope:{slope}")
print(f"intercept:{int}")
print(f"rvalue (Correlation coefficient):{r}")
print(f"pandas (Correlation coefficient):{p}")
print(f"stderr:{std_err}")

line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(int,2))
print(line_eq)

fig1, ax1 = plt.subplots(figsize=(7, 5))
plt.scatter(x_values,y_values,s=175, color="blue")
plt.plot(mouse_weight,fit,"--")
plt.title('Regression Plot of Mouse Weight Versus Average Tumor Volume',fontsize =20)
plt.xlabel('Weight(g)',fontsize =12)
plt.ylabel('Average Tumore Volume (mm3)',fontsize =12)
ax1.annotate(line_eq, xy=(20, 40), xycoords='data',xytext=(0.8, 0.95), textcoords='axes fraction',horizontalalignment='right', verticalalignment='top',fontsize=30,color="red")

print(f"The r-squared is: {rvalue**2}")

plt.savefig("../Images/linear_regression.png", bbox_inches = "tight")
plt.show()


# In[ ]:




