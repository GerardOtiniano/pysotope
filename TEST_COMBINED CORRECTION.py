#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 13:37:55 2025

@author: gerard
"""

# %%
import statsmodels.api as sm
import pandas as pd
import matplotlib.pyplot as plt


bl = pd.read_csv('/Users/gerard/Downloads/example_run2.csv')
bl = bl[bl['Identifier 1'].isin(['BL SEDIMENT'])]
bl = bl[bl['Component']=="N2"]


df = pd.read_csv('/Users/gerard/Downloads/drift3.csv')
df = df[df['Identifier 1'].isin(['ACETANILIDE', 'SORGHUM', 'WHEAT FLOUR', 'L-ALANINE'])]
df = df[df['Component']=="N2"]
# df['drift_rel'] = [0]*len(df)
for standard in ['ACETANILIDE', 'SORGHUM', 'WHEAT FLOUR']:
    df.loc[df['Identifier 1'] == standard, 'drift_rel'] = (
        df.loc[df['Identifier 1'] == standard, 'd 15N/14N'] -
        df.loc[df['Identifier 1'] == standard, 'd 15N/14N'].mean())
    
    df.loc[df['Identifier 1'] == standard, 'ampl_rel'] = (
        df.loc[df['Identifier 1'] == standard, 'Amt%'] -
        df.loc[df['Identifier 1'] == standard, 'Amt%'].mean())
    
fig = plt.figure()
for standard, marker in zip(['ACETANILIDE', 'SORGHUM', 'WHEAT FLOUR', 'L-ALANINE'],['o','d','s', 'X']):
    temp = df[df['Identifier 1']==standard]
    plt.scatter(temp["Seconds Since Start"], temp["d 15N/14N"], marker=marker, c=temp['Ampl. 28'], label=standard) #c=temp['ampl_rel']
    # plt.scatter(temp["Seconds Since Start"], temp["drift_rel"], label=standard)
plt.legend()
# plt.scatter(df['Seconds Since Start'], df['drift_rel'])
plt.show()

# %%
df = pd.read_csv('/Users/gerard/Downloads/vpdb.csv')
df = df[df['Identifier 1'].isin(['ACETANILIDE', 'SORGHUM', 'WHEAT FLOUR'])]
df = df[df['Component']=="N2"]
standard = 'ACETANILIDE'
for name in ['ACETANILIDE', 'SORGHUM', 'WHEAT FLOUR']:
    print(df.loc[df['Identifier 1'] == name, 'd 15N/14N'].mean())
    print(df.loc[df['Identifier 1'] == name, 'd 15N/14N'])
    print((df.loc[df['Identifier 1'] == name, 'd 15N/14N'])-df.loc[df['Identifier 1'] == name, 'd 15N/14N'].mean())
    print("---")

# %%
fig = plt.figure()
for standard in ['ACETANILIDE', 'SORGHUM', 'WHEAT FLOUR']:
    temp = df[df['Identifier 1']==standard]
    plt.scatter(temp["Seconds Since Start"], temp["d 15N/14N"], label=standard)
plt.legend()
# plt.scatter(df['Seconds Since Start'], df['drift_rel'])
plt.show()
# %%
test = df[df['Identifier 1'].isin(['ACETANILIDE', 'SORGHUM', 'WHEAT FLOUR'])]
# test = df[df['Identifier 1'].isin(['SORGHUM'])]
test = test[test['Component']=="N2"]
X = test[["Seconds Since Start", "Ampl. 28"]]
y = test[['d 15N/14N']]

# Drift
X = sm.add_constant(X)
results = sm.OLS(y, X).fit()
pred = results.predict(X)
relative = y.mean()

df = pd.read_csv('/Users/gerard/Downloads/vpdb.csv')
df = df[df['Identifier 1'].isin(['ACETANILIDE', 'SORGHUM', 'WHEAT FLOUR'])]
df = df[df['Component']=="N2"]
pred_x = df[["Seconds Since Start", "Ampl. 28"]]
pred_x = sm.add_constant(pred_x)
pred_x = results.predict(pred_x)

bl_drift = sm.add_constant(bl[["Seconds Since Start", "Ampl. 28"]])
bl_drift = bl["d 15N/14N"]-results.predict(bl_drift)-relative.values

# %%
# VPDB
df = pd.read_csv('/Users/gerard/Downloads/vpdb.csv')
df = df[df['Identifier 1'].isin(['ACETANILIDE', 'SORGHUM', 'WHEAT FLOUR'])]
df = df[df['Component']=="N2"]
x_v = pred_x
y = df['std']
x_v = sm.add_constant(x_v)
results = sm.OLS(y, x_v).fit()
pred = results.predict(x_v)

df = pd.read_csv('/Users/gerard/Downloads/vpdb.csv')
df = df[df['Identifier 1'].isin(['BL SEDIMENT'])]
df = df[df['Component']=="N2"]
bl_predict = sm.add_constant(df[''])

fig = plt.figure()
plt.scatter(pred_x, pred, c='k')

plt.show()


# %%
fig = plt.figure()
plt.scatter(df['Seconds Since Start'],df['d 15N/14N'],c='k')
plt.scatter(df['Seconds Since Start'],pred_x,c='red')
plt.show()