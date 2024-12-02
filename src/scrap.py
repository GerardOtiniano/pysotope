import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('/Users/gerard/Desktop/UB/Projects/RAW/Emanda/dD Results/20231127_GAO_Emanda_FAME_d2H_run_2.csv')
df = df[df.type=="standard"]
#mask = (df['Identifier 1'].str.contains("C18") & df['Identifier 1'].str.contains("C24"))
#df = df[mask]
fig = plt.figure()
for i in df["Identifier 1"].unique():
    temp = df[df["Identifier 1"]== i]
    #temp = temp[temp.Component.isin(["C18", "C24"])]
    for j in ["C18", "C20", "C24", "C28"]:
        temp = temp[temp.Component == j]
        if temp is not None:
            plt.scatter(temp["Area All"], temp["d 2H/1H"]-temp["d 2H/1H"].mean())
plt.show()

# %%
# Import necessary libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.stats.weightstats import DescrStatsW
from statsmodels.regression.linear_model import WLS
from statsmodels.sandbox.regression.predstd import wls_prediction_std

# Set random seed for reproducibility
np.random.seed(42)

# Generate synthetic data
# Let's assume x ranges from 1 to 100
x = np.linspace(1, 100, 100)

# Define a true linear relationship: y = 2.5x + 10
true_slope = 2.5
true_intercept = 10
y_true = true_slope * x + true_intercept

# Add some noise to y
noise = np.random.normal(0, 25, size=x.shape)
y_noisy = y_true + noise

# Assume weights inversely proportional to variance (for WLS)
# Let's say higher x has higher variance
weights = 1 / (0.5 * x + 5)  # Example weights

# Create a DataFrame for convenience
data = pd.DataFrame({
    'x': x,
    'y': y_noisy,
    'weights': weights
})

# Display the first few rows
data.head()

# Define independent variable (with constant term for intercept)
X = sm.add_constant(data['x'])

# Define dependent variable
Y = data['y']

# Define weights
W = data['weights']

# Fit the WLS model
wls_model = WLS(Y, X, weights=W)
wls_results = wls_model.fit()

# Print the regression summary
print(wls_results.summary())

# Create a scatter plot of the original data
plt.figure(figsize=(10, 6))
plt.scatter(data['x'], data['y'], label='Noisy Data', alpha=0.6)

# Plot the regression line
plt.plot(data['x'], wls_results.predict(X), color='red', label='WLS Regression Line')

plt.title('Weighted Least Squares Regression on Original Data')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.grid(True)
plt.show()