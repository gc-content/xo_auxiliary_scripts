import pandas as pd
import numpy as np
import sys
from sklearn.linear_model import LogisticRegression

# Check for the input file argument
if len(sys.argv) != 2:
    print("Usage: python logit_model.py <input_file>")
    sys.exit(1)

input_file = sys.argv[1]

# Load data
data = pd.read_csv(input_file, sep="\t")

# Convert necessary columns to numeric
data['SwitchScore'] = pd.to_numeric(data['SwitchScore'], errors='coerce')
data['LocalSwitchScore'] = pd.to_numeric(data['LocalSwitchScore'], errors='coerce')

# Drop rows with missing values in relevant columns
data = data.dropna(subset=['SwitchScore', 'LocalSwitchScore'])

# Define logistic regression model
X = data[['SwitchScore', 'LocalSwitchScore']]
y = data['on_target']
log_reg_model = LogisticRegression()

# Fit the model
log_reg_model.fit(X, y)

# Extract coefficients and intercept
beta_0 = log_reg_model.intercept_[0]  # Intercept
beta_1 = log_reg_model.coef_[0][0]    # Coefficient for switch
beta_2 = log_reg_model.coef_[0][1]    # Coefficient for local

# Output coefficients in a format suitable for bash/awk
print(f"beta0={beta_0}")
print(f"beta1={beta_1}")
print(f"beta2={beta_2}")
