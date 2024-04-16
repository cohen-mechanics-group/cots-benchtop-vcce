# Brendan M Unikewicz, PhD Student
# Andre M Pincot, MSc Student
# Tal Cohen, Asc. Professor
# MIT, Dept. Mechanical Engineering
# MIT, Dept. Civil & Environmental Engineering
# Date of Creation: 03/26/2024
# Code Purpose: checking boundary effects for neo-Hookean material models

import numpy as np
import matplotlib.pyplot as plt

# Define Colors
font_size = 14
colors = [(0, 0.4470, 0.7410), (0.8500, 0.3250, 0.0980), 
          (0.9290, 0.6940, 0.1250), (0.4940, 0.1840, 0.5560), 
          (0.4660, 0.6740, 0.1880), (0.3010, 0.7450, 0.9330), (0.6350, 0.0780, 0.1840)]

# Start Boundary Effects
# Define stretch range
lambda_range = np.linspace(0.05, 50, 10000)

# Define B/A ratios
BA_ratios = [5.0, 10.0, 15.0, 20.0, 1e6] # Including infinity
BA_ratios_legend = ['5', '10', '15', '20', r'$\infty$']

# Calculate P/E for each B/A ratio
P_E = np.zeros((len(BA_ratios), len(lambda_range)))
for i, BA in enumerate(BA_ratios):
    if np.isinf(BA):
        # For B/A -> infinity, P/E is always 0
        P_E[i, :] = np.zeros(len(lambda_range))
    else:
        # Calculate lambda_b using the provided equation
        A_over_B = 1/BA # Since B/A is provided, we invert it
        lambda_b = (1 + (lambda_range**3 - 1)*(A_over_B**3))**(1/3)
        # Calculate P/E using the values of lambda_b and lambda
        P_E[i, :] = (1/6) * (lambda_b**-4 + 4*lambda_b**-1 - lambda_range**-4 - 4*lambda_range**-1)
    plt.plot(lambda_range, P_E[i, :], color=colors[i], linewidth=1.5)
 
asymptote = np.repeat(5/6, len(lambda_range))
plt.plot(lambda_range, asymptote, '--k', linewidth=1.5)

txt_AL = r'Asymptotic Limit: 5/6'
plt.text(4.75, 0.916, txt_AL, ha='right', fontdict={'size': 13, 'family': 'serif'})

plt.xlabel(r'Stretch ($\lambda$)', fontdict={'size': font_size, 'family': 'serif'})
plt.xlim([1, 5])
plt.ylabel(r'p/E', fontdict={'size': font_size, 'family': 'serif'})
plt.ylim([0, 1])
plt.show()

