# Preparing the data for the heatmap using calculated values from equation 7
import pandas as pd
import matplotlib.pyplot as plt


## Birthweight
# Creating a DataFrame from the data
df = pd.read_csv('../results/selectionBoundsBwt.csv', sep=';', decimal=',', thousands='.')
print(df.head())
#df = pd.DataFrame(data)

# Preparing the data for the heatmap
x_unique = sorted(df['tau_m'].unique())
y_unique = sorted(df['tau_A'].unique())
z_matrix = df.pivot(index='tau_A', columns='tau_m', values='DELTA').values

# Creating the heatmap
plt.figure(figsize=(10, 8))
plt.imshow(z_matrix, cmap='viridis', aspect='auto', origin='lower', extent=[min(x_unique), max(x_unique), min(y_unique), max(y_unique)])
cbar = plt.colorbar(label=r'$\Delta$ (Birthweight)')
cbar.set_label(r'$\Delta$ (Birthweight)', fontsize=16)  
cbar.ax.tick_params(labelsize=14)


plt.xlabel(r'$\tau^A$', fontsize=20)
plt.ylabel(r'$\tau^M$', fontsize=20)

# Increase axis values font size
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.annotate(
    'East et al. (2023)',                  # Text to display
    xy=(71, 71),                      # Coordinates for the point
    xytext=(100, 100),                # Position of the text
    arrowprops=dict(facecolor='white', shrink=0.05),  # Arrow properties
    fontsize=12, color='white', bbox=dict(boxstyle="round,pad=0.3", edgecolor='white', facecolor='black', alpha=0.5)
)

plt.savefig('../results/selection_BWT.pdf', format='pdf', bbox_inches='tight')

plt.clf()

## VLBW
# Creating a DataFrame from the data
df = pd.read_csv('../data/selectionBoundsVLBW.csv', sep=';', decimal=',', thousands='.')
print(df.head())
#df = pd.DataFrame(data)

# Preparing the data for the heatmap
x_unique = sorted(df['tau_m'].unique())
y_unique = sorted(df['tau_A'].unique())
z_matrix = df.pivot(index='tau_A', columns='tau_m', values='DELTA').values

# Creating the heatmap
plt.figure(figsize=(10, 8))
plt.imshow(z_matrix, cmap='viridis', aspect='auto', origin='lower', extent=[min(x_unique), max(x_unique), min(y_unique), max(y_unique)])
cbar = plt.colorbar(label=r'$\Delta$ (Proportion VLBW)')
cbar.set_label(r'$\Delta$ (Proportion VLBW)', fontsize=16)  
cbar.ax.tick_params(labelsize=14)
plt.xlabel(r'$\tau^A$', fontsize=20)
plt.ylabel(r'$\tau^M$', fontsize=20)

# Increase axis values font size
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.annotate(
    'East et al. (2023)',                  # Text to display
    xy=(-0.012, -0.012),                      # Coordinates for the point
    xytext=(-0.024, -0.027),                # Position of the text
    arrowprops=dict(facecolor='white', shrink=0.05),  # Arrow properties
    fontsize=12, color='white', bbox=dict(boxstyle="round,pad=0.3", edgecolor='white', facecolor='black', alpha=0.5)
)
plt.savefig('../data/selection_VLBW.pdf', format='pdf', bbox_inches='tight')

