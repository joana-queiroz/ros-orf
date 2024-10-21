import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def linear_regression_plot(df, xvar, yvar, title):
    """
    Perform linear regression on the DataFrame and plot the results.
    
    """
    
    # Check if required columns exist
    if xvar not in df.columns or yvar not in df.columns:
        print(f'The DataFrame must contain {xvar} and/or {yvar} columns.')
        return None

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(df[xvar], df[yvar])
    
    # Generate linear fit line
    df['Y_linear_fit'] = slope * df['Length of sequences'] + intercept
    
    # Plot linear fit
    plt.scatter(df[xvar], df[yvar], label='Data points',)
    plt.plot(df[xvar], df['Y_linear_fit'], color='red', label='Linear Fit', linewidth=1)
    plt.text(0.05, 0.95, f'f(x) = {slope:.4f}x + {intercept:.4f}\nR-squared = {r_value**2:.4f}', transform=plt.gca().transAxes, 
             fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    plt.xlabel(xvar)
    plt.ylabel(yvar)
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.show()

    return