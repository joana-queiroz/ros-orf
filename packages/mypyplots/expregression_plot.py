import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

def exponential_function(x, a, b, c):
        """Exponential model function."""
        return a * (np.exp(b * x) - 1) + c

def exponential_regression_fit(df, xvar, yvar, title=None,  p0=None):
    """Fits the exponential model to the data and returns optimal parameters and their covariance."""
    
    print(title)
    
    # Check if required columns exist in the DataFrame
    if xvar not in df.columns or yvar not in df.columns:
        print(f'The DataFrame must contain {xvar} and/or {yvar} columns.\n')
        return
    
    # Extract data
    xdata = df[xvar]
    ydata = df[yvar]
    
    # Use default title if not provided
    if title is None:
        title = 'Exponential regression'

    # Use default initial guesses if not provided
    if p0 is None:
        p0 = [1, 0.02, 0]
    
    # Fit the exponential model
    try:
        popt, pcov = curve_fit(exponential_function, xdata, ydata, p0=p0)
        perr = np.sqrt(np.diag(pcov))  # Standard deviations of parameters
    except RuntimeError:
        # Error handling at the 'exponential_regression_plot' level
        print(f"Optimal parameters could not be found. Returning None.\n")
        return None, None, None, None
    
    # If fitting failed, exit function
    if popt is None:
        return  

    # If fittin successful...
    # Extract optimal parameters, Ppedict y values using the fitted model and calculate R-squared
    a_opt, b_opt, c_opt = popt   
    y_pred = exponential_function(xdata, a_opt, b_opt, c_opt) 
    r2 = r2_score(ydata, y_pred)

    print(f'Optimal parameters: a = {a_opt}, b = {b_opt}, c = {c_opt}')
    print(f'Parameter uncertainties: {perr}')
    print(f'R-squared: {r2}\n')
    
    return a_opt, b_opt, c_opt, r2

    

def exponential_regression_plot(df, xvar, yvar, title, a_opt, b_opt, c_opt, r2, xrange=(0, 10000)):
    """Plots data and fitted exponential model."""
    
    # Extract data
    xdata = df[xvar]
    ydata = df[yvar]

    # Generate x values for plotting the fitted curve
    xrange = np.linspace(xrange[0], xrange[1], 500)  # Use 500 points for smooth curve

    plt.scatter(xdata, ydata, label='Data points')
    plt.plot(xrange, exponential_function(xrange, a_opt, b_opt, c_opt), color='red', label='Exponential Fit', linewidth=1)
    plt.text(0.05, 0.95, f'f(x)={a_opt:.4f} * ({b_opt:.4f})^x + {c_opt:.4f}\nR-squared = {r2:.4f}', 
             transform=plt.gca().transAxes, fontsize=12, verticalalignment='top',
             bbox=dict(facecolor='white', alpha=0.5))
    plt.xlabel(xvar)
    plt.ylabel(yvar)
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.show()