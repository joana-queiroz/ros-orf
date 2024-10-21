import seaborn as sns
import matplotlib.pyplot as plt

# Function to create the FacetGrid and colorbar
def heatmap_grid(data, xaxis, yaxis, heatval, subpvar, gridtitle, heatbartitle):
    
    # Find min and max Running Time
    vmin = data['Running Time'].min()
    vmax = data['Running Time'].max()
    

    # Function to create heatmap for each subset of data
    def draw_heatmap(data, vmin, vmax, **kwargs):
        # Pivot the data for heatmap
        heatmap_data = data.pivot_table(values=heatval, 
                                         index=yaxis, 
                                         columns=xaxis, 
                                         aggfunc='mean')  # Aggregate if there are multiple values

        heatmap_data = heatmap_data[::-1]

        # Draw the heatmap without individual colorbars
        sns.heatmap(heatmap_data, 
                     cmap='Greens',  # Choose a different palette if desired
                     annot=True, 
                     fmt=".1f", 
                     cbar=False,  # Disable individual colorbars
                     linewidths=.5,
                     vmin=vmin,  
                     vmax=vmax,  
                     **kwargs)
        # Set title for the heatmap
        subptitle = data[subpvar].iloc[0]
        plt.title(subptitle, fontweight='bold', fontsize=14, pad=15)
        # Set x and y labels with bold font
        plt.xlabel(xaxis, fontsize=10, labelpad=10)
        plt.ylabel(yaxis, fontsize=10, labelpad=10)

    # Create a FacetGrid for heatmaps separated by subpvar
    g = sns.FacetGrid(data, col=subpvar, height=3, aspect=1.0)
    g.figure.suptitle(f'{gridtitle}', fontsize=18, y=1.05)

    # Create a shared colorbar axis
    cbar_ax = g.figure.add_axes([1.05, 0.2, 0.02, 0.6])

    # Map the heatmap function to the FacetGrid
    g.map_dataframe(draw_heatmap, vmin=vmin, vmax=vmax)

    # Create a ScalarMappable for the colorbar
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap='Greens', norm=norm)
    plt.colorbar(sm, cax=cbar_ax, label=heatbartitle)

    # Show the plot
    plt.show()