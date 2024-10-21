import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plotting
import matplotlib.ticker as ticker

import matplotlib
print(matplotlib.get_backend())  # Check current backend
matplotlib.use('TkAgg')
#matplotlib.use('Agg') 





import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D





def plot_3d_map_interactive_save(csv_path, title):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_path)

    # Set 'W/K' as the index
    df.set_index('W/K', inplace=True)

    # Convert all data to numeric, coerce errors to NaN
    df = df.apply(pd.to_numeric, errors='coerce')

    # Multiply each cell by (w+1) of the corresponding row index
    for w in df.index:
        df.loc[w] = df.loc[w] * (w + 1)

    # Get x, y, z data
    x = df.columns.values.astype(float)
    y = df.index.values.astype(float)
    z = df.values

    # Create meshgrid
    X, Y = np.meshgrid(x, y)

    # Mask NaN values for plotting
    Z = np.ma.array(z, mask=np.isnan(z))

    # Create the figure and 3D axes
    fig = plt.figure(figsize=(14, 10))  # Increased figure size for better visibility
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface with increased alpha for better visibility
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.9)

    # Set labels with increased padding
    ax.set_xlabel('k', labelpad=20, fontsize=25)
    ax.set_ylabel('w', labelpad=20, fontsize=25)
    ax.set_zlabel('Expected density factor', labelpad=25, fontsize=25)
    ax.zaxis.set_rotate_label(False)  # Disable automatic rotation
    ax.set_zlabel('Expected density factor', rotation=90, labelpad=25, fontsize=25)  # Set custom rotation

    # make the font larger
    ax.tick_params(axis='both', which='major', labelsize=20)

    # Adjust the viewing angle for better label visibility
    ax.view_init(elev=30, azim=-60)  # Adjust elevation and azimuth as needed

    # add some padding for the tick labels
    #ax.xaxis.set_tick_params(pad=20)
    #ax.yaxis.set_tick_params(pad=20)
    ax.zaxis.set_tick_params(pad=15)

    # Manually adjust the margins to prevent cropping
   # plt.subplots_adjust(left=0.1, right=0.85, bottom=0.1, top=0.9)

    # Display the plot interactively
    plt.show()  # This will open an interactive window where you can adjust the view

    # After closing the interactive window, save the plot with the adjusted view
    fig.savefig('figures/' + title + '.png', bbox_inches='tight', dpi=400)
    plt.close(fig)  # Close the figure to free up memory


from scipy.interpolate import griddata
def plot_3d_map_2(csv_path, title):
     # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_path)

    # Set 'W/K' as the index
    df.set_index('W/K', inplace=True)

    # Convert all data to numeric, coerce errors to NaN
    df = df.apply(pd.to_numeric, errors='coerce')
    
    # Multiply each cell by (w+1) of the corresponding row index
    for w in df.index:
        df.loc[w] = df.loc[w] * (w + 1)

    # Get x, y, z data
    x = df.columns.values.astype(float)
    y = df.index.values.astype(float)
    z = df.values

    # Create meshgrid
    X, Y = np.meshgrid(x, y)

    # Mask NaN values for plotting
    Z = np.ma.array(z, mask=np.isnan(z))

    # Create the figure and 3D axes
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface with transparency (remove depthshade)
    surf = ax.plot_surface(X, Y, Z, cmap='twilight_shifted', edgecolor='none',
                           alpha=0.9)

    # Plot the line for w = k
    # Generate a range of values covering the data range
    line_range = np.linspace(max(x.min(), y.min()), min(x.max(), y.max()), num=200)
    line_k = line_range
    line_w = line_range

    # Flatten the grid data for interpolation
    points = np.array([X.flatten(), Y.flatten()]).T
    values = Z.flatten()

    # Remove masked points
    mask = ~Z.mask.flatten()
    points = points[mask]
    values = values[mask]

    # Prepare the line points
    line_points = np.array([line_k, line_w]).T

    # Interpolate z values along the line
    line_z = griddata(points, values, line_points, method='linear')

    # Add a small offset to z to plot the line above the surface
    offset = (np.nanmax(Z) - np.nanmin(Z)) * 0 # Increased to 5% of the z-range
    line_z += offset

    # Remove any NaN values resulting from interpolation
    valid_mask = ~np.isnan(line_z)
    line_k = line_k[valid_mask]
    line_w = line_w[valid_mask]
    line_z = line_z[valid_mask]

    # Plot the line (plot before the surface)
    ax.plot(line_k, line_w, line_z, color='red', linewidth=2, label='w=k')

    # Bring the line to the front
    ax.collections[-1].set_zorder(10)

    # Set labels and title with increased padding
    ax.set_xlabel('k', labelpad=15, fontsize=16)
    ax.set_ylabel('w', labelpad=15, fontsize=16)
    ax.set_zlabel('Expected density factor', labelpad=15, fontsize=16)

    # Add legend
    ax.legend()

    # Adjust the viewing angle for better label visibility
    ax.view_init(elev=20, azim=-60)

    # Adjust the plot margins to prevent cropping
    plt.tight_layout()

    # Save the plot
    plt.savefig('figures/' + title + '.png', bbox_inches='tight', dpi=400)

    # Show the plot
    plt.show()


def plot_heatmap(csv_path, title):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(csv_path)

    # Set 'W/K' as the index
    df.set_index('W/K', inplace=True)

    # Convert all data to numeric, coerce errors to NaN (handles missing values)
    df = df.apply(pd.to_numeric, errors='coerce')

    # Multiply each cell by (w+1) of the corresponding row index
    for w in df.index:
        df.loc[w] = df.loc[w] * (w + 1)

    # Create a custom colormap that displays NaN values as gray
    cmap = plt.get_cmap('viridis')
    cmap.set_bad(color='gray')

    # Plot the heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(df, cmap=cmap, annot=True, fmt=".4f", cbar_kws={'label': 'Expected density factor'}, annot_kws={"size": 8})

    # Set plot labels and title
    plt.xlabel('k')
    plt.ylabel('w')

    plt.savefig('figures/' + title + '.png', bbox_inches='tight', dpi=600)

    # Show the plot
    plt.show()


# def plot_k_extended(w, path, title, extension_k, show_density_factor=True):
#     if path.endswith('.tsv'):
#         df = pd.read_csv(path, sep='\t')
#     else:
#         df = pd.read_csv(path)

#     # Extract x values from the column headers (excluding the 'method' column)
#     x_values = df.columns[1:].astype(float)

#     # Make them int
#     x_values = x_values.astype(int)

#     k_values = x_values


#     # Initialize figure
#     plt.figure(figsize=(10, 6))

#     # Store handles for creating custom legend
#     legend_handles = []
#     legend_labels = []

#      # Calculate density factor for each column (excluding 'method')
#     for col in df.columns[1:]:
#         df[col] = df[col] * (w + 1)  # Multiply all values in the column by (W + 1)

#     # Plot each method as a separate line
#     for index, row in df.iterrows():
#         method_name = row['method']
#         y_values = row[1:]  # Exclude 'method' column to get y-values
#         if method_name == 'GreedyMini':
#             # Identify the extension point index
#             extension_index = np.where(x_values == extension_k)[0][0]
            
#             # Plot the segment before and including extension_k
#             line1, = plt.plot(
#                 x_values[:extension_index + 1], y_values[:extension_index + 1], 
#                 linestyle='-', label=method_name, linewidth=2
#             )

#             # Plot the segment after extension_k with a different line style
#             line2, = plt.plot(
#                 x_values[extension_index:], y_values[extension_index:], 
#                 linestyle='--', color=line1.get_color(), linewidth=2, label='_nolegend_'
#             )

#             # Plot X markers for points before extension_k
#             plt.plot(
#                 x_values[:extension_index], y_values[:extension_index], 
#                 marker='x', linestyle='', label='_nolegend_', markersize=6, color=line1.get_color()
#             )

#             # Special star marker for the extension_k point
#             plt.scatter(
#                 extension_k, y_values[extension_index], 
#                 color='red', s=100, marker='*', label='_nolegend_', zorder=5
#             )

#             # Add the first segment and extension point to the custom legend
#             legend_handles.append(line1)
#             legend_labels.append(method_name)
#             legend_handles.append(plt.Line2D([0], [0], color='red', marker='*', markersize=10, lw=0))
#             legend_labels.append('extension point')
#         else:
#             # Plot the line for other methods
#             line, = plt.plot(x_values, y_values, linestyle='-', label=method_name, marker='x', markersize=6)
            
#             # Get the color of the line and add it to the legend
#             line_color = line.get_color()
#             legend_handles.append(line)
#             legend_labels.append(method_name)

#     # Lower bound calculations
#     lower_bound_1 = np.ceil((w + k_values) / w) / (w + k_values)
#     k_values_new = np.array([find_first_mod_1_after_x(x, w) for x in k_values])
#     lower_bound_2 = np.ceil((w + k_values_new) / w) / (w + k_values_new)
#     lower_bound = np.maximum(lower_bound_1, lower_bound_2)

#     lower_bound = lower_bound * (w + 1)

#     # Plot lower bound with a dashed line
#     line, = plt.plot(x_values, lower_bound, linestyle='--', color='black', label='lower bound')
#     legend_handles.append(line)
#     legend_labels.append('lower bound')

#     # Plot formatting
#     plt.xlabel('k', fontsize=16)
#     if show_density_factor:
#         plt.ylabel('Expected density factor', fontsize=16)

#     # limit plot to y=2.2, keeping the lower limit unchanged if it is lower than current limit
#     current_ylim = plt.ylim()
#     if current_ylim[1] > 2.2:
#         plt.ylim(current_ylim[0], 2.2)

#     # Increase the size of the x and y axis ticks
#     plt.xticks(fontsize=16)
#     plt.yticks(fontsize=16)

#     ax = plt.gca()
#     ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
#     ax.tick_params(axis='x', which='minor', length=5, width=1)

#     plt.grid(True)

#     # Create a figures folder if it doesn't exist
#     if not os.path.exists('figures'):
#         os.makedirs('figures')

     
#     # CURRENT CONSTANT LIMIT:
#     plt.ylim(1.3, 2.05)   

#     # Save the plot without the legend
#     plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')

#     # Show the plot without the legend
#     plt.show()

#     # Create a separate figure for the legend
#     fig_legend = plt.figure(figsize=(10, 2))  # Adjust size for legend
#     plt.figlegend(legend_handles, legend_labels, loc='center', fontsize=12, ncol=3)  # 3-column legend

#     # Save the legend as a separate file
#     fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')

#     # Show the legend figure
#     plt.show()


# def plot_w_extended(k, path, title, extension_w, max_w_for_focus, show_density_factor=True):
#     if path.endswith('.tsv'):
#         df = pd.read_csv(path, sep='\t')
#     else:
#         df = pd.read_csv(path)

#     # Extract x values from the column headers (excluding the 'method' column)
#     x_values = df.columns[1:].astype(float)

#     # Make them int
#     x_values = x_values.astype(int)

#     w_values = x_values

#     # Calculate density factor for each column (excluding 'method')
#     for col in df.columns[1:]:
#         W = int(col)  # The first cell in the column is W
#         df[col] = df[col] * (W + 1)  # Multiply all values in the column by (W + 1)

#     # Calculate the lower bound
#     lower_bound_1 = np.ceil((w_values + k) / w_values) / (w_values + k)
#     new_k_values = np.array([find_first_mod_1_after_x(y, k) for y in x_values])
#     lower_bound_2 = np.ceil((w_values + new_k_values) / w_values) / (w_values + new_k_values)
#     lower_bound = np.maximum(lower_bound_1, lower_bound_2)

#     # Multiply the lower bound by (w + 1), where w = x_values
#     lower_bound = lower_bound * (x_values + 1)

#     #print(lower_bound[:10])

#     # First Plot: Up to max_w_for_focus with markers on every point
#     fig1, ax1 = plt.subplots(figsize=(10, 6))

#     # Filter x_values for the first plot
#     indices_focus = x_values <= max_w_for_focus
#     x_values_focus = x_values[indices_focus]
#     lower_bound_focus = lower_bound[indices_focus]


#     # Plot each method
#     legend_handles = []
#     legend_labels = []
#     for index, row in df.iterrows():
#         method_name = row['method']
#         y_values_full = row[1:].values.astype(float)
#         y_values = y_values_full[indices_focus]

#         if method_name == 'GreedyMini':
#             # Split line before and after extension_w
#             extension_index_focus = np.where(x_values_focus == extension_w)[0][0]

#             # Plot the line before extension_w
#             line1, = ax1.plot(
#                 x_values_focus[:extension_index_focus + 1], y_values[:extension_index_focus + 1],
#                 linestyle='-', label=method_name, linewidth=2
#             )

#             # Plot the line after extension_w with a different style (e.g., dashed)
#             line2, = ax1.plot(
#                 x_values_focus[extension_index_focus:], y_values[extension_index_focus:], 
#                 linestyle='--', color=line1.get_color(), linewidth=2, label='_nolegend_'
#             )

#             # Special star marker for the extension_w point if within focus range
#             ax1.scatter(
#                 extension_w, y_values[extension_index_focus],
#                 color='red', s=100, marker='*', label='_nolegend_', zorder=5
#             )

#             # Add to legend
#             legend_handles.append(line1)
#             legend_labels.append(method_name)
#             legend_handles.append(plt.Line2D([0], [0], color='red', marker='*', markersize=10, lw=0))
#             legend_labels.append('extension point')

#         else:
#             line, = ax1.plot(x_values_focus, y_values, linestyle='-', label=method_name)
#             legend_handles.append(line)
#             legend_labels.append(method_name)

#     # Plot the lower bound
#     ax1.plot(x_values_focus, lower_bound_focus, linestyle='--', color='black', label='lower bound')
#     legend_handles.append(plt.Line2D([0], [0], color='black', lw=2, linestyle='--'))
#     legend_labels.append('lower bound')

#     # Plot formatting
#     ax1.set_xlabel('w', fontsize=16)
#     if show_density_factor:
#         ax1.set_ylabel('Expected density factor', fontsize=16)
#     ax1.tick_params(axis='both', which='major', labelsize=14)
#     ax1.grid(True)

#     ax = plt.gca()
#     ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
#     ax.tick_params(axis='x', which='minor', length=5, width=1)

    
#     # CURRENT CONSTANT LIMIT:
#     plt.ylim(1.4, 2.01)

#     # Save the first plot
#     plt.savefig('figures/' + title + '_focus.svg', bbox_inches='tight')



#     # Second Plot: Zoom-in around the extension_w point
#     fig2, ax2 = plt.subplots(figsize=(10, 6))

#     # Define a zoom window around extension_w (e.g., +/- 5 units)
#     zoom_window = 5
#     indices_zoom = (x_values >= extension_w - zoom_window) & (x_values <= extension_w + zoom_window)
#     x_values_zoom = x_values[indices_zoom]
#     lower_bound_zoom = lower_bound[indices_zoom]

#     # Plot each method
#     for index, row in df.iterrows():
#         method_name = row['method']
#         y_values_zoom = row[1:].values[indices_zoom].astype(float)

#         if method_name == 'GreedyMini':
#             extension_index_zoom = np.where(x_values_zoom == extension_w)[0][0]

#             # Plot the line before extension_w
#             line1, = ax2.plot(
#                 x_values_zoom[:extension_index_zoom + 1], y_values_zoom[:extension_index_zoom + 1],
#                 linestyle='-', label=method_name, linewidth=2
#             )

#             # Plot the line after extension_w with a different style (e.g., dashed)
#             line2, = ax2.plot(
#                 x_values_zoom[extension_index_zoom:], y_values_zoom[extension_index_zoom:], 
#                 linestyle='--', color=line1.get_color(), linewidth=2, label='_nolegend_'
#             )

#             # Special star marker for the extension_w point
#             ax2.scatter(
#                 extension_w, y_values_zoom[extension_index_zoom],
#                 color='red', s=100, marker='*', label='_nolegend_', zorder=5
#             )
#         else:
#             ax2.plot(x_values_zoom, y_values_zoom, linestyle='-', label=method_name)

#     # Plot the lower bound
#     ax2.plot(x_values_zoom, lower_bound_zoom, linestyle='--', color='black', label='lower bound')

#     # Plot formatting
#     ax2.set_xlabel('w', fontsize=16)
#     if show_density_factor:
#         ax2.set_ylabel('Expected density factor', fontsize=16)
#     ax2.tick_params(axis='both', which='major', labelsize=14)
#     ax2.grid(True)

#     plt.ylim(1.5,1.7)

#     # Save the zoomed-in plot
#     plt.savefig('figures/' + title + '_zoom.svg', bbox_inches='tight')



#     # Create a separate figure for the legend
#     fig_legend = plt.figure(figsize=(10, 2))  # Adjust size for legend
#     plt.figlegend(legend_handles, legend_labels, loc='center', fontsize=12, ncol=3)  # 3-column legend

#     # Save the legend as a separate file
#     fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')

#     # Show the legend figure
#     plt.show()


def find_first_mod_1_after_x(x, w):
    remainder = x % w
    if remainder == 0:
        to_add = 1
    else:
        to_add = w - remainder + 1
    return x + to_add


# def plot_w_is_k(path, title, extension_w, show_density_factor=True):
#     if path.endswith('.tsv'):
#         df = pd.read_csv(path, sep='\t')
#     else:
#         df = pd.read_csv(path)

#     # Extract x values (w values) from the column headers (excluding the 'method' column)
#     x_values = df.columns[1:].astype(float)
#     x_values = x_values.astype(int)  # Convert to int

#     # Since k = w
#     k_values = x_values.copy()

#     # Initialize figure
#     plt.figure(figsize=(10, 6))

#     # Store handles for custom legend
#     legend_handles = []
#     legend_labels = []

#     # Plot each method
#     for index, row in df.iterrows():
#         method_name = row['method']
#         y_values = row[1:].astype(float)
#         # Multiply y-values by (k + 1)
#         y_values = y_values * (k_values + 1)

#         if method_name == 'GreedyMini':
#             # Plot the line without markers first and get the line color
#             line, = plt.plot(x_values, y_values, linestyle='-', label=method_name, linewidth=2)
#             line_color = line.get_color()

#             # Plot markers on all points except the extension_w point
#             all_indices = np.arange(len(x_values))
#             if extension_w in x_values:
#                 extension_index = np.where(x_values == extension_w)[0][0]
#                 other_indices = np.delete(all_indices, extension_index)
#                 plt.plot(x_values[other_indices], y_values[other_indices], marker='x', linestyle='', markersize=6, color=line_color)
#                 # Plot the star at the extension_w point
#                 plt.scatter(
#                     extension_w, y_values[extension_index],
#                     color='red', s=100, marker='*', zorder=5
#                 )
#                 # Add to legend
#                 legend_handles.append(line)
#                 legend_labels.append(method_name)
#                 legend_handles.append(plt.Line2D([0], [0], color='red', marker='*', markersize=10, lw=0))
#                 legend_labels.append('extension point')
#             else:
#                 # Plot markers on all points
#                 plt.plot(x_values, y_values, marker='x', linestyle='', markersize=6, color=line_color)
#                 # Add to legend
#                 legend_handles.append(line)
#                 legend_labels.append(method_name)
#         else:
#             # Plot other methods with markers
#             line, = plt.plot(x_values, y_values, linestyle='-', label=method_name, marker='x', markersize=6)
#             legend_handles.append(line)
#             legend_labels.append(method_name)

#     # lower bound Calculations for k = w
#     w = x_values
#     k = k_values

#     # Using the given formula for the lower bound
#     lower_bound_1 = np.ceil((w + k) / w) / (w + k)
#     x_values_2 = np.array([find_first_mod_1_after_x(k_i, w_i) for k_i, w_i in zip(k, w)])
#     lower_bound_2 = np.ceil((w + x_values_2) / w) / (w + x_values_2)
#     lower_bound = np.maximum(lower_bound_1, lower_bound_2)

#     # Multiply the lower bound by (k + 1)
#     lower_bound = lower_bound * (k + 1)

#     # Plot the lower bound
#     line, = plt.plot(x_values, lower_bound, linestyle='--', color='black', label='lower bound')
#     legend_handles.append(line)
#     legend_labels.append('lower bound')

#     # Plot formatting
#     plt.xlabel('w', fontsize=16)
#     if show_density_factor:
#         plt.ylabel('Expected density factor', fontsize=16)

#     # Increase the size of the x and y axis ticks
#     plt.xticks(fontsize=16)
#     plt.yticks(fontsize=16)

#     plt.grid(True)

#     # Create figures directory if it doesn't exist
#     if not os.path.exists('figures'):
#         os.makedirs('figures')

#     # limit plot to y=2.2, keeping the lower limit unchanged if it is lower than current limit
#     # current_ylim = plt.ylim()
#     # if current_ylim[1] > 2.2:
#     #     plt.ylim(current_ylim[0], 2.2)

#     # CURRENT CONSTANT LIMIT:
#     plt.ylim(1.53, 2.05)


#     ax = plt.gca()
#     ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
#     ax.tick_params(axis='x', which='minor', length=5, width=1)

#     # Save the plot without the legend
#     plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')

#     # Show the plot without the legend
#     plt.show()

#     # Create a separate figure for the legend
#     fig_legend = plt.figure(figsize=(10, 2))  # Adjust size for legend
#     plt.figlegend(legend_handles, legend_labels, loc='center', fontsize=12, ncol=3)  # 3-column legend

#     # Save the legend as a separate file
#     fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')

#     # Show the legend figure
#     plt.show()


# def plot_w_k_is_w_plus_1(path, title, extension_w, show_density_factor=True):
#     if path.endswith('.tsv'):
#         df = pd.read_csv(path, sep='\t')
#     else:
#         df = pd.read_csv(path)

#     # Extract x values (w values) from the column headers (excluding the 'method' column)
#     x_values = df.columns[1:].astype(float)
#     x_values = x_values.astype(int)  # Convert to int

#     # Since k = w + 1
#     w = x_values
#     k = w + 1

#     # Initialize figure
#     plt.figure(figsize=(10, 6))

#     # Store handles for custom legend
#     legend_handles = []
#     legend_labels = []

#     # Plot each method
#     for index, row in df.iterrows():
#         method_name = row['method']
#         y_values = row[1:].astype(float)
#         # Multiply y-values by (w + 1)
#         y_values = y_values * (w + 1)

#         if method_name == 'GreedyMini':
#             # Plot the line without markers first and get the line color
#             line, = plt.plot(w, y_values, linestyle='-', label=method_name, linewidth=2)
#             line_color = line.get_color()

#             # Plot markers on all points except the extension_w point
#             all_indices = np.arange(len(w))
#             if extension_w in w:
#                 extension_index = np.where(w == extension_w)[0][0]
#                 other_indices = np.delete(all_indices, extension_index)
#                 plt.plot(w[other_indices], y_values[other_indices], marker='x', linestyle='', markersize=6, color=line_color)
#                 # Plot the star at the extension_w point
#                 plt.scatter(
#                     extension_w, y_values[extension_index],
#                     color='red', s=100, marker='*', zorder=5
#                 )
#                 # Add to legend
#                 legend_handles.append(line)
#                 legend_labels.append(method_name)
#                 legend_handles.append(plt.Line2D([0], [0], color='red', marker='*', markersize=10, lw=0))
#                 legend_labels.append('extension point')
#             else:
#                 # Plot markers on all points
#                 plt.plot(w, y_values, marker='x', linestyle='', markersize=6, color=line_color)
#                 # Add to legend
#                 legend_handles.append(line)
#                 legend_labels.append(method_name)
#         else:
#             # Plot other methods with markers
#             line, = plt.plot(w, y_values, linestyle='-', label=method_name, marker='x', markersize=6)
#             legend_handles.append(line)
#             legend_labels.append(method_name)

#     # lower bound Calculations for k = w + 1
#     lower_bound_1 = np.ceil((w + k) / w) / (w + k)
#     k_values_2 = np.array([find_first_mod_1_after_x(k_i, w_i) for k_i, w_i in zip(k, w)])
#     lower_bound_2 = np.ceil((w + k_values_2) / w) / (w + k_values_2)
#     lower_bound = np.maximum(lower_bound_1, lower_bound_2)

#     # Multiply the lower bound by (w + 1)
#     lower_bound = lower_bound * (w + 1)

#     # Plot the lower bound
#     line, = plt.plot(w, lower_bound, linestyle='--', color='black', label='lower bound')
#     legend_handles.append(line)
#     legend_labels.append('lower bound')

#     # Plot formatting
#     plt.xlabel('w', fontsize=16)
#     if show_density_factor:
#         plt.ylabel('Expected density factor', fontsize=16)

#     # Increase the size of the x and y axis ticks
#     plt.xticks(fontsize=16)
#     plt.yticks(fontsize=16)

#     plt.grid(True)

#     # Create figures directory if it doesn't exist
#     if not os.path.exists('figures'):
#         os.makedirs('figures')

       
#     # CURRENT CONSTANT LIMIT:
#     plt.ylim(1.53, 2.05) 

#     ax = plt.gca()
#     ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
#     ax.tick_params(axis='x', which='minor', length=5, width=1)

#     # Save the plot without the legend
#     plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')

#     # Show the plot without the legend
#     plt.show()

#     # Create a separate figure for the legend
#     fig_legend = plt.figure(figsize=(10, 2))  # Adjust size for legend
#     plt.figlegend(legend_handles, legend_labels, loc='center', fontsize=12, ncol=3)  # 3-column legend

#     # Save the legend as a separate file
#     fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')

#     # Show the legend figure
#     plt.show()



# def plot_k_specific(w, path, title, show_density_factor=True):
#     if path.endswith('.tsv'):
#         df = pd.read_csv(path, sep='\t')
#     else:
#         df = pd.read_csv(path)

#     # Extract x values from the column headers (excluding the 'method' column)
#     x_values = df.columns[1:].astype(float)
    
#     # Make them int
#     x_values = x_values.astype(int)

#     # Calculate density factor for each column (excluding 'method')
#     for col in df.columns[1:]:
#         df[col] = df[col] * (w + 1)


#     # Plot each method as a separate line
#     plt.figure(figsize=(10, 6))
#     lines = []  # To store lines for legend
#     labels = []  # To store method names for legend
#     for index, row in df.iterrows():
#         method_name = row['method']
#         y_values = row[1:]  # Exclude 'method' column to get y-values
#         line, = plt.plot(x_values, y_values, linestyle='-', label=method_name, linewidth=2)

#         # Get the color of the line to apply it to the markers
#         line_color = line.get_color()

#         # Plot X marks for all points
#         plt.plot(x_values, y_values, marker='x', linestyle='', label='_nolegend_', markersize=6, color=line_color)

#         # Store lines for legend
#         lines.append(line)
#         labels.append(method_name)

#     # Update font size for axis labels
#     plt.xlabel('k', fontsize=16)
#     if show_density_factor:
#         plt.ylabel('Particular density factor', fontsize=16)

#     # Increase the size of the x and y axis ticks
#     plt.xticks(fontsize=16)
#     plt.yticks(fontsize=16)

#     plt.grid(True)

#     # Create a figures folder if it doesn't exist
#     if not os.path.exists('figures'):
#         os.makedirs('figures')

#     # limit plot to y=2.2, keeping the lower limit unchanged if it is lower than current limit
#     current_ylim = plt.ylim()
#     if current_ylim[1] > 2.2:
#         plt.ylim(current_ylim[0], 2.2)  

#     ax = plt.gca()
#     ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
#     ax.tick_params(axis='x', which='minor', length=5, width=1)  

#     # Save the plot without the legend
#     plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')
    
#     # Show the plot without the legend
#     plt.show()

#     # Create a separate figure for the legend
#     fig_legend = plt.figure(figsize=(10, 2))  # Adjust size for legend
#     plt.figlegend(lines, labels, loc='center', fontsize=12, ncol=3)  # 3-column legend

#     # Save the legend as a separate file
#     fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')

#     # Show the legend figure
#     plt.show()



# def plot_w_specific(w, path, title, show_density_factor=True):
#     if path.endswith('.tsv'):
#         df = pd.read_csv(path, sep='\t')
#     else:
#         df = pd.read_csv(path)

#     # Extract x values from the column headers (excluding the 'method' column)
#     x_values = df.columns[1:].astype(float)
#     x_values = x_values.astype(int)  # Make them int

#     # Calculate density factor for each column (excluding 'method')
#     for col in df.columns[1:]:
#         W = int(col)  # The first cell in the column is W
#         df[col] = df[col] * (W + 1)  # Multiply all values in the column by (W + 1)

#     # Plot each method as a separate line
#     plt.figure(figsize=(10, 6))
#     lines = []  # To store lines for legend
#     labels = []  # To store method names for legend
#     for index, row in df.iterrows():
#         method_name = row['method']
#         y_values = row[1:]  # Exclude 'method' column to get y-values
#         line, = plt.plot(x_values, y_values, linestyle='-', label=method_name, linewidth=2)

#         # Get the color of the line to apply it to the markers
#         line_color = line.get_color()

#         # Plot X marks for all points
#         plt.plot(x_values, y_values, marker='x', linestyle='', label='_nolegend_', markersize=6, color=line_color)

#         # Store lines for legend
#         lines.append(line)
#         labels.append(method_name)

#     # Update font size for axis labels
#     plt.xlabel('w', fontsize=16)
#     if show_density_factor:
#         plt.ylabel('Particular density factor', fontsize=16)

#     # Increase the size of the x and y axis ticks
#     plt.xticks(fontsize=16)
#     plt.yticks(fontsize=16)

#     plt.grid(True)

#     # Create a figures folder if it doesn't exist
#     if not os.path.exists('figures'):
#         os.makedirs('figures')


#     ax = plt.gca()
#     ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
#     ax.tick_params(axis='x', which='minor', length=5, width=1)

#     # Save the plot without the legend
#     plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')
    
#     # Show the plot without the legend
#     plt.show()

#     # Create a separate figure for the legend
#     fig_legend = plt.figure(figsize=(10, 2))  # Adjust size for legend
#     plt.figlegend(lines, labels, loc='center', fontsize=12, ncol=3)  # 2-column legend

#     # Save the legend as a separate file
#     fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')

#     # Show the legend figure
#     plt.show()


# NEW VERSION STARTS HERE

def plot_k_extended(w, path, title, extension_k, show_density_factor=True):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import pandas as pd
    import numpy as np
    import os

    if path.endswith('.tsv'):
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.read_csv(path)

    x_values = df.columns[1:].astype(int)
    k_values = x_values

    for col in df.columns[1:]:
        df[col] = df[col] * (w + 1)

    # Create figure and axes with fixed size and position
    fig = plt.figure(figsize=(6, 4))
    left_margin = 0.15
    bottom_margin = 0.15
    width = 0.8
    height = 0.8
    ax = fig.add_axes([left_margin, bottom_margin, width, height])

    legend_handles = []
    legend_labels = []

    for index, row in df.iterrows():
        method_name = row['method']
        y_values = row[1:]

        if method_name == 'GreedyMini':
            extension_index = np.where(x_values == extension_k)[0][0]

            line1, = ax.plot(
                x_values[:extension_index + 1], y_values[:extension_index + 1],
                linestyle='-', label=method_name, linewidth=2
            )

            line2, = ax.plot(
                x_values[extension_index:], y_values[extension_index:],
                linestyle='--', color=line1.get_color(), linewidth=2, label='_nolegend_'
            )

            ax.plot(
                x_values[:extension_index], y_values[:extension_index],
                marker='x', linestyle='', label='_nolegend_', markersize=6, color=line1.get_color()
            )

            ax.scatter(
                extension_k, y_values[extension_index],
                color='red', s=100, marker='*', label='_nolegend_', zorder=5
            )

            legend_handles.append(line1)
            legend_labels.append(method_name)
            legend_handles.append(plt.Line2D([0], [0], color='red', marker='*', markersize=10, lw=0))
            legend_labels.append('extension point')
        else:
            line, = ax.plot(x_values, y_values, linestyle='-', label=method_name, marker='x', markersize=6)
            legend_handles.append(line)
            legend_labels.append(method_name)

    lower_bound_1 = np.ceil((w + k_values) / w) / (w + k_values)
    k_values_new = np.array([find_first_mod_1_after_x(x, w) for x in k_values])
    lower_bound_2 = np.ceil((w + k_values_new) / w) / (w + k_values_new)
    lower_bound = np.maximum(lower_bound_1, lower_bound_2) * (w + 1)

    line, = ax.plot(x_values, lower_bound, linestyle='--', color='black', label='lower bound')
    legend_handles.append(line)
    legend_labels.append('lower bound')

    ax.set_xlabel('k', fontsize=12)
    if show_density_factor:
        ax.set_ylabel('Expected density factor', fontsize=12)

    ax.set_ylim(1.3, 2.05)
    ax.tick_params(axis='both', which='major', labelsize=10)

    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis='x', which='minor', length=5, width=1)
    ax.grid(True)

    if not os.path.exists('figures'):
        os.makedirs('figures')

    plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')
    plt.show()

    fig_legend = plt.figure(figsize=(6, 1))
    fig_legend.legend(legend_handles, legend_labels, loc='center', fontsize=10, ncol=3)

    fig_legend.patch.set_visible(False)
    fig_legend.canvas.draw()
    #fig_legend.axes[0].set_axis_off()

    fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')
    plt.show()

def plot_w_extended(k, path, title, extension_w, max_w_for_focus, show_density_factor=True):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import pandas as pd
    import numpy as np
    import os

    if path.endswith('.tsv'):
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.read_csv(path)

    x_values = df.columns[1:].astype(int)
    w_values = x_values

    for col in df.columns[1:]:
        W = int(col)
        df[col] = df[col] * (W + 1)

    lower_bound_1 = np.ceil((w_values + k) / w_values) / (w_values + k)
    new_k_values = np.array([find_first_mod_1_after_x(k, w) for w in w_values])
    lower_bound_2 = np.ceil((w_values + new_k_values) / w_values) / (w_values + new_k_values)
    lower_bound = np.maximum(lower_bound_1, lower_bound_2) * (w_values + 1)

    # First Plot
    fig = plt.figure(figsize=(6, 4))
    left_margin = 0.15
    bottom_margin = 0.15
    width = 0.8
    height = 0.8
    ax1 = fig.add_axes([left_margin, bottom_margin, width, height])

    indices_focus = x_values <= max_w_for_focus
    x_values_focus = x_values[indices_focus]
    lower_bound_focus = lower_bound[indices_focus]

    legend_handles = []
    legend_labels = []
    for index, row in df.iterrows():
        method_name = row['method']
        y_values_full = row[1:].astype(float)
        y_values = y_values_full[indices_focus]

        if method_name == 'GreedyMini':
            extension_index_focus = np.where(x_values_focus == extension_w)[0][0]

            line1, = ax1.plot(
                x_values_focus[:extension_index_focus + 1], y_values[:extension_index_focus + 1],
                linestyle='-', label=method_name, linewidth=2
            )

            line2, = ax1.plot(
                x_values_focus[extension_index_focus:], y_values[extension_index_focus:],
                linestyle='--', color=line1.get_color(), linewidth=2, label='_nolegend_'
            )

            ax1.scatter(
                extension_w, y_values[extension_index_focus],
                color='red', s=100, marker='*', label='_nolegend_', zorder=5
            )

            legend_handles.append(line1)
            legend_labels.append(method_name)
            legend_handles.append(plt.Line2D([0], [0], color='red', marker='*', markersize=10, lw=0))
            legend_labels.append('extension point')
        else:
            line, = ax1.plot(x_values_focus, y_values, linestyle='-', label=method_name)
            legend_handles.append(line)
            legend_labels.append(method_name)

    ax1.plot(x_values_focus, lower_bound_focus, linestyle='--', color='black', label='lower bound')
    legend_handles.append(plt.Line2D([0], [0], color='black', lw=2, linestyle='--'))
    legend_labels.append('lower bound')

    ax1.set_xlabel('w', fontsize=12)
    if show_density_factor:
        ax1.set_ylabel('Expected density factor', fontsize=12)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.grid(True)

    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax1.tick_params(axis='x', which='minor', length=5, width=1)
    ax1.set_ylim(1.4, 2.01)

    if not os.path.exists('figures'):
        os.makedirs('figures')

    plt.savefig('figures/' + title + '_focus.svg', bbox_inches='tight')
    plt.show()

    # The zoomed-in plot can remain unchanged or adjusted similarly if needed.



# **Mini-Plot: Zoomed around the extension point**
    fig_mini = plt.figure(figsize=(2, 2))
    ax_mini = fig_mini.add_axes([0, 0, 1, 1])  # Use full figure area for compactness

    # Define the zoom window
    x_min = extension_w - 3
    x_max = extension_w + 32
    indices_mini = (x_values >= x_min) & (x_values <= x_max)
    x_values_mini = x_values[indices_mini]
    lower_bound_mini = lower_bound[indices_mini]

    for index, row in df.iterrows():
        method_name = row['method']
        y_values_full = row[1:].astype(float)
        y_values = y_values_full[indices_mini]

        if method_name == 'GreedyMini':
            if extension_w in x_values_mini:
                extension_index_mini = np.where(x_values_mini == extension_w)[0][0]

                ax_mini.plot(
                    x_values_mini[:extension_index_mini + 1], y_values[:extension_index_mini + 1],
                    linestyle='-', linewidth=2
                )

                ax_mini.plot(
                    x_values_mini[extension_index_mini:], y_values[extension_index_mini:],
                    linestyle='--', color='blue', linewidth=2
                )

                ax_mini.scatter(
                    extension_w, y_values[extension_index_mini],
                    color='red', s=50, marker='*', zorder=5
                )
            else:
                ax_mini.plot(x_values_mini, y_values, linestyle='-', linewidth=2)
        else:
            ax_mini.plot(x_values_mini, y_values, linestyle='-', linewidth=2)

    ax_mini.plot(x_values_mini, lower_bound_mini, linestyle='--', color='black')

    # Set the desired limits
    ax_mini.set_xlim(x_min, x_max)
    ax_mini.set_ylim(1.58, 1.66)


    # Remove x and y labels and ticks for compactness
    #ax_mini.set_xticks([])
    #ax_mini.set_yticks([])
    #ax_mini.axis('off')  # Hide the axes frame

    # Optional: Add a grid if desired
    ax_mini.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Save the mini-plot
    plt.savefig('figures/' + title + '_mini.svg', bbox_inches='tight')
    plt.show()












    fig_legend = plt.figure(figsize=(6, 1))
    fig_legend.legend(legend_handles, legend_labels, loc='center', fontsize=10, ncol=3)

    fig_legend.patch.set_visible(False)
    fig_legend.canvas.draw()
    #fig_legend.axes[0].set_axis_off()

    fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')
    plt.show()

def plot_w_is_k(path, title, extension_w, show_density_factor=True):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import pandas as pd
    import numpy as np
    import os

    if path.endswith('.tsv'):
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.read_csv(path)

    x_values = df.columns[1:].astype(int)
    k_values = x_values.copy()

    # Create figure and axes with fixed size and position
    fig = plt.figure(figsize=(6, 4))
    left_margin = 0.15
    bottom_margin = 0.15
    width = 0.8
    height = 0.8
    ax = fig.add_axes([left_margin, bottom_margin, width, height])

    legend_handles = []
    legend_labels = []

    for index, row in df.iterrows():
        method_name = row['method']
        y_values = row[1:].astype(float) * (k_values + 1)

        if method_name == 'GreedyMini':
            line, = ax.plot(x_values, y_values, linestyle='-', label=method_name, linewidth=2)
            line_color = line.get_color()

            all_indices = np.arange(len(x_values))
            if extension_w in x_values:
                extension_index = np.where(x_values == extension_w)[0][0]
                other_indices = np.delete(all_indices, extension_index)
                ax.plot(x_values[other_indices], y_values[other_indices],
                        marker='x', linestyle='', markersize=6, color=line_color)
                ax.scatter(
                    extension_w, y_values[extension_index],
                    color='red', s=100, marker='*', zorder=5
                )
                legend_handles.append(line)
                legend_labels.append(method_name)
                legend_handles.append(plt.Line2D([0], [0], color='red',
                                                 marker='*', markersize=10, lw=0))
                legend_labels.append('extension point')
            else:
                ax.plot(x_values, y_values, marker='x', linestyle='',
                        markersize=6, color=line_color)
                legend_handles.append(line)
                legend_labels.append(method_name)
        else:
            line, = ax.plot(x_values, y_values, linestyle='-',
                            label=method_name, marker='x', markersize=6)
            legend_handles.append(line)
            legend_labels.append(method_name)

    w = x_values
    k = k_values

    lower_bound_1 = np.ceil((w + k) / w) / (w + k)
    x_values_2 = np.array([find_first_mod_1_after_x(k_i, w_i) for k_i, w_i in zip(k, w)])
    lower_bound_2 = np.ceil((w + x_values_2) / w) / (w + x_values_2)
    lower_bound = np.maximum(lower_bound_1, lower_bound_2) * (k + 1)

    line, = ax.plot(x_values, lower_bound, linestyle='--', color='black', label='lower bound')
    legend_handles.append(line)
    legend_labels.append('lower bound')

    ax.set_xlabel('w', fontsize=12)
    if show_density_factor:
        ax.set_ylabel('Expected density factor', fontsize=12)

    ax.set_ylim(1.53, 2.05)
    ax.tick_params(axis='both', which='major', labelsize=10)

    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis='x', which='minor', length=5, width=1)
    ax.grid(True)

    if not os.path.exists('figures'):
        os.makedirs('figures')

    plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')
    plt.show()

    fig_legend = plt.figure(figsize=(6, 1))
    fig_legend.legend(legend_handles, legend_labels, loc='center', fontsize=10, ncol=3)

    fig_legend.patch.set_visible(False)
    fig_legend.canvas.draw()
    #fig_legend.axes[0].set_axis_off()

    fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')
    plt.show()

def plot_w_k_is_w_plus_1(path, title, extension_w, show_density_factor=True):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import pandas as pd
    import numpy as np
    import os

    if path.endswith('.tsv'):
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.read_csv(path)

    x_values = df.columns[1:].astype(int)
    w = x_values
    k = w + 1

    # Create figure and axes with fixed size and position
    fig = plt.figure(figsize=(6, 4))
    left_margin = 0.15
    bottom_margin = 0.15
    width = 0.8
    height = 0.8
    ax = fig.add_axes([left_margin, bottom_margin, width, height])

    legend_handles = []
    legend_labels = []

    for index, row in df.iterrows():
        method_name = row['method']
        y_values = row[1:].astype(float) * (w + 1)

        if method_name == 'GreedyMini':
            line, = ax.plot(w, y_values, linestyle='-', label=method_name, linewidth=2)
            line_color = line.get_color()

            all_indices = np.arange(len(w))
            if extension_w in w:
                extension_index = np.where(w == extension_w)[0][0]
                other_indices = np.delete(all_indices, extension_index)
                ax.plot(w[other_indices], y_values[other_indices],
                        marker='x', linestyle='', markersize=6, color=line_color)
                ax.scatter(
                    extension_w, y_values[extension_index],
                    color='red', s=100, marker='*', zorder=5
                )
                legend_handles.append(line)
                legend_labels.append(method_name)
                legend_handles.append(plt.Line2D([0], [0], color='red',
                                                 marker='*', markersize=10, lw=0))
                legend_labels.append('extension point')
            else:
                ax.plot(w, y_values, marker='x', linestyle='',
                        markersize=6, color=line_color)
                legend_handles.append(line)
                legend_labels.append(method_name)
        else:
            line, = ax.plot(w, y_values, linestyle='-',
                            label=method_name, marker='x', markersize=6)
            legend_handles.append(line)
            legend_labels.append(method_name)

    lower_bound_1 = np.ceil((w + k) / w) / (w + k)
    k_values_2 = np.array([find_first_mod_1_after_x(k_i, w_i) for k_i, w_i in zip(k, w)])
    lower_bound_2 = np.ceil((w + k_values_2) / w) / (w + k_values_2)
    lower_bound = np.maximum(lower_bound_1, lower_bound_2) * (w + 1)

    line, = ax.plot(w, lower_bound, linestyle='--', color='black', label='lower bound')
    legend_handles.append(line)
    legend_labels.append('lower bound')

    ax.set_xlabel('w', fontsize=12)
    if show_density_factor:
        ax.set_ylabel('Expected density factor', fontsize=12)

    ax.set_ylim(1.53, 2.05)
    ax.tick_params(axis='both', which='major', labelsize=10)

    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis='x', which='minor', length=5, width=1)
    ax.grid(True)

    if not os.path.exists('figures'):
        os.makedirs('figures')

    plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')
    plt.show()

    fig_legend = plt.figure(figsize=(6, 1))
    fig_legend.legend(legend_handles, legend_labels, loc='center', fontsize=10, ncol=3)

    fig_legend.patch.set_visible(False)
    fig_legend.canvas.draw()
    #fig_legend.axes[0].set_axis_off()

    fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')
    plt.show()

def plot_k_specific(w, path, title, show_density_factor=True):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import pandas as pd
    import numpy as np
    import os

    if path.endswith('.tsv'):
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.read_csv(path)

    x_values = df.columns[1:].astype(int)

    for col in df.columns[1:]:
        df[col] = df[col] * (w + 1)

    # Create figure and axes with fixed size and position
    fig = plt.figure(figsize=(6, 4))
    left_margin = 0.15
    bottom_margin = 0.15
    width = 0.8
    height = 0.8
    ax = fig.add_axes([left_margin, bottom_margin, width, height])

    lines = []
    labels = []
    for index, row in df.iterrows():
        method_name = row['method']
        y_values = row[1:]
        line, = ax.plot(x_values, y_values, linestyle='-', label=method_name, linewidth=2)
        line_color = line.get_color()
        ax.plot(x_values, y_values, marker='x', linestyle='', label='_nolegend_', markersize=6, color=line_color)
        lines.append(line)
        labels.append(method_name)

    ax.set_xlabel('k', fontsize=12)
    if show_density_factor:
        ax.set_ylabel('Particular density factor', fontsize=12)

    current_ylim = ax.get_ylim()
    if current_ylim[1] > 2.2:
        ax.set_ylim(current_ylim[0], 2.2)

    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis='x', which='minor', length=5, width=1)
    ax.grid(True)

    if not os.path.exists('figures'):
        os.makedirs('figures')

    plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')
    plt.show()

    fig_legend = plt.figure(figsize=(6, 1))
    fig_legend.legend(lines, labels, loc='center', fontsize=10, ncol=3)

    fig_legend.patch.set_visible(False)
    fig_legend.canvas.draw()
    #fig_legend.axes[0].set_axis_off()

    fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')
    plt.show()

def plot_w_specific(w, path, title, show_density_factor=True):
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import pandas as pd
    import numpy as np
    import os

    if path.endswith('.tsv'):
        df = pd.read_csv(path, sep='\t')
    else:
        df = pd.read_csv(path)

    x_values = df.columns[1:].astype(int)

    for col in df.columns[1:]:
        W = int(col)
        df[col] = df[col] * (W + 1)

    # Create figure and axes with fixed size and position
    fig = plt.figure(figsize=(6, 4))
    left_margin = 0.15
    bottom_margin = 0.15
    width = 0.8
    height = 0.8
    ax = fig.add_axes([left_margin, bottom_margin, width, height])

    lines = []
    labels = []
    for index, row in df.iterrows():
        method_name = row['method']
        y_values = row[1:]
        line, = ax.plot(x_values, y_values, linestyle='-', label=method_name, linewidth=2)
        line_color = line.get_color()
        ax.plot(x_values, y_values, marker='x', linestyle='', label='_nolegend_', markersize=6, color=line_color)
        lines.append(line)
        labels.append(method_name)

    ax.set_xlabel('w', fontsize=12)
    if show_density_factor:
        ax.set_ylabel('Particular density factor', fontsize=12)

    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis='x', which='minor', length=5, width=1)
    ax.grid(True)

    if not os.path.exists('figures'):
        os.makedirs('figures')

    plt.savefig('figures/' + title + '_without_legend.svg', bbox_inches='tight')
    plt.show()

    fig_legend = plt.figure(figsize=(6, 1))
    fig_legend.legend(lines, labels, loc='center', fontsize=10, ncol=3)

    fig_legend.patch.set_visible(False)
    fig_legend.canvas.draw()
    #fig_legend.axes[0].set_axis_off()

    fig_legend.savefig('figures/' + title + '_legend.svg', bbox_inches='tight')
    plt.show()

