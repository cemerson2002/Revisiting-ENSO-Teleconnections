import cartopy.mpl.ticker as cticker
import matplotlib.ticker as mticker
import numpy as np
import cartopy.crs as ccrs

def apply_multipanel_map_styling(fig, axes_3d, contour, bounds, col_labels=False, row_labels=False, label_font_size=6, colorbar_fraction=0.02):
    """
    Applies unified map styling, gridlines, custom black outer ticks, 
    row/column labels, and a vertical colorbar to a multipanel Cartopy grid.
    """
    if fig.get_layout_engine():
        fig.get_layout_engine().set(rect=[0.02, 0, 0.98, 1])
    

#    if row_labels:
    num_rows, num_cols = axes_3d.shape
    axes0=axes_3d[0]
#    else: 
#        num_rows = 1
#        num_cols = len(axes_3d)
#        axes0=axes_3d

    # 1. Row and Column Header Annotations
    if col_labels:
        for ax, col_title in zip(axes0, col_labels):
            ax.set_title(col_title, fontsize=label_font_size+5)

    if row_labels:
      for ax, row_title in zip(axes_3d[:, 0], row_labels): # add [:, 0] to axes_3d
          ax.text(-0.2, 0.5, row_title, transform=ax.transAxes, rotation='vertical',
                  va='center', ha='center', fontsize=label_font_size+5)

    # 2. Iterate through individual axes to format gridlines and custom tick boundaries
    for i in range(num_rows):
        for j in range(num_cols):
        
            if len(np.shape(axes_3d))>1:
                ax = axes_3d[i,j]
            else:
                ax = axes_3d[j]
            
            # Interior Gray Gridlines
            ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False, 
                         linewidth=0.3, color='gray', alpha=0.5)
            
            # Outer Tick Text Labels Config
            gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                              linewidth=0, color='black')
            gl.top_labels = False
            gl.right_labels = False
            gl.bottom_labels = (i == num_rows - 1)
            gl.left_labels = (j == 0)

            gl.xlabel_style = {'size': label_font_size, 'color': 'black'}
            gl.ylabel_style = {'size': label_font_size, 'color': 'black'}

            gl.xformatter = cticker.LongitudeFormatter(degree_symbol='')
            gl.yformatter = cticker.LatitudeFormatter(degree_symbol='')

            gl.xpadding = 4
            gl.ypadding = 4
            gl.xlines = False
            gl.ylines = False
            gl.xlocator = mticker.FixedLocator([-120, -60, 0, 60, 120])

            # Force Black Native Tick Stubs (without overriding text labels)
            ax.tick_params(axis='both', which='both', direction='out',
                           width=0.5, length=3, color='black',
                           bottom=(i == num_rows - 1), top=False,
                           left=(j == 0), right=False,
                           labelbottom=False, labelleft=False,
                           labeltop=False, labelright=False)

            if i == num_rows - 1:
                ax.set_xticks([-120, -60, 0, 60, 120], crs=ccrs.PlateCarree())
            if j == 0:
                ax.set_yticks(np.arange(-60, 61, 30), crs=ccrs.PlateCarree())

    # 3. Master Colorbar Construction
    if row_labels:
        cbar = fig.colorbar(contour, ax=axes_3d.ravel(), orientation='vertical', fraction=colorbar_fraction, pad=0.03, 
                        boundaries=bounds, ticks=bounds)
    else:
        cbar = fig.colorbar(contour, ax=axes_3d.ravel(), orientation='horizontal', fraction=colorbar_fraction, pad=0.03, 
                        boundaries=bounds, ticks=bounds)
    cbar.ax.tick_params(labelsize=label_font_size)
    cbar.set_label("Correlation", fontsize=label_font_size+2)

