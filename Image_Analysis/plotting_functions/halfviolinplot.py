import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import warnings


def half_violin_box_plot(excel_path, sheet_id=0, nrows=None, skip=0, columns=None, column_names=None, Y_Value="Y Value",
                                  save_path=None, figsize=(4, 8)):
    """
    Create a paired plot with half-violin, box and strip plots from a single Excel file.

    Returns:
    --------
    fig, ax : matplotlib figure and axes objects
    """

    # Filter out the specific FutureWarning
    warnings.filterwarnings('ignore', category=FutureWarning,
                           message="Passing `palette` without assigning `hue` is deprecated")

    # Set style
    sns.set_style("ticks", {'axes.grid': False})

    if columns is None:
        df = pd.read_excel(excel_path, sheet_name=sheet_id, skiprows=skip, nrows=nrows)
    else:
        df = pd.read_excel(excel_path, sheet_name=sheet_id, usecols=columns, skiprows=skip, nrows=nrows)


    df = df.dropna(how='all', axis=0)
    print(df)

    use_columns = df.columns.tolist()


    if column_names is None:
        condition_names = use_columns
    else:
        condition_names = column_names


    fig, ax = plt.subplots(figsize=figsize)

    # Create DataFrame for plotting with multiple columns
    data_parts = []
    for i, col in enumerate(use_columns):
        data_parts.append(pd.DataFrame({
            'Value': df[col],
            'Condition': [condition_names[i] if column_names else col] * len(df)
        }))

    data = pd.concat(data_parts, axis=0).reset_index(drop=True)

    # Create visible and invisible datasets for violin plot
    data_visible = pd.DataFrame({
        'Value': pd.concat([df[col] for col in use_columns], axis=0),
        'Condition': [name for name in condition_names for _ in range(len(df))],
        'Plotting': ['Visible'] * (len(df) * len(use_columns))
    })

    data_invisible = pd.DataFrame({
        'Value': pd.concat([df[col] for col in use_columns], axis=0),
        'Condition': [name for name in condition_names for _ in range(len(df))],
        'Plotting': ['Invisible'] * (len(df) * len(use_columns))
    })

    combined_data_visible = pd.concat([data_visible, data_invisible], axis=0).reset_index(drop=True)
    flierprops = {'marker': 'o', 'markerfacecolor': 'none', 'markeredgecolor': 'm',
                  'markersize': 6, 'linestyle': 'none'}

    # Set color palette - extended for multiple conditions
    num_conditions = len(condition_names)
    if num_conditions == 2:
        palet = sns.color_palette(["lightsteelblue", "white"])
        palet2 = sns.color_palette(["gainsboro", "paleturquoise"])
    else:
        # Generate a palette for more than 2 conditions
        palet = sns.color_palette(["lightsteelblue", "white"])
        palet2 = sns.color_palette("pastel", num_conditions)

    # Create plots
    sns.violinplot(x='Condition', y='Value', hue='Plotting', data=combined_data_visible,
                   split=True, inner=None, palette=palet, linewidth=0, alpha=0.59,
                   legend=False, ax=ax, width=1.1)

    # Modified boxplot call 
    sns.boxplot(x='Condition', y='Value', hue='Condition', data=data, width=0.20, linewidth=1.8,
                palette=palet2, flierprops=flierprops, legend=False, ax=ax)

    sns.stripplot(x='Condition', y='Value', data=data, color='r', size=6,
                  alpha=0.3, jitter=0.02, ax=ax)

    # Customize plot
    sns.despine()
    ax.set_xlabel("")
    ax.set_ylabel(Y_Value, fontsize=16)
    ax.tick_params(axis='both', labelsize=16)
    plt.tick_params(axis='both', which='major', length=8, width=0.5)
    
    # plt.ylim(0, 0.5)
    # plt.yticks(np.arange(0, 0.51, 0.1))  # Longer and thicker ticks
    
    plt.tight_layout()

    # Save if path provided
    if save_path:
        plt.savefig(save_path)

    return fig, ax


def main():
    """Main function to run when script is executed directly."""

    import matplotlib
    matplotlib.use('TkAgg')  # Using TkAgg backend instead of %matplotlib tk
    matplotlib.rcParams['font.family'] = 'Arial'


    fig, ax = half_violin_box_plot(
        '/home/tatsatb/data.xlsx',
        columns=[0, 1],  
        column_names=['-DOX', '+DOX'],
        Y_Value ="Migration Speed (μm/min)",
        sheet_id=0
    )
    plt.show()



if __name__ == "__main__":
    main()
