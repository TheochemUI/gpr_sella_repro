import matplotlib.pyplot as plt
import re
import click
import gzip
import numpy as np
from matplotlib.ticker import PercentFormatter
import cmcrameri.cm as cmc
from pathlib import Path
from plotly.graph_objects import Sankey, Figure

def plot_sankey(all_iterations, all_opt_times, all_max_values, plot_title):
    """Generates a Sankey diagram showing optimization time and max value for each iteration."""

    # Create labels for Sankey diagram
    labels = [f"Iteration {i}" for i in all_iterations] + ["End"]
    source = list(range(len(all_iterations)))
    target = [len(all_iterations)] * len(all_iterations)  # All iterations link to "End"
    values = all_opt_times

    # Add max values as text to labels
    for i, label in enumerate(labels[:-1]):
        labels[i] = f"{label}\nTime: {all_opt_times[i]:.2f}s\nMax: {all_max_values[i]:.2f}"

    # Create Sankey diagram
    fig = Figure(data=[Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            color="blue"  # Customize node color if needed
        ),
        link=dict(
            source=source,
            target=target,
            value=values,
            color="lightblue"  # Customize link color if needed
        ))])

    fig.update_layout(title_text=f"{plot_title} - Sankey Diagram of Optimization Steps",
                      font_size=10)
    fig.show()

@click.command()
@click.argument('logfile_path', type=click.Path(exists=True))
@click.option('--cumulative', is_flag=True, help='Show cumulative optimization time as a filled area.')
@click.option('--sankey', is_flag=True, help='Generate a Sankey diagram of optimization steps.')
def plot_optimization_progress(logfile_path, cumulative, sankey):
    """
    Plots optimization progress data (optimization time, current max, total time)
    from a log file as a function of outer relaxation iteration.

    LOGFILE_PATH: Path to the log file.
    """

    logfile_path = Path(logfile_path)

    if logfile_path.suffix == '.gz':
        with gzip.open(logfile_path, "rb") as f:
            fcontents = f.read().decode('utf-8')
    else:
        with open(logfile_path, 'r') as f:
            fcontents = f.read()

    # Regular expressions to find the relevant data
    re_otime = re.compile(r'optimize time:\s*(\d+\.\d+)s')  # Optimization time
    re_max = re.compile(r'Current Max is\s*(\d+\.\d+)')  # Current Max
    re_iter = re.compile(r'Outer relaxation iteration:\s*(\d+)')  # Outer iteration number
    re_total_time = re.compile(r'real\s*(\d+\.\d+)\s*seconds')  # Total time
    re_title = re.compile(r'gprd/(.*)') # Title from directory structure

    # Extract data using the regular expressions
    all_opt_times = [float(x) for x in re_otime.findall(fcontents)]
    all_max_values = [float(x) for x in re_max.findall(fcontents)]
    all_iterations = [int(x) for x in re_iter.findall(fcontents)]
    total_time_matches = re_total_time.findall(fcontents)
    try:
        title_matches = re_title.findall(fcontents)[0].split('/')
    except IndexError:
        title_matches = None

    # Extract total time (take the last match if multiple lines are present)
    total_time = float(total_time_matches[-1]) if total_time_matches else None

    # Extract title information
    if title_matches:
        spin, mol_id = title_matches
        plot_title = f"GPRD {spin.capitalize()} {mol_id}"
    else:
        plot_title = "Optimization Progress"

    # Check if data lengths match
    if not (len(all_iterations) == len(all_opt_times) == len(all_max_values)):
        print("Warning: Data lengths do not match. Plotting might be inaccurate.")
        print(f"Iterations: {len(all_iterations)}, Opt Times: {len(all_opt_times)}, Max Values: {len(all_max_values)}")

        # Find minimum length for plotting
        min_len = min(len(all_iterations), len(all_opt_times), len(all_max_values))
        all_iterations = all_iterations[:min_len]
        all_opt_times = all_opt_times[:min_len]
        all_max_values = all_max_values[:min_len]

    # Generate Sankey diagram if requested
    if sankey:
        plot_sankey(all_iterations, all_opt_times, all_max_values, plot_title)
        return

    # Create the main plot
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Use Batlow colormap
    plt.rcParams['axes.prop_cycle'] = plt.cycler('color', cmc.batlow.colors)

    # Plot optimization times
    color = 'tab:blue'
    ax1.set_xlabel('Outer Relaxation Iteration')
    ax1.set_ylabel('Optimization Time (s)', color=color)
    if cumulative and total_time is not None:
        # Calculate cumulative optimization time
        cumulative_opt_times = np.cumsum(all_opt_times)
        # Normalize by total time to get percentages
        cumulative_opt_times_percentage = 100 * cumulative_opt_times / total_time
        # Fill the area under the cumulative optimization time curve
        ax1.fill_between(all_iterations, cumulative_opt_times_percentage, color=color, alpha=0.5)

        # Format the y-axis as percentages
        ax1.set_ylim(0, 100)
        ax1.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
        ax1.set_ylabel('Cumulative Optimization Time / Total Time (%)', color=color)
    else:
        ax1.plot(all_iterations, all_opt_times, color=color, marker='o')

    ax1.tick_params(axis='y', labelcolor=color)

    # Create a second y-axis for current max
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Current Max', color=color)
    ax2.plot(all_iterations, all_max_values, color=color, marker='x')
    ax2.tick_params(axis='y', labelcolor=color)

    # Add total time and optimization time difference to the plot (if available)
    if total_time is not None:
        total_opt_time = sum(all_opt_times)
        time_diff = total_time - total_opt_time

        plt.text(0.95, 0.95, f"Total Time: {(total_time/60):.2f} min\n"
                 f"Opt Time: {(total_opt_time/60):.2f} min\n"
                 f"Difference: {(time_diff/60):.2f} min",
                 horizontalalignment='right',
                 verticalalignment='top',
                 transform=ax1.transAxes,
                 bbox=dict(facecolor='white', alpha=0.8))

    # Title
    plt.title(plot_title)
    fig.tight_layout()
    plt.show()

if __name__ == '__main__':
    plot_optimization_progress()
