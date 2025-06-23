

import sys, os, shutil, glob, argparse
# Import performance_analysis module
file_loc = os.path.dirname(__file__)
sys.path.append(os.path.join(file_loc, "../devtools"))
import performance_analysis as pf
import thicket as th
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
from matplotlib.legend_handler import HandlerTuple

parser = argparse.ArgumentParser()

parser.add_argument("--out-dir", default="./time_doc")
parser.add_argument("--bench", default="/usr/workspace/sduser/Spheral/benchmarks")
args = parser.parse_args()

bench = args.bench
out_dir = args.out_dir

# TODO: This is a very hacky way to handle this
config_shorthand = {"toss_4_x86_64_ib_clang_14.0.6_mvapich2_2.3.6": "clang",
                    "toss_4_x86_64_ib_gcc_10.3.1_mvapich2_2.3.6": "gcc",
                    "toss_4_x86_64_ib_cray_rocmcc_6.2.0_cray-mpich_8.1.31": "~rocm",
                    "toss_4_x86_64_ib_cray_rocmcc_6.2.0_cray-mpich_8.1.31_rocm": "+rocm"}

metric = "Avg time/rank"
region = "advance"
# Which tests to display
disp_tests = ["3DCONV", "3DTAYLORFSI", "3DTAYLORSOLIDSPH",
              "3DTAYLORCRK", "NS3DCRKSPH", "NS3DSPH", "NS3DFSISPH",
              "NC2DSPH"]
def plot_hist_times(test_name, bench, clusters, region, metric, savefile=None):
    num_clusters = len(clusters)
    figs, axes = plt.subplots(num_clusters, 1, sharex=True)
    if (num_clusters == 1):
        axes = [axes]
    marker_style = dict(alpha=0.5, markersize=5, fillstyle="none")
    for i, cluster in enumerate(clusters):
        lgd_tups = []
        lgd_names = []
        data = pf.get_hist_times(test_name, bench, cluster, region)
        ax = axes[i]
        for key, tk in data.items():
            assert (tk.metadata["cluster"].iloc[0] == cluster)
            date_group = pf.group_dates(tk)
            avgtimes = []
            lotimes = []
            hitimes = []
            dates = []
            alltimes = []
            alldates = []
            for cdate, ctest in date_group.items():
                new_metric = th.stats.mean(ctest, [metric])[0]
                # Region is checked for in get_hist_times, no need to check again
                avgtimes.append(pf.get_times(ctest.statsframe, region, new_metric)[0])
                vals = pf.get_times(ctest, region, metric)
                dates.append(cdate[0])
                lotimes.append(min(vals))
                hitimes.append(max(vals))
                alldates.extend([cdate for x in range(len(vals))])
                alltimes.extend(vals)
            install_config = tk.metadata["install_config"].iloc[0]
            lgd_entry = f"{cluster} {config_shorthand[install_config]}"
            p1, = ax.plot(dates, avgtimes)
            p2 = ax.fill_between(dates, lotimes, hitimes, color=p1.get_color(), alpha=0.2)
            p3, = ax.plot(alldates, alltimes, "o", color=p1.get_color(), **marker_style)
            lgd_tups.append((p1, p2, p3))
            lgd_names.append(lgd_entry)
        ax.xaxis.set_major_formatter(mdate.DateFormatter('%Y-%b'))
        for label in ax.get_xticklabels(which='major'):
            label.set(rotation=30, horizontalalignment='right')
        ax.legend(lgd_tups, lgd_names, fancybox=True,
                  handler_map={tuple: HandlerTuple(ndivide=None)})
        ax.set_ylabel(f"{metric} (s)")
        if (i == 0):
            ax.set_title(f"{test_name}, region: {region}")
    figs.tight_layout(pad=1.1)
    figs.set_figheight(7.5)
    plt.tight_layout()
    if (savefile):
        print(f"Writing {savefile}")
        plt.savefig(savefile)
    else:
        plt.show()

if (os.path.exists(out_dir)):
    shutil.rmtree(out_dir)
os.mkdir(out_dir)
conf_loc = os.path.join(file_loc, "../../docs/conf.py")
shutil.copy(conf_loc, out_dir)
indx_file = """
Historical Timings
##################

"""
# Retrieve cluster names from benchmark directories
bdirs = glob.glob(os.path.join(bench, "**/**"))
cluster_names = set(os.path.basename(x) for x in bdirs)
# Retrieve all test names from benchmark file names
all_files = glob.glob(os.path.join(bench, "**/*.cali"), recursive=True)
test_names = set(os.path.basename(x).split("_")[0] for x in all_files)
plot_images = {}
for test in test_names:
    if (test not in disp_tests):
        continue
    print(f"Getting data for {test}")
    pfile = f"{test}.png"
    plot_file = os.path.join(out_dir, pfile)
    plot_hist_times(test, bench, cluster_names, region, metric, savefile=plot_file)
    plot_images.update({test: pfile})

plot_images = dict(sorted(plot_images.items()))
for test, pfile in plot_images.items():
    indx_file += f"""
.. dropdown:: {test}

"""
    indx_file += f"""
    .. image:: {pfile}

"""

with open(os.path.join(out_dir, "index.rst"), "w") as f:
    f.write(indx_file)
