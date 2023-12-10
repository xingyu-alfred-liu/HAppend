import argparse

from HAppend.plot import PlotAgent

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--rootpath", type=str, help="root path point to all JSON")

    args = parser.parse_args()
    if args.rootpath:
        root_path = args.rootpath
        print("root path: %s" % root_path)
    
    plot_agent = PlotAgent(root_path=root_path)
    plot_agent.plot_histogram(name="bse_Es", savefig_path="singlet_optical_gap")
    plot_agent = PlotAgent(root_path=root_path)
    plot_agent.plot_histogram(name="fundamental_gap", savefig_path="fundamental_gap")
    plot_agent = PlotAgent(root_path=root_path)
    plot_agent.plot_histogram(name="bse_Es_bind", savefig_path="bse_Es_bind")
    plot_agent = PlotAgent(root_path=root_path)
    plot_agent.plot_histogram(name="bse_Et_bind", savefig_path="bse_Et_bind")
