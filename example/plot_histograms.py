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
