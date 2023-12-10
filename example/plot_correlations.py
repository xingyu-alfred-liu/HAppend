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
    plot_agent.plot_correlation(name_list=["bse_Es", "fundamental_gap"], savefig_path="bse_es_fundamental_gap.png", color="b")
    plot_agent = PlotAgent(root_path=root_path)
    plot_agent.plot_correlation(name_list=["bse_Et", "fundamental_gap"], savefig_path="bse_et_fundamental_gap.png", color="r")
    plot_agent = PlotAgent(root_path=root_path)
    plot_agent.plot_correlation(name_list=["fundamental_gap", "bse_Es_bind"], savefig_path="bse_es_bse_es_bind.png", color="g")
    plot_agent = PlotAgent(root_path=root_path)
    plot_agent.plot_correlation(name_list=["fundamental_gap", "bse_Et_bind"], savefig_path="bse_es_bse_et_bind.png", color="m")
    plot_agent = PlotAgent(root_path=root_path)
    plot_agent.plot_correlation(name_list=["fundamental_gap", "bandgap"], savefig_path="fundamental_band_gap.png", color="r")
