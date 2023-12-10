import argparse

from HAppend.plot import PlotAgent

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--rootpath", type=str, help="root path point to all JSON")
    parser.add_argument("-id", "--struct_id", type=str, help="struct_id match with materials name")
    parser.add_argument("-b", "--bandstructure", type=bool, help="if plot bandstructure")
    parser.add_argument("-a", "--absorption", type=bool, help="if plot absorption spectrum")
    parser.add_argument("-bp", "--bandstructure_path", type=str, help="target path to save bandstructure plot")
    parser.add_argument("-ap", "--absorption_path", type=str, help="target path to save absorption spectrum")

    args = parser.parse_args()
    if args.rootpath:
        root_path = args.rootpath
        print("root path: %s" % root_path)
    if args.struct_id:
        struct_id = args.struct_id
        print("struct id: %s" % struct_id)
    if args.bandstructure:
        plot_band = True
    else:
        plot_band = False
    if args.absorption:
        plot_absorp = True
    else:
        plot_absorp = False
    if args.bandstructure_path:
        bandstructure_path = args.bandstructure_path
    if args.absorption_path:
        absorption_path = args.absorption_path
    
    # root_path = "/Users/alfred/Desktop/work/PAH101/PAH101Plot/data"
    plot_agent = PlotAgent(root_path=root_path)
    if plot_band:
        print("Plot Bandstructure...")
        bandstruct_fig = plot_agent.plot_bandstructure(
            struct_id=struct_id, 
            savefig_path=bandstructure_path,
            band_type="MF",
        )
    if plot_absorp:
        print("Plot Absorption Spectrum...")
        absorption_fig = plot_agent.plot_absorption(
            struct_id=struct_id,
            savefig_path=absorption_path
        )
