"""
This script plots the printout file that contains the convergence history
of the Jacobi-Davidson eigenvalue solver
"""
import enum
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import re


class JDHistoryPlotter:
    def __init__(self, filename):

        # Name string for stdout printout file
        self.filename = filename

        # List of info for each optimization iteration
        self.history = []

        # Regex patterns
        self.pattern_header = re.compile(
            r"Iter     JD Residual      Ritz value      toler"
        )

        # This matches the following formatted string:
        # "4d %15.5e %15.5e %10.2e"
        self.pattern_content = re.compile(
            r"[\s\d]{4}\s{4}[-\s]\d\.[\d]{5}e[+-]\d\d\s{4}"
            r"[-\s]\d\.[\d]{5}e[+-]\d\d\s{2}[\s-]\d\.[\d]{2}e[+-]\d\d"
        )

        return

    def read_file(self):
        """
        Get JD residual, Ritz value and tolerance
        """
        with open(self.filename, "r") as fp:

            entry = None

            for line in fp:
                is_header = self.pattern_header.search(line)
                is_content = self.pattern_content.search(line)

                # If header line, append previous JD history to self.history
                # and reset entry
                if is_header:
                    if entry is not None:
                        self.history.append(entry)
                    entry = []

                # If content line, append it to entry
                if is_content:
                    itr, JDr, ritz, toler = line.split()
                    entry.append(
                        {
                            "itr": int(itr),
                            "JDr": float(JDr),
                            "ritz": float(ritz),
                            "toler": float(toler),
                        }
                    )

            # Append last entry to self.history
            self.history.append(entry)
        return

    def plot(self, plot_type="all", savefig=False, savename="JD_history.png", dpi=800):
        """
        Plot the JD residual, Ritz Value and tolerance

        Args:
            plot_type: 'all' or 'last', plotting entire history or
                       only the history for the last optimization step
        """
        # Set up plotting environment
        # try:
        #     mpl_style_path = os.path.dirname(os.path.realpath(__file__)) \
        #         + '/paper.mplstyle'
        #     plt.style.use(mpl_style_path)
        # except:
        #     print("[Warning] cannot load matplotlib style: paper.mplstyleFalse")
        fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8))
        ax2 = ax.twinx()

        if plot_type == "last":
            entry = self.history[-1]
            itr = [d["itr"] for d in entry]
            JDr = [d["JDr"] for d in entry]
            ritz = [d["ritz"] for d in entry]
            toler = [d["toler"] for d in entry]

        elif plot_type == "all":
            offset = 0
            vl_locs = []
            itr = []
            JDr = []
            ritz = []
            toler = []
            for entry in self.history:
                _itr = [d["itr"] + offset for d in entry]
                _JDr = [d["JDr"] for d in entry]
                _ritz = [d["ritz"] for d in entry]
                _toler = [d["toler"] for d in entry]

                itr.extend(_itr)
                JDr.extend(_JDr)
                ritz.extend(_ritz)
                toler.extend(_toler)

                offset = _itr[-1]
                vl_locs.append(offset)

        itr_converged = [itr[i] for i in range(len(itr)) if JDr[i] < toler[i]]
        JDr_converged = [JDr[i] for i in range(len(itr)) if JDr[i] < toler[i]]
        ritz_converged = [ritz[i] for i in range(len(itr)) if JDr[i] < toler[i]]

        # Plot vertical lines
        if plot_type == "all":
            iL = iR = 0
            xtic = []
            for vl_loc in vl_locs:
                ax.axvline(x=vl_loc, linestyle="-", color="gray", lw=0.5)
                iL = iR
                iR = vl_loc
                xtic.append(int(iL + iR) / 2)

            ax.set_xticks(xtic)
            labels = [i + 1 for i, _ in enumerate(xtic)]
            ax.set_xticklabels(labels)
            ax.tick_params(axis="x", length=0)

        l1 = ax.semilogy(itr, JDr, ".", color="cornflowerblue", label="JD residual")
        l2 = ax.semilogy(
            itr_converged,
            JDr_converged,
            ".",
            color="red",
            label="JD residual - converged",
        )
        l3 = ax.semilogy(itr, toler, lw=1.0, color="orange", label="toler")

        l4 = ax2.plot(itr, ritz, "+", color="gray", label="Ritz value")
        l5 = ax2.plot(
            itr_converged,
            ritz_converged,
            "+",
            color="green",
            label="Ritz value - converged",
        )

        ax2.axhline(y=0.0, linestyle="-", color="gray", lw=1.0)

        ax.set_xlabel("JD iteration")
        ax.set_ylabel("Residual")
        ax2.set_ylabel("Ritz value")
        ax2.set_ylim(np.min(ritz_converged), np.max(ritz_converged))

        lns = l1 + l2 + l3 + l4 + l5
        labs = [l.get_label() for l in lns]
        ax.legend(lns, labs, loc="upper right")

        if not savefig:
            plt.show()
            return

        else:
            fig.savefig(savename, dpi=dpi)
            return


def plot_JD_history(
    filename, plot_type, savefig=True, savename="JD_history.png", dpi=800
):
    plotter = JDHistoryPlotter(filename)
    plotter.read_file()
    plotter.plot(plot_type=plot_type, savefig=savefig, savename=savename, dpi=dpi)
    return


if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument("filename", type=str)
    p.add_argument("--plot-type", type=str, default="all", choices=["all", "last"])
    p.add_argument("--savefig", action="store_true")
    p.add_argument("--savename", type=str, default="JD_history.png")
    p.add_argument("--dpi", type=int, default=800)
    args = p.parse_args()

    plot_JD_history(
        args.filename,
        args.plot_type,
        savefig=args.savefig,
        savename=args.savename,
        dpi=args.dpi,
    )
