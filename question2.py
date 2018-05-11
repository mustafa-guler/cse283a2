import sys
import argparse
import tabulate
from numpy import dot
import matplotlib.pyplot as plt

residue_masses = { "G": 57.02147, 
        "A": 71.03712,
        "S": 87.03203,
        "P": 97.05277,
        "V": 99.06842,
        "T": 101.04768,
        "C": 103.00919,
        "I": 113.08407,
        "L": 113.08407,
        "N": 114.04293,
        "D": 115.02695,
        "Q": 128.05858,
        "K": 128.09497,
        "E": 129.04260,
        "M": 131.04049,
        "H": 137.05891,
        "F": 147.06842,
        "R": 156.10112,
        "Y": 163.06333,
        "W": 186.07932
}

H = 1.01
O = 16
C = 12.01
N = 14.01
missing_H2O = 2*H + O




def calc_theoretical_mass(peptide):
    return sum([residue_masses[aa] for aa in peptide]) + missing_H2O

def ppm_error(pre_mass, can_mass):
    return (abs(pre_mass - can_mass) / pre_mass) * 1e6

def non_ppm_error(pre_mass, can_mass):
    return (abs(pre_mass - can_mass))

def calc_res_mass(peptide):
    return sum([residue_masses[x] for x in peptide])

def q2b(args):
    print(calc_theoretical_mass(args.peptide))

def q2c(args):
    precursor = args.precursor_mass
    with open(args.fasta) as f:
        f.readline() #header
        seq = "".join(map(lambda x:x.strip(), f.readlines()))

    #tryptophan has largest residue mass
    #finds minimum required length of peptide
    min_len = int((precursor - 18) / residue_masses["W"]) 
    seq_len = len(seq)
    results = []

    for i in range(seq_len - min_len + 1):
        curr = calc_theoretical_mass(seq[i:i+min_len])
        err = ppm_error(precursor, curr)
        prev_err = float('inf')
        j = i
        while j < seq_len-1 and (err <= prev_err or err <= 50):
            if err <= 50:
                results.append([i,j+min_len,
                    calc_theoretical_mass(seq[i:j+min_len]),
                    err,seq[i:j+min_len]])

            j += 1
            prev_err = err
            curr = calc_theoretical_mass(seq[i:j+min_len])
            err = ppm_error(precursor, curr)

    print(tabulate.tabulate(results, 
        headers=["start index", "end index", "mass", "ppm error", "peptide"]))

def q2d(args):
    candidates = args.candidates
    for candidate in candidates:
        spectrum = []
        for i in range(1,len(candidate)):
            first = calc_res_mass(candidate[:i])
            second = calc_res_mass(candidate[i:])

            # different ions 
            ions = {}
            ions["a"] = first - C - O + H
            ions["b"] = first + H
            ions["c"] = first + N + 2*H
            ions["x"] = second + C + 2*O + H
            ions["y"] = second + O + H
            ions["z"] = second - N - H + O + H

            # each possible charge for an ion
            for ion in ions:
                for charge in range(1,7):
                    if ion in {"a","b","c"}:
                        seq = candidate[:i] 
                    else:
                        seq = candidate[i:]
                    spectrum.append(((ions[ion]+charge)/charge,ion,charge,seq))

        spectrum.sort()
        with open("{}_spectrum.csv".format(candidate), "w") as f:
            f.write("m/z,ion_type\n")
            for x in spectrum:
                f.write("{:.7f},{}\n".format(x[0],x[1]))

        print("CSV with m/z and ion types for {} written to {}_spectrum.csv".
                format(candidate, candidate))

def q2e(args):
    with open(args.experimental) as f:
        while not f.readline().startswith("C"):
            pass
        exp_mz = []
        S = []
        for line in f:
            if line.startswith("END"):
                break
            mz,intensity = map(float,line.strip().split())
            exp_mz.append(mz)
            S.append(intensity)

    for theoretical in args.theoreticals:
        with open(theoretical) as f:
            f.readline() # header
            theo_mz = []
            for line in f:
                theo_mz.append(float(line.strip().split(",")[0]))

        B = []
        for exp in exp_mz:
            for theo in theo_mz:
                curr = 0
                if non_ppm_error(exp, theo) <= args.E:
                    curr = 1
                    break
            B.append(curr)
            
        print(theoretical.split('_')[0], dot(B,S)/sum(S))
        with open("{}.mz".format(theoretical), "w") as f:
            [f.write("{}\n".format(exp_mz[x])) for x in range(len(exp_mz)) if B[x]]

def q2f(args):
    with open(args.experimental) as f:
        while not f.readline().startswith("C"):
            pass
        exp_mz = []
        exp_intense = []
        for line in f:
            if line.startswith("END"):
                break
            mz,intensity = map(float,line.strip().split())
            exp_mz.append(mz)
            exp_intense.append(intensity)

    plt.figure(figsize=(15,15))
    for theoretical in args.theoreticals:
        bars = plt.bar(exp_mz, exp_intense)
        plt.xlabel("m/z")
        plt.ylabel("intensity")
        with open(theoretical) as f:
            good = set(map(float,f.readlines()))

        make_red = [x in good for x in exp_mz]

        for x,y in zip(bars,make_red):
            if y:
                x.set_color('r')

        plt.savefig("{}.png".format(theoretical.split('.')[0]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="all scripts for quesiton 2 in assignment 3")
    subs = parser.add_subparsers(help="question part to run")

    b = subs.add_parser("2b", help="calculates the theoretical mass of a peptide")
    b.add_argument("peptide", type=str, help="target peptide")
    b.set_defaults(func=q2b)

    c = subs.add_parser("2c", help="generates candidate peptides from a FASTA file, 50ppm error tolerance")
    c.add_argument("fasta", help="path to fasta file to search from")
    c.add_argument("precursor_mass", help="molecular mass of precursor", type=float)
    c.set_defaults(func=q2c)

    d = subs.add_parser("2d", help="generates theoretical spectra for candidate peptides, charges from +1 to +6")
    d.add_argument("candidates", nargs="+", type=str, help="each peptide sequence to generate to theoretical spectrum for")
    d.set_defaults(func=q2d)

    e = subs.add_parser("2e", help="generates spectral similarity scores according to algorithm described")
    e.add_argument("experimental", help="path to mgf file containing experimental spectrum")
    e.add_argument("theoreticals", nargs="+", help="paths to theoretical spectra generated from 2d")
    e.add_argument("-E", type=float, help="error tolerance in ppm, DEFAULT=0.5", default=0.3)
    e.set_defaults(func=q2e)

    f = subs.add_parser("2f", help="plots experimental spectrum")
    f.add_argument("experimental", help="path to mgf file containing experimental spectrum")
    f.add_argument("theoreticals", nargs="+", help="paths to .mz files from 2e")
    f.set_defaults(func=q2f)

    args = parser.parse_args()
    args.func(args)
