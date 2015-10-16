//
// Created by Jelmer Mulder on 13/04/15.
//

#include <iostream>
#include "lana.h"
#include "config.h"
#include <lemon/smart_graph.h>
#include <lemon/arg_parser.h>
#include "verbose.h"
#include "output/output.h"
#include "bronkerboschconnected.h"

using namespace lemon;
using namespace nina;
using namespace nina::gna;


typedef SmartGraph Graph;
typedef SmartBpGraph BpGraph;
typedef Lana<Graph, BpGraph> LanaType;
typedef Parser<Graph> ParserType;
typedef BpParser<Graph, BpGraph> BpParserType;
typedef Output<Graph, BpGraph> OutputType;


// SIGUSR1 can be used to indicate that the user wants to check the status of the algorithm.
// When SIGUSR1 is received, the volatile variable sig_caught is set to 1.
// In every invocation of bk_pivot in BronKerboschConnected and BronKerboschConnectedRelaxed
// this variable is checked. If it is set to 1, the histogram of the sizes of the alignments
// found in the current component are printed to the file specified by the '-hist' parameter.
// Note that ONLY the alignments found in the current component are reported. Also note that
// the reported alignements can still contain automorphisms.
// To view the histogram of all alignments across all componenets, without automorphisms, one
// must wait for entire algorithm to terminate.

volatile int sig_caught = 0;

void handle_sigusr1(int signum) {
    if (signum == SIGUSR1) {
        sig_caught = 1;
    }
}

int main(int argc, char** argv)
{
    signal(SIGUSR1, handle_sigusr1);


    // TODO: Update readme instructions for OSX.

    ArgParser ap(argc, (char const *const *) argv);
    Options options;

    std::string g1, g2, gm, outputFile, hist_file_name, prod_file_name, comp_hist_file_name;
    int inputFormatG1 = static_cast<int>(LanaType::IN_STRING);
    int inputFormatG2 = static_cast<int>(LanaType::IN_STRING);
    int inputFormatGm = static_cast<int>(LanaType::BP_IN_BLAST);
    int verbosityLevel = static_cast<int>(VERBOSE_NON_ESSENTIAL);
    int outputType = static_cast<int>(OutputType::ORIG_EDGES);
    double eValCutOff = 10.;
    int nSamples = 0;
    std::string outputFormat = "3";

    ap
            .boolOption("version", "Show version number")
            .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false)
            .synonym("-verbosity", "v")
            .refOption("hist", "File name for the histogram of alignment sizes.", hist_file_name, false)
            .refOption("comp_hist", "File name for the histogram of connected-component sizes.", comp_hist_file_name, false)
            .refOption("prod", "File name for the product graph DOT-file output.", prod_file_name, false)
            .refOption("g1", "File name of input graph G_1", g1, false)
            .refOption("g2", "File name of input graph G_2", g2, false)
            .refOption("gm", "File name in which matching edges of G_m are defined;\n"
                    "     if omitted, the complete graph is used", gm, false)
            .refOption("if1", "Specifies the input file format for G_1:\n"
                    "     0 - GML format\n"
                    "     1 - GraphML format\n"
                    "     2 - STRING format (default)\n"
                    "     3 - LGF format\n"
                    "     4 - CSV format\n"
                    "     5 - LEDA format\n"
                    "     6 - Edge list format", inputFormatG1, false)
            .refOption("if2", "Specifies the input file format for G_2:\n"
                    "     0 - GML format\n"
                    "     1 - GraphML format\n"
                    "     2 - STRING format (default)\n"
                    "     3 - LGF format\n"
                    "     4 - CSV format\n"
                    "     5 - LEDA format\n"
                    "     6 - Edge list format", inputFormatG2, false)
            .refOption("ifm", "Specifies the input file format for G_m:\n"
                    "     0 - Candidate list\n"
                    "     1 - BLAST (default)\n"
                    "     2 - LGF", inputFormatGm, false)
            .refOption("o", "Output file name",
                       outputFile, false)
            .refOption("of", "Specifies the output file format:\n"
                    "     0 - DOT format\n"
                    "     1 - GML format\n"
                    "     2 - LGF format\n"
                    "     3 - SIF format (default)\n"
                    "     4 - JSON format\n"
                    "     5 - NEATO format\n"
                    "     6 - CSV (matched) format\n"
                    "     7 - CSV (unmatched in G_1) format\n"
                    "     8 - CSV (unmatched in G_2) format\n"
                    "     9 - CSV (alignment) format\n"
                    "    10 - EDA format\n"
                    "    11 - NOA format\n"
                    "    12 - ANALYSE (SIF, EDA and NOA) format", outputFormat, false)
            .refOption("op", "Specifies parts of the matching graph to output:\n"
                    "     0 - Nodes and matching edges present in the solution\n"
                    "     1 - Nodes, matching edges and original edges present\n"
                    "         in the solution (default)\n"
                    "     2 - Nodes and matching edges present in the solution\n"
                    "         as well as all original edges", outputType, false)
            .refOption("e", "Specifies e value cutoff (only used in conjunction\n"
                                        "     with -ifm 1 and -ifm 2, default: 10)",
                                                eValCutOff, false)
            .refOption("mcs", "Minimum clique size. Smaller cliques are not processed. (default: 0)",
                       options._minCliqueSize, false)
            .refOption("mse", "Maximum s-edges. The higher this number, the more\n"
                    "     s-edges can be used, and thus the more relaxation is applied.\n"
                    "     (default: 0)",
                       options._nMaxSEdges, false)
            .refOption("p", "Compute empirical p-value using specified number of samples (default: 0)",
                       nSamples, false)
            .boolOption("a", "Do not remove automorphisms.")
            .boolOption("sol", "Print all found alignments in Human-Readable format to STDOUT.\n"
                    "     Should probably only be used in combination with -mcs");

    ap.parse();

    options._removeAutomorphisms = !ap.given("a");
    options._printProductVector = ap.given("sol");
    options._hist_file_name = &hist_file_name;
    options._comp_hist_file_name = &comp_hist_file_name;
    options._prod_file_name = &prod_file_name;

    if (ap.given("version"))
    {
        std::cout << "Version number: " << LANA_VERSION << std::endl;
        return 0;
    }

    if (!ap.given("g1") || !ap.given("g2"))
    {
        std::cerr << "Both -g1 and -g2 need to be specified" << std::endl;
        return 1;
    }


    outputFormat += ",";

    g_verbosity = static_cast<VerbosityLevel>(verbosityLevel);

    LanaType lana(options);

    ParserType* pParserG1 =
            LanaType::createParser(g1,
                                      static_cast<LanaType::InputFormatEnum>(inputFormatG1));

    ParserType* pParserG2 =
            LanaType::createParser(g2,
                                      static_cast<LanaType::InputFormatEnum>(inputFormatG2));

    BpParserType* pParserGm =
            LanaType::createBpParser(gm,
                                        static_cast<LanaType::BpInputFormatEnum>(inputFormatGm),
                                        pParserG1,
                                        pParserG2,
                                        eValCutOff);

    if (!lana.init(pParserG1, pParserG2, pParserGm))
        return 1;

    if (lana.solve() == 0)
    {
//        lana.parseOutputString(outputFormat);
        // TODO: Test output formats other than .sif
//        lana.generateOutput(static_cast<OutputType::OutputType>(outputType), outputFile);
    }


    // Prints the frequencies of the sizes of the found alignments.
    if (!hist_file_name.empty())
    {
        std::ofstream hist_file;
        hist_file.open(hist_file_name.c_str());
        lana.printCliqueSizeFrequencies(hist_file);
        hist_file.close();
    }




}

