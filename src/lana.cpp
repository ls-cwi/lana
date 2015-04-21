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

using namespace lemon;
using namespace nina;
using namespace nina::gna;


typedef SmartGraph Graph;
typedef SmartBpGraph BpGraph;
typedef Lana<Graph, BpGraph> LanaType;
typedef Parser<Graph> ParserType;
typedef BpParser<Graph, BpGraph> BpParserType;
typedef Output<Graph, BpGraph> OutputType;

int main(int argc, char** argv)
{
    ArgParser ap(argc, (char const *const *) argv);

    std::string g1, g2, gm, outputFile;
    int inputFormatG1 = static_cast<int>(LanaType::IN_STRING);
    int inputFormatG2 = static_cast<int>(LanaType::IN_STRING);
    int inputFormatGm = static_cast<int>(LanaType::BP_IN_BLAST);
    int verbosityLevel = static_cast<int>(VERBOSE_NON_ESSENTIAL);
    int outputType = static_cast<int>(OutputType::ORIG_EDGES);
    std::string outputFormat = "3";

    ap
            .boolOption("version", "Show version number")
            .refOption("v", "Specifies the verbosity level:\n"
                    "     0 - No output\n"
                    "     1 - Only necessary output\n"
                    "     2 - More verbose output (default)\n"
                    "     3 - Debug output", verbosityLevel, false)
            .synonym("-verbosity", "v")
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
                    "         as well as all original edges", outputType, false);
    ap.parse();

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

    LanaType lana;


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
                                        pParserG2);

    if (!lana.init(pParserG1, pParserG2, pParserGm))
        return 1;

    if (lana.solve() == 0)
    {
        lana.parseOutputString(outputFormat);
        lana.generateOutput(static_cast<OutputType::OutputType>(outputType), outputFile);
    }

}

