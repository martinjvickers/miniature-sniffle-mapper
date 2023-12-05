#ifndef COMMANDLINEPARSE_H
#define COMMANDLINEPARSE_H

#include <iostream>
#include <seqan/sequence.h>
#include <string>

using namespace seqan;

struct ModifyStringOptions
{
   CharString inputFileName = NULL;
   CharString referenceFileName;
   CharString inputLeftFileName = NULL;
   CharString inputRightFileName = NULL;
   CharString outputBamFileName;
   CharString ambigOutputBamFileName = NULL;
   int cpu_cores = 1;
   bool stranded = true;
   bool debug = false;
   bool keep_unmapped = false;
   bool non_directional = false;
   bool singleEnd;
   bool compress = false;
   //bool outputambig = false;
   CharString ambig;
   std::string command;
};

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options,
                                             int argc, char const ** argv)
{
   ArgumentParser parser("miniature-sniffle-mapper");
   addOption(parser, ArgParseOption("i", "input-file", "Input filename",
                                    ArgParseArgument::INPUT_FILE, "IN"));
   addOption(parser, ArgParseOption("o", "output-file", "Output filename",
                                    ArgParseArgument::OUTPUT_FILE, "OUT"));
   setRequired(parser, "output-file");
   addOption(parser, ArgParseOption("r", "reference-file", "Reference file",
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "reference-file");
   addOption(parser, ArgParseOption("af", "ambig-file", "Ambigious reads file. \
                                    This is the output file where you would put\
                                    all the ambigious reads if not in the main \
                                    output file.",
                                    ArgParseArgument::OUTPUT_FILE, "OUT"));
   addOption(parser, ArgParseOption("d", "debug", "Print out debug information.\
                                    There is a lot of text, you may wish to \
                                    redirect it to a file."));
   addOption(parser, ArgParseOption("k", "keep-unmapped", "Keep unmapped reads.\
                                    If ambig-file is set then they will go \
                                    there, else they will go into the output\
                                    file."));
   addOption(parser, ArgParseOption("nd", "non-directional",
                                    "non-directional"));
   addOption(parser, ArgParseOption("gz", "compress",
                                    "Compress temporary files, this will \
                                    increase the run time"));
   addOption(parser, ArgParseOption("1", "input-file-left",
                                    "Path to the left reads",
                                    ArgParseArgument::INPUT_FILE, "IN"));
   addOption(parser, ArgParseOption("2", "input-file-right",
                                    "Path to the right reads",
                                    ArgParseArgument::INPUT_FILE, "IN"));
   addOption(parser, ArgParseOption("c", "num-cores", "Number of Bowtie CPU \
                                    cores to use. NOTE: This is per bowtie2 \
                                    instance (x2 for directional and x4 for \
                                    non-directional).",
                                    ArgParseArgument::INTEGER, "INT"));
   addOption(parser, ArgParseOption("a", "ambig", "Define the way in which we \
                                    will deal with ambiguous reads.",
                                    ArgParseArgument::STRING, "STR"));
   setValidValues(parser, "ambig", "ignore random all");
   setDefaultValue(parser, "ambig", "ignore");

   setDefaultValue(parser, "num-cores", "1");
   setShortDescription(parser, "BS Mapper");
   setVersion(parser, "0.0.6");
   setDate(parser, "November 2018");
   addUsageLine(parser, "-i sequence.fastq [\\fIOPTIONS\\fP] ");

   addDescription(parser, "Perform the mapping.");
   ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

   if(res != seqan::ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.outputBamFileName, parser, "output-file");
   getOptionValue(options.ambigOutputBamFileName, parser, "ambig-file");
   getOptionValue(options.inputLeftFileName, parser, "input-file-left");
   getOptionValue(options.inputRightFileName, parser, "input-file-right");
   getOptionValue(options.referenceFileName, parser, "reference-file");
   options.non_directional = isSet(parser, "non-directional");
   options.debug = isSet(parser, "debug");
   options.keep_unmapped = isSet(parser, "keep-unmapped");
   getOptionValue(options.ambig, parser, "ambig");
   // if it's single end, then true
   options.singleEnd = isSet(parser, "input-file");
   getOptionValue(options.cpu_cores, parser, "num-cores");
   options.compress = isSet(parser, "compress");

   return seqan::ArgumentParser::PARSE_OK;
}

#endif
