#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/gff_io.h>
#include <math.h>       /* sqrt */
#include <seqan/store.h> /* FragmentStore */
#include <queue>
#include <vector>
#include <ctime>
#include "boost/multi_array.hpp"
#include <cassert>
#include <boost/unordered_map.hpp>
#include <string>
#include <thread>
#include <mutex>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <iomanip>

#include <sstream>
#include <fstream>
#include <string>

using namespace seqan;
using namespace std;

struct ModifyStringOptions
{
   CharString inputFileName;
   CharString outputFileName;
   CharString outputReadsFileName;
   unsigned int percentage;
   unsigned int reads;
   unsigned int readLength;
   bool methylate = false;
   bool nonDirectional = false;
};

ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, 
                                             int argc, char const ** argv)
{
   ArgumentParser parser("methylation_simulator");
   addOption(parser, ArgParseOption("i", "input-file", 
                                    "Path to the input file", 
                                    ArgParseArgument::INPUT_FILE, "IN"));
   setRequired(parser, "input-file");
   addOption(parser, ArgParseOption("o", "output-file", 
                                    "Path to the output file", 
                                    ArgParseArgument::OUTPUT_FILE, "OUT"));
   setRequired(parser, "output-file");
   addOption(parser, ArgParseOption("r", "reads-output",
                                    "Path to the reads output file",
                                    ArgParseArgument::OUTPUT_FILE, "OUT"));
   setRequired(parser, "reads-output");
   addOption(parser, ArgParseOption("p", "percentage", 
                                    "Percent of C's methylated",
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "percentage", "70");
   addOption(parser, ArgParseOption("n", "reads",
                                    "Number of reads to produce, million.",
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "reads", "30");
   addOption(parser, ArgParseOption("l", "read-length",
                                    "Simulated read length",
                                    ArgParseArgument::INTEGER, "INT"));
   setDefaultValue(parser, "read-length", "75");
   addOption(parser, ArgParseOption("DN", "no-methyl", 
                                    "Do not methylate the genome."));
   addOption(parser, ArgParseOption("nd", "non-directional", 
                                    "Create non-directional reads."));

   setShortDescription(parser, "BS Mapper");
   setVersion(parser, "0.0.1");
   setDate(parser, "March 2018");
   addUsageLine(parser, "-i sequence.fa -o output.fa [\\fIOPTIONS\\fP] ");

   addDescription(parser, "Create a methylated genome and reads.");
   ArgumentParser::ParseResult res = parse(parser, argc, argv);

   // If parsing was not successful then exit with code 1 if there were errors.
   // Otherwise, exit with code 0 (e.g. help was printed).
   if(res != seqan::ArgumentParser::PARSE_OK)
      return res;

   getOptionValue(options.inputFileName, parser, "input-file");
   getOptionValue(options.outputFileName, parser, "output-file");
   getOptionValue(options.outputReadsFileName, parser, "reads-output");
   getOptionValue(options.percentage, parser, "percentage");
   getOptionValue(options.reads, parser, "reads");
   getOptionValue(options.readLength, parser, "read-length");
   options.methylate = isSet(parser, "no-methyl");
   options.nonDirectional = isSet(parser, "non-directional");

   return seqan::ArgumentParser::PARSE_OK;
}

/*
   What do we need to achieve
   THe basics
      Produce a methylated genome file
		
      Make reads from this file
         * this could just be wgsim but we need to account for 
           partial methylation

*/
int main(int argc, char const ** argv)
{
   //parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options, 
                                                      argc, argv);

   if(res != ArgumentParser::PARSE_OK)
      return res == ArgumentParser::PARSE_ERROR;

   // Store the current contig
   CharString id;
   IupacString seq;

   // Open our input file
   SeqFileIn seqFileIn(toCString(options.inputFileName));

   string CGgffout = toCString(options.outputFileName) + string(".CG.w1.gff");
   GffFileOut CGgffOut(toCString(CGgffout));
   string CHGgffout = toCString(options.outputFileName) + string(".CHG.w1.gff");
   GffFileOut CHGgffOut(toCString(CHGgffout));
   string CHHgffout = toCString(options.outputFileName) + string(".CHH.w1.gff");
   GffFileOut CHHgffOut(toCString(CHHgffout));

   // Open the output file
   SeqFileOut seqFileOut(toCString(options.outputFileName));

   // methylate the genome
   if(options.methylate == false)
   {
      while(!atEnd(seqFileIn))
      {
         try
         {
            readRecord(id, seq, seqFileIn);

            // Loop through contigs
            for(unsigned int i = 1; i < length(seq); i++)
            {
               if(seq[i-1] == 'C')
               {
                  int num = rand() % 100 + 1;
                  if(num > options.percentage)
                     seq[i-1] = 'T';

                  GffRecord record;
                  record.ref = id;
                  record.beginPos = i - 1;
                  record.endPos = i;
                  record.source = "xftools";
                  if(seq[i-1] == 'C')
                     record.score = 1;
                  else
                     record.score = 0;

                  if(seq[i] == 'G')
                  {
                     record.type = "CG";
                     writeRecord(CGgffOut, record);
                  }
                  else if(seq[i+1] == 'G')
                  {
                     record.type = "CHG";
                     writeRecord(CHGgffOut, record);
                  }
                  else
                  {
                     record.type = "CHH";
                     writeRecord(CHHgffOut, record);
                  }
               }

               if(seq[i-1] != 'A' && seq[i-1] != 'C' && seq[i-1] != 'G' && seq[i-1] != 'T')
                  seq[i-1] = 'N';
            
            }
            writeRecord(seqFileOut, id, seq);
         }
         catch (Exception const & e)
         {
            std::cout << "ERROR: " << e.what() << std::endl;
            return 1;
         }
      }
   }
   else
   {
      cout << "Not methylating the genome" << endl;
   }
   close(seqFileIn);
   close(seqFileOut);

   // read in the fasta file again and then extract reads
   SeqFileIn seqFileRefIn;

   if(options.methylate == true)
   {
      open(seqFileRefIn, toCString(options.inputFileName));
      cout << "Opening " << options.outputFileName << endl;
   }
   else
   {
      open(seqFileRefIn, toCString(options.outputFileName));
      cout << "Opening " << options.outputFileName << endl;
   }

   cout << "Finished opening genome" << endl;

   SeqFileOut seqFileReadsOut(toCString(options.outputReadsFileName));
   vector<IupacString> meh; // store the genome in MEM

   while(!atEnd(seqFileRefIn))
   {
      try
      {
         readRecord(id, seq, seqFileRefIn);
         meh.push_back(seq);
         cout << "Pushed back " << id << endl;
      }
      catch (Exception const & e)
      {
         std::cout << "ERROR: " << e.what() << std::endl;
         return 1;
      }
   }

   close(seqFileRefIn);

   srand(time(NULL));

   // now let's make reads
   for(unsigned int i = 0; i < (options.reads*1000000); i++)
   {
      unsigned int contig = rand() % meh.size();
      unsigned int pos = rand() % length(meh[contig]);
      CharString id = "simdata_" + to_string(i+1) + "_" + to_string(contig+1) + "_" + to_string(pos+1); // all +1 because that is the way that we would see it in SAM tools.
      Dna5String read;
      read = infix(meh[contig], pos, pos+options.readLength);

      unsigned int flip = rand() % 2;

      // Randomly reverse compliment 1/2 the reads
      // if non-directional
      if(options.nonDirectional == true)
      {
         if(flip == 1)
         {
            reverseComplement(read);
         }
      }

      if(pos+options.readLength < length(meh[contig]))
      {
         writeRecord(seqFileReadsOut, id, read);
      }
   }
 
   close(seqFileReadsOut);
   return 0;
}
