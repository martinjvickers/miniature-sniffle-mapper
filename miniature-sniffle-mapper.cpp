#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/store.h>
#include <string>
#include <seqan/bam_io.h>
#include <map>
#include <vector>
#include "tbb/tbb.h"
#include "tbb/atomic.h"
#include <stdlib.h>
#include "commandLineParse.h"
#include "prepareReads.h"
#include "calcMethylTag.h"
#include "bowtieCommands.h"
#include "processAlignments.h"

using namespace seqan;
using namespace std;
using namespace tbb;

int mapping(ModifyStringOptions options, BamFileOut &newBamFileOut,
            BamFileOut &ambigBamFileOut, FaiIndex &faiIndex)
{
   FILE *in_GA_GA, *in_CT_GA, *in_GA_CT, *in_CT_CT;

   // these are the BamIOContexts required for each stream
   BamHeader header_GA_GA, header_GA_CT, header_CT_GA, header_CT_CT;
   StringSet<CharString> referenceNameStore_GA_GA;
   NameStoreCache<StringSet<CharString>> referenceNameStoreCache_GA_GA(referenceNameStore_GA_GA);
   BamIOContext<StringSet<CharString>> bamIOContext_GA_GA(referenceNameStore_GA_GA, referenceNameStoreCache_GA_GA);
   StringSet<CharString> referenceNameStore_GA_CT;
   NameStoreCache<StringSet<CharString>> referenceNameStoreCache_GA_CT(referenceNameStore_GA_CT);
   BamIOContext<StringSet<CharString>> bamIOContext_GA_CT(referenceNameStore_GA_CT, referenceNameStoreCache_GA_CT);
   StringSet<CharString> referenceNameStore_CT_GA;
   NameStoreCache<StringSet<CharString>> referenceNameStoreCache_CT_GA(referenceNameStore_CT_GA);
   BamIOContext<StringSet<CharString>> bamIOContext_CT_GA(referenceNameStore_CT_GA, referenceNameStoreCache_CT_GA);
   StringSet<CharString> referenceNameStore_CT_CT;
   NameStoreCache<StringSet<CharString>> referenceNameStoreCache_CT_CT(referenceNameStore_CT_CT);
   BamIOContext<StringSet<CharString>> bamIOContext_CT_CT(referenceNameStore_CT_CT, referenceNameStoreCache_CT_CT);

   string bowtie_GA_GA, bowtie_GA_CT, bowtie_CT_GA, bowtie_CT_CT;

   SeqFileIn seqFileIn;
   SeqFileIn seqFileInR1;
   SeqFileIn seqFileInR2;

   /*The original data*/
   if(options.inputFileName != NULL)
   {
      open(seqFileIn, toCString(options.inputFileName));
   }
   else if(options.inputLeftFileName != NULL && 
           options.inputRightFileName != NULL &&
           options.inputFileName == NULL)
   {
      open(seqFileInR1, toCString(options.inputLeftFileName));
      open(seqFileInR2, toCString(options.inputRightFileName));
   }

   // get the bowtie commands to run
   prepareBowtie2Commands(options, bowtie_GA_GA, bowtie_GA_CT, 
                          bowtie_CT_GA, bowtie_CT_CT);

   // run bowtie!
   if(options.non_directional == true)
   {
      if(!(in_GA_GA = popen(bowtie_GA_GA.c_str(), "r")))
      {
         cerr << "Error:\t" << in_GA_GA << "\tfailed to start." << endl;
         return 1;
      }
      if(!(in_GA_CT = popen(bowtie_GA_CT.c_str(), "r")))
      {
         cerr << "Error:\t" << in_GA_CT << "\tfailed to start." << endl;
         return 1;
      }
   }

   if(!(in_CT_GA = popen(bowtie_CT_GA.c_str(), "r")))
   {
      cerr << "Error:\t" << in_CT_GA << "\tfailed to start." << endl;
      return 1;
   }
   if(!(in_CT_CT = popen(bowtie_CT_CT.c_str(), "r")))
   {
      cerr << "Error:\t" << in_CT_CT << "\tfailed to start." << endl;
      return 1;
   }

   // note, this is a globally accessible output file
   newBamFileOut.context = bamIOContext_CT_CT;

   if(options.ambigOutputBamFileName != NULL)
      ambigBamFileOut.context = bamIOContext_CT_CT;

   parse_buffer(in_GA_GA, in_CT_GA, in_GA_CT, in_CT_CT, options, header_GA_GA, 
                header_GA_CT, header_CT_GA, header_CT_CT, bamIOContext_GA_GA, 
                bamIOContext_CT_GA, bamIOContext_GA_CT, bamIOContext_CT_CT, 
                ambigBamFileOut, newBamFileOut, faiIndex, seqFileIn, 
                seqFileInR1, seqFileInR2);

   // close buffers and open files
   if(options.non_directional == true)
      pclose(in_GA_CT), pclose(in_GA_GA);

   pclose(in_CT_GA), pclose(in_CT_CT);
   close(newBamFileOut);

   if(options.ambigOutputBamFileName != NULL)
      close(ambigBamFileOut);

   return 0;
}

int openFiles(ModifyStringOptions options, BamFileOut &newBamFileOut, 
              BamFileOut &ambigBamFileOut, FaiIndex &faiIndex)
{
   // open output BAM file
   if(!open(newBamFileOut, toCString(options.outputBamFileName)))
   {
      cerr << "ERROR: Could not open the output file: ";
      cerr << options.outputBamFileName << endl;
      return 1;
   }

   // open output ambig output file if needed
   if(options.ambigOutputBamFileName != NULL)
   {
      if(!open(ambigBamFileOut, toCString(options.ambigOutputBamFileName)))
      {
         cerr << "ERROR: Could not open the ambig output file: ";
         cerr << options.ambigOutputBamFileName << endl;
         return 1;
      }
   }

   // index the reference file
   if(!build(faiIndex, toCString(options.referenceFileName)))
   {
      cerr << "ERROR: Could not build FAI index for file ";
      cerr << options.referenceFileName << ".\n";
      return 1;
   }

   return 0;
}

int main(int argc, char const ** argv)
{
   // parse our options
   ModifyStringOptions options;
   ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
   if(res != seqan::ArgumentParser::PARSE_OK)
      return res == seqan::ArgumentParser::PARSE_ERROR;

   // store the user input for the BAM file header
   string command;
   for(int m = 0; m < argc; m++)
   {
      command.append(argv[m]);
      command.append(" ");
   }
   options.command = command;

   // open necessary files
   BamFileOut newBamFileOut, ambigBamFileOut;
   FaiIndex faiIndex;
   // TODO: catch the return code here too
   openFiles(options, newBamFileOut, ambigBamFileOut, faiIndex);

   // process reads
   if(options.singleEnd == true)
   {
      cout << "Single End Reads" << endl;
      if(process_reads(options, options.inputFileName) != 0)
      {
         cerr << "ERROR: Reads could not be processed" << endl;
         return 1;
      }
   }
   else
   {
      cout << "Paired-end Reads" << endl;
      if(process_reads(options, options.inputLeftFileName) != 0 || 
         process_reads(options, options.inputRightFileName) != 0)
      {
         cerr << "ERROR: Reads could not be processed" << endl;
         return 1;
      }
   }

   // map the reads
   // TODO: there is something odd with this return. !mapping(options)
   // should work but it doesn't. mjv08
   if(mapping(options, newBamFileOut, ambigBamFileOut, faiIndex))
   {
      cerr << "ERROR: Failed to map" << endl;
      return 1;
   }

   // clean up
   if(options.singleEnd == true)
   {
      cout << "Cleaning up temp files" << endl;
      if(cleanup(options, options.inputFileName) != 0)
      {
         cerr << "ERROR: Temp files were not deleted" << endl;
         return 1;
      }
   }

   // print out the stats
   cout << "Number of Reads:	" << runstats.reads << endl;
   cout << "Unique Mapped:		" << runstats.uniq_mapped << "(";
   cout << ((float)runstats.uniq_mapped / (float)runstats.reads) * 100 << "%)";
   cout << endl;
   cout << "Ambigious Reads:	" << runstats.ambig << "(";
   cout << ((float)runstats.ambig / (float)runstats.reads) * 100 << "%)";
   cout << endl;
   cout << "Did not map:		";
   cout << runstats.reads - (runstats.uniq_mapped + runstats.ambig) << "(";
   cout << ((float)(runstats.reads - (runstats.uniq_mapped + runstats.ambig)) / (float)runstats.reads) * 100;
   cout << "%)"<< endl;

   return 0;
}
