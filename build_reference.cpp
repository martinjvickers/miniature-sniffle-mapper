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

//overload the SeqFileBuffer_ so that it uses Iupac String. In this way 
//the input file is checked against Iupac and any non-A/C/G/T is silently 
//converted into a N.
namespace seqan {
	template <typename TString, typename TSSetSpec, typename TSpec>
	struct SeqFileBuffer_<StringSet<TString, TSSetSpec>, TSpec>
	{
		typedef String<Iupac> Type;
	};
}

struct ModifyStringOptions
{
        CharString inputFileName;
};

seqan::ArgumentParser::ParseResult parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
	seqan::ArgumentParser parser("build_reference");
	addOption(parser, seqan::ArgParseOption("i", "input-file", "Path to the input file", seqan::ArgParseArgument::INPUT_FILE, "IN"));
	setRequired(parser, "input-file");
	setShortDescription(parser, "BS Mapper");
	setVersion(parser, "0.0.5");
	setDate(parser, "September 2018");
	addUsageLine(parser, "-i sequence.fastq [\\fIOPTIONS\\fP] ");

	addDescription(parser, "Create the unbiased BS reference.");
	seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

	// If parsing was not successful then exit with code 1 if there were errors.
	// Otherwise, exit with code 0 (e.g. help was printed).
	if (res != seqan::ArgumentParser::PARSE_OK)
		return res;

	getOptionValue(options.inputFileName, parser, "input-file");

	return seqan::ArgumentParser::PARSE_OK;

}

struct ConvertCT :
    public std::unary_function<Iupac, Iupac>
{
    inline Iupac operator()(Iupac x) const
    {
        if (x == 'C') return 'T';

        return x;
    }

};

struct ConvertGA :
    public std::unary_function<Iupac, Iupac>
{
    inline Iupac operator()(Iupac x) const
    {
        if (x == 'G') return 'A';

        return x;
    }
};

/*
*/
int main(int argc, char const ** argv)
{

	//parse our options
	ModifyStringOptions options;
	seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

	if (res != seqan::ArgumentParser::PARSE_OK)
		return res == seqan::ArgumentParser::PARSE_ERROR;

	//read in all fastqfile
	CharString id;
	//Dna5String seq;
	IupacString seq;
	SeqFileIn seqFileIn(toCString(options.inputFileName));
	
	//open two fileout handlers for C2T and G2A
	SeqFileOut seqFileOutC2T;
	string c2t_outname = toCString(options.inputFileName) + string("_C2T.fa");
        if (!open(seqFileOutC2T, toCString(c2t_outname)))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }

        SeqFileOut seqFileOutG2A;
	string g2a_outname = toCString(options.inputFileName) + string("_G2A.fa");
        if (!open(seqFileOutG2A, toCString(g2a_outname)))
        {
                std::cerr << "ERROR: Could not open the file.\n";
                return 1;
        }
	while(!atEnd(seqFileIn))
	{
        	try
        	{
        	        readRecord(id, seq, seqFileIn);
        	}
        	catch (Exception const & e)
        	{
        	        std::cout << "ERROR: " << e.what() << std::endl;
        	        return 1;
        	}
		//convert C's to T's
		typedef ModifiedString<IupacString, ModView<ConvertCT> > TModCT;
		TModCT modCT(seq);	
		writeRecord(seqFileOutC2T, id, modCT);

		//convert G's to A's
		typedef ModifiedString<IupacString, ModView<ConvertGA> > TModGA;
		TModGA modGA(seq);
		writeRecord(seqFileOutG2A, id, modGA);

	}

	close(seqFileOutG2A);
	close(seqFileOutC2T);
	//open each file and run bowtie-build on it;
	string const bowtie_cmd_G2A = string("bowtie2-build  ") + g2a_outname + " " + toCString(options.inputFileName) + string("_G2A");
	system(bowtie_cmd_G2A.c_str());
	string const bowtie_cmd_C2T = string("bowtie2-build  ") + c2t_outname + " " + toCString(options.inputFileName) + string("_C2T");
	system(bowtie_cmd_C2T.c_str());
	return 0;
}
