#ifndef PREPAREREADS_H 
#define PREPAREREADS_H

#include <string>

using namespace std;
using namespace seqan;

struct ConvertCT : public std::unary_function<Dna5, Dna5>
{
   inline Dna5 operator()(Dna5 x) const
   {
      if(x == 'C')
         return 'T';
      return x;
   }

};

struct ConvertGA : public std::unary_function<Dna5, Dna5>
{
   inline Dna5 operator()(Dna5 x) const
   {
      if(x == 'G')
         return 'A';
      return x;
   }
};

int cleanup(ModifyStringOptions options, CharString fileToProcess)
{
   string ct_out = toCString(fileToProcess) + string("_ref.temp.C2T.fastq");
   string ga_out = toCString(fileToProcess) + string("_ref.temp.G2A.fastq");

   if(options.compress == true)
   {
      ct_out = ct_out + string(".gz");
      ga_out = ga_out + string(".gz");
   }

   cout << "Deleting " << ga_out << endl;
   
   if(options.non_directional == true)
   {
      if(remove(ga_out.c_str()) != 0)
         return 1;
   }

   cout << "Deleting " << ct_out << endl;

   if(remove(ct_out.c_str()) != 0)
      return 1;

   return 0;
}

int process_thread(ModifyStringOptions options, CharString fileToProcess, string conversion)
{
   //read in the fastq file in one go
   CharString updatedid, newqual;
   String<Dna5> newseq;
   SeqFileIn seqFileIn;
   SeqFileOut seqFileOut;
   unsigned int line_number = 0;

   string file_out;

   if(conversion == "CT")
      file_out = toCString(fileToProcess) + string("_ref.temp.C2T.fastq");
   else if(conversion == "GA")
      file_out = toCString(fileToProcess) + string("_ref.temp.G2A.fastq");
   
   if(options.compress == true)
   {
      file_out = file_out + string(".gz");
   }

   if(!open(seqFileIn, toCString(fileToProcess), seqan::OPEN_RDONLY) ||
      !open(seqFileOut, toCString(file_out)))
   {
      cerr << "ERROR: Could not open the file.\n";
      return 1;
   }

   while(!atEnd(seqFileIn))
   {
      try
      {
         readRecord(updatedid, newseq, newqual, seqFileIn);
      }
      catch (Exception const & e)
      {
         cout << "ERROR: " << e.what() << std::endl;
         return 1;
      }

      typedef ModifiedString<String<Dna5>, ModView<ConvertCT> > TModCT;
      typedef ModifiedString<String<Dna5>, ModView<ConvertGA> > TModGA;

      /*
         I'm not especially keen on this method but I'm not sure how else to 
         tackle it. add line number to id string in order to get that info 
         later on to fix an odd bug (between SAM/BAM and fastq records, fastq 
         recognises spaces in IDS, SAM/BAMs dont)
         I need to take only the first part of the ID and clean the ID now 
         rather than later as at this point the line_number is added to a part
         that the SAM/BAM ignores. (ARGH)
      */
      StringSet<CharString> split;
      strSplit(split, updatedid, seqan::EqualsChar<' '>());
      string newid = toCString(split[0]) + string("_");
      newid = newid + to_string(line_number);
      line_number++;

      //convert G's to A's and C's to T's
      if(conversion == "CT")
      {
         TModCT modCT(newseq);
         writeRecord(seqFileOut, newid, modCT, newqual);
      }
      else if(conversion == "GA")
      {
         TModGA modGA(newseq);
         writeRecord(seqFileOut, newid, modGA, newqual);
      }
   } 

   close(seqFileOut), close(seqFileIn);

   return 0;
}

int process_reads(ModifyStringOptions options, CharString fileToProcess)
{
   vector<thread> vectorOfThreads;

   vectorOfThreads.push_back(thread(process_thread, options, fileToProcess, "CT"));
   
   if(options.non_directional == true)
      vectorOfThreads.push_back(thread(process_thread, options, fileToProcess, "GA"));

   for(auto &thread : vectorOfThreads)
      thread.join();

   return 0;
}

#endif
