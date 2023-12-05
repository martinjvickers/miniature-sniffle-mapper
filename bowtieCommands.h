#ifndef BOWTIECOMMANDS_H 
#define BOWTIECOMMANDS_H

int prepareBowtie2Commands(ModifyStringOptions options, string &bowtie_GA_GA, 
                           string &bowtie_GA_CT, string &bowtie_CT_GA, 
                           string &bowtie_CT_CT)
{
   CharString converted_CT = toCString(options.inputFileName) 
                           + string("_ref.temp.C2T.fastq");
   CharString converted_GA = toCString(options.inputFileName) 
                           + string("_ref.temp.G2A.fastq");

   if(options.compress == true)
   {
      converted_CT = toCString(converted_CT) + string(".gz");
      converted_GA = toCString(converted_GA) + string(".gz");
   }

   // single End
   if(options.singleEnd == true)
   {
      /* If non-directional*/
      if(options.non_directional == true)
      {
         bowtie_GA_GA = string("bowtie2 -q --score-min L,0,-0.2 --ignore-quals")
                      + (" --norc")
                      + (" -p ") + to_string(options.cpu_cores)
                      + (" -x ") + toCString(options.referenceFileName)+("_G2A")
                      + (" -U ") + toCString(converted_GA);

         bowtie_GA_CT = string("bowtie2 -q --score-min L,0,-0.2 --ignore-quals")
                      + (" --nofw")
                      + (" -p ") + to_string(options.cpu_cores)
                      + (" -x ") + toCString(options.referenceFileName)+("_C2T")
                      + (" -U ") + toCString(converted_GA);
      }

      bowtie_CT_GA = string("bowtie2 -q --score-min L,0,-0.2 --ignore-quals")
                   + (" --nofw")
                   + (" -p ") + to_string(options.cpu_cores)
                   + (" -x ") + toCString(options.referenceFileName) + ("_G2A")
                   + (" -U ") + toCString(converted_CT);

      bowtie_CT_CT = string("bowtie2 -q --score-min L,0,-0.2 --ignore-quals")
                   + (" --norc")
                   + (" -p ") + to_string(options.cpu_cores)
                   + (" -x ") + toCString(options.referenceFileName) + ("_C2T")
                   + (" -U ") + toCString(converted_CT);
   }
   else // paired end
   {
      if(options.non_directional == true)
      {
         bowtie_GA_GA;
         bowtie_GA_CT;
      }

      bowtie_CT_GA = string("bowtie2 -q --score-min L,0,-0.2 --ignore-quals")
                   + (" --no-mixed --no-discordant --dovetail")
                   + (" --maxins 500 --nofw")
                   + (" -p ") + to_string(options.cpu_cores)
                   + (" -x ") + toCString(options.referenceFileName) + ("_G2A")
                   + (" -1 ") + toCString(options.inputLeftFileName)
                   + ("_ref.temp.C2T.fastq")
                   + (" -2 ") + toCString(options.inputRightFileName)
                   + ("_ref.temp.G2A.fastq");

      bowtie_CT_CT = string("bowtie2 -q --score-min L,0,-0.2 --ignore-quals")
                   + (" --no-mixed --no-discordant --dovetail")
                   + (" --maxins 500 --norc")
                   + (" -p ") + to_string(options.cpu_cores)
                   + (" -x ") + toCString(options.referenceFileName) + ("_C2T")
                   + (" -1 ") + toCString(options.inputLeftFileName)
                   + ("_ref.temp.C2T.fastq")
                   + (" -2 ") + toCString(options.inputRightFileName)
                   + ("_ref.temp.G2A.fastq");
   }

   cout << "Begining to map " << endl;
   cout << bowtie_GA_GA << endl;
   cout << bowtie_GA_CT << endl;
   cout << bowtie_CT_GA << endl;
   cout << bowtie_CT_CT << endl;

   return 0;
}

#endif
