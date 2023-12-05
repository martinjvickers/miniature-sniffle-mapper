#ifndef PROCESSALIGNMENTS_H
#define PROCESSALIGNMENTS_H

using namespace tbb;

struct OurScore
{
   int trans;
   int as, nm, xs;
   int id;
   bool hasXS = false;
   CharString conversion;
};

struct stats
{
   int reads = 0;
   int uniq_mapped = 0;
   int ambig = 0;
};

stats runstats;

std::map<int,Dna5String> original_reads;
std::map<int,Dna5String> original_reads_R1;
std::map<int,Dna5String> original_reads_R2;
int stepSE = 0;
int stepPE_R1 = 0;
int stepPE_R2 = 0;

// this allows me to compare the CharString BAM result IDs with oneanother
struct MyHashCompare {
    static size_t hash( const CharString& x ) {
        size_t h = 0;
        for(int i = 0; i < length(x); i++)
                h =(h*17)^(char)x[i];
        return h;
    }
    //! True if strings are equal
    static bool equal( const CharString& x, const CharString& y ) {
        return x==y;
    }
};

struct Result
{
   bool done_CT_CT, done_CT_GA, done_GA_CT, done_GA_GA;
   map<int, vector<BamAlignmentRecord>> alignments;
};

typedef concurrent_hash_map<int, Result> TempResultsStore;
TempResultsStore resultsTable;

int getLine(CharString qName)
{
   StringSet<CharString> split;
   strSplit(split, qName, EqualsChar<'_'>());
   int line = strtol(toCString(split[length(split)-1]), NULL, 0);
   return line;
}

void get_mapQ(BamAlignmentRecord &alignment)
{
   float score_min_intercept = 0.0;
   float score_min_slope = -0.2;
   float scMin = score_min_intercept + score_min_slope * length(alignment.seq);
   int diff = abs(scMin);

   BamTagsDict tagsDict(alignment.tags);
   unsigned tagIdx = 0;
   if(!findTagKey(tagIdx, tagsDict, "AS"))
      cerr << "ERROR: Unknown key! There is no AS! " << alignment.qName;

   int AS_best = 0;
   if(!extractTagValue(AS_best, tagsDict, tagIdx))
      cerr << "ERROR: There was an error extracting AS from tags!\n";

   float bestover = AS_best - scMin;

   if(bestover >= diff * 0.8)
      alignment.mapQ = 42;
   else if(bestover >= diff * 0.7)
      alignment.mapQ = 40;
   else if(bestover >= diff * 0.6)
      alignment.mapQ = 24;
   else if(bestover >= diff * 0.5)
      alignment.mapQ = 23;
   else if(bestover >= diff * 0.4)
      alignment.mapQ = 8;
   else if(bestover >= diff * 0.3)
      alignment.mapQ = 3;
   else
      alignment.mapQ = 0;
}

// Given an alignment, return relevant OurScore struct
OurScore get_score(BamAlignmentRecord alignment)
{
   OurScore currscore;
   BamTagsDict tagsDict(alignment.tags);
   unsigned tag_NM_id, tag_XS_id, tag_AS_id, tag_XG_id;
//   currscore.id = id;

   // keys we should 100% have these TAGs
   if (!findTagKey(tag_NM_id, tagsDict, "NM") ||
       !findTagKey(tag_AS_id, tagsDict, "AS") ||
       !findTagKey(tag_XG_id, tagsDict, "XG"))
   {
      cerr << "ERROR: There was an error finding essential tag IDs!\n";
   }

   // extract essential keys
   if(!extractTagValue(currscore.nm, tagsDict, tag_NM_id) ||
      !extractTagValue(currscore.as, tagsDict, tag_AS_id) ||
      !extractTagValue(currscore.conversion, tagsDict, tag_XG_id))
   {
      cerr << "ERROR: There was an error extracting essential tags!\n";
   }

   // extract XS if we have it and set 
   if(findTagKey(tag_XS_id, tagsDict, "XS"))
   {
      extractTagValue(currscore.xs, tagsDict, tag_XS_id);
      currscore.hasXS = true;
   }
   return currscore;
}

bool process_ambig(vector<BamAlignmentRecord> alignments, int &randid)
{
   int rand_gen;

   do
   {
      rand_gen = (rand() % length(alignments));
   }
   while (hasFlagUnmapped(alignments[rand_gen]));

   randid = rand_gen;
   return true;
}

/*
   Bismark best definition!

   1) If there is only a single result, then it's the best.
      1a) But you need to make sure bowtie2 didn't find other results but
          not report them. So check to see if the XS and AS flags are the
          same. If so, then bowtie2 found other results.
      1b) Also, if you have the best mapQ score, you need to check that
          no other result has an AS score that is as good.
      1c) BUT!!! despite 1b), if the mapQ and the AS score are identical
          and postition in the genome is identical you keep that.
          So that is fine. Not 100% sure if it matters which you choose.
*/
int getBest(ModifyStringOptions options, 
            vector<BamAlignmentRecord> &results, 
            TempResultsStore::accessor &a)
{
   // the accessor contains an ordered_map of vectors<BamAlignmentRecord>
   map<int,vector<BamAlignmentRecord>>::iterator iter;
   // the best one with have the score closest to zero
   iter = a->second.alignments.end(); 
   iter--;

   if(options.singleEnd == true)
   {
      // and there is only one entry for this best score
      if(iter->second.size() == 1)
      {
         // now check that bowtie2 didn't actually find secondary results
         OurScore best = get_score(iter->second[0]);
         if(best.hasXS == true && (best.xs != best.as))
         {
            results.push_back(iter->second[0]);
         }
         else if(best.hasXS == false)
         {
            results.push_back(iter->second[0]);
         }
         else
         {
            //cout << "Not entering " << a->first << "\t" << iter->first << endl;
            return 1;
         }
      }
      else if(iter->second.size() == 2)
      {
         // if there is more than one with the same score, it's okay to keep them
         // if it's at the same position
         if(iter->second[0].beginPos == iter->second[1].beginPos &&
            iter->second[0].rID == iter->second[1].rID
           )
         {
            OurScore best = get_score(iter->second[0]);
            if(best.hasXS == true && (best.xs != best.as))
            {
               results.push_back(iter->second[0]);
            }
            else if(best.hasXS == false)
            {
               results.push_back(iter->second[0]);
            }
         } 
         else if(options.ambig == "random")
         {
            int rand = rand % iter->second.size();
            BamAlignmentRecord ambigRead = iter->second[rand];
            BamTagsDict tagsDict(ambigRead.tags);
            setTagValue(tagsDict, "AM", "True");
            results.push_back(ambigRead);
         }
      }
      else
      {
         //cout << "There are more than two entries with the same AS score" << endl;
         // these are the ambig reads so here, 
         // if I simply select a random one, TAG it as ambig, and then return it.
         // The calling function can then read the tag and put it in the correct
         // file.
         if(options.ambig == "random")
         {
            int rand = rand % iter->second.size();
            BamAlignmentRecord ambigRead = iter->second[rand];
            BamTagsDict tagsDict(ambigRead.tags);
            setTagValue(tagsDict, "AM", "True");
            results.push_back(ambigRead);
         }
// This works fine, but needs to be handled correctly when prepping the
// duplicate alignments for writing out
/*         else if(options.ambig == "all")
         {
            for(auto r : iter->second)
            {
               BamTagsDict tagsDict(r.tags);
               setTagValue(tagsDict, "AM", "True");
               results.push_back(r);
               cout << "Pushing " << r.qName << "\tis ambig" << endl;
            }
         }
*/
         return 1; // we don't actually do anything with this
      }
   }
   else if(options.singleEnd == false)
   {

      // find out the number of alignments
      unsigned int count = 0;
      for( auto i : a->second.alignments )
      {
         count = count + i.second.size();
      }

      // the most common alignment case, where a single pair is found.
      if(count == 2)
      {
         //cout << "Count 2 " <<  iter->second[0].qName << endl;

         // horrible but we know there are only two and I can rewrite this
         // at a later point
         OurScore scoreR1, scoreR2;
         unsigned c = 0;

         for(auto i : a->second.alignments)
         {
            for(auto j : i.second)
            {
               if(c == 0)
               {
                  scoreR1 = get_score(j);
               }
               else if(c == 1)
               {
                  scoreR2 = get_score(j);
               }
               else
               {
                  cerr << "WARN: Should never have happened" << endl;
                  return 1; // should never have happened
               }
               c++;
            }
         }

         // now determine if all is well
         if(  ( (scoreR1.hasXS == true && (scoreR1.as > scoreR1.xs)) &&
                (scoreR2.hasXS == true && (scoreR2.as > scoreR2.xs)) ) || 
                (scoreR1.hasXS == false && scoreR2.hasXS == false) ||

              ( (scoreR1.hasXS == true && scoreR2.hasXS == true) &&
                (scoreR1.as + scoreR2.as > scoreR1.xs + scoreR2.xs) ) ||

              ( (scoreR1.hasXS == true && scoreR2.hasXS == false) &&
                (scoreR1.as + scoreR2.as > scoreR1.xs + scoreR2.as) ) ||

              ( (scoreR1.hasXS == false && scoreR2.hasXS == true) &&
                (scoreR1.as + scoreR2.as > scoreR1.as + scoreR2.xs) ) ||

              ( (scoreR1.hasXS == true && scoreR2.hasXS == true) &&
                (scoreR1.as == scoreR1.xs && scoreR2.as < scoreR2.xs) ) ||

              ( (scoreR1.hasXS == true && scoreR2.hasXS == true) &&
                (scoreR1.as < scoreR1.xs && scoreR2.as == scoreR2.xs) )
           )
         {
            for(auto i : a->second.alignments)
               for(auto j : i.second)
                  results.push_back(j);

            return 0;
         }
         else
         {
            return 1;
         }
      }
      else if (count == 4)
      {
         // basically I need them in pairs
         //cout << iter->second[0].qName << " we have four mapped reads!" << endl;
         vector<BamAlignmentRecord> reads;
         
         // This is going to be even worse than count == 2. But again, I will 
         // rewrite this. We know there are exactly four.
         for(auto i : a->second.alignments)
         {
            for(auto j : i.second)
            {
               if(reads.size() < 2)
               {
                  reads.push_back(j);
               }
               else
               {
                  unsigned int c = 1;
                  bool inserted = false;
                  // insert this new read into the location it should be
                  for(auto read : reads)
                  {
                     if(read.beginPos == j.pNext && read.rID == j.rID)
                     {
                        // may need to consider here if first/second made order
                        // is important
                        reads.insert(reads.begin()+c, j);
                        inserted = true;
                     }
                     c++;
                  }
                  
                  // if it didn't match anything already in the vector, 
                  // append to the end
                  if(inserted == false)
                  {
                     reads.push_back(j);
                  }
               }
            }
         }

         //cout << "NUMBER of reads " << reads.size() << endl;

         if(reads.size() == 4)
         {
            OurScore pair1ScoreR1 = get_score(reads[0]);
            OurScore pair1ScoreR2 = get_score(reads[1]);
            OurScore pair2ScoreR1 = get_score(reads[2]);
            OurScore pair2ScoreR2 = get_score(reads[3]);

            if( (pair1ScoreR1.as + pair1ScoreR2.as > pair2ScoreR1.as + pair2ScoreR2.as) &&
                (
                   (pair1ScoreR1.hasXS == false && pair1ScoreR2.hasXS == false) ||
                   
                   ( (pair1ScoreR1.hasXS == true && pair1ScoreR2.hasXS == true) && 
                     (pair1ScoreR1.as > pair1ScoreR1.xs && pair1ScoreR2.as > pair1ScoreR2.xs) ) ||
                
                   ( (pair1ScoreR1.hasXS == true && pair1ScoreR2.hasXS == true) &&
                     (pair1ScoreR1.as + pair1ScoreR2.as > pair1ScoreR1.xs + pair1ScoreR2.xs) ) ||

                   ( (pair1ScoreR1.hasXS == true && pair1ScoreR2.hasXS == false) &&
                     (pair1ScoreR1.as + pair1ScoreR2.as > pair1ScoreR1.xs + pair1ScoreR2.as) ) ||

                   ( (pair1ScoreR1.hasXS == false && pair1ScoreR2.hasXS == true) &&
                     (pair1ScoreR1.as + pair1ScoreR2.as > pair1ScoreR1.as + pair1ScoreR2.xs) )
                )
              )
            {
               results.push_back(reads[0]);
               results.push_back(reads[1]);
               if(options.debug == true)
               {
                  cout << "Pushing Read" << endl;
               }
               return 0;
            }
            else if(  (pair1ScoreR1.as + pair1ScoreR2.as < pair2ScoreR1.as + pair2ScoreR2.as) &&
                      (
                         (pair2ScoreR1.hasXS == false && pair2ScoreR2.hasXS == false) ||

                         ( (pair2ScoreR1.hasXS == true && pair2ScoreR2.hasXS == true) &&
                           (pair2ScoreR1.as > pair2ScoreR1.xs && pair2ScoreR2.as > pair2ScoreR2.xs) ) ||
                        
                         ( (pair2ScoreR1.hasXS == true && pair2ScoreR2.hasXS == true) &&  
                           (pair2ScoreR1.as + pair2ScoreR2.as > pair2ScoreR1.xs + pair2ScoreR2.xs)  ) ||

                         ( (pair2ScoreR1.hasXS == true && pair2ScoreR2.hasXS == false) &&
                           (pair2ScoreR1.as + pair2ScoreR2.as > pair2ScoreR1.xs + pair2ScoreR2.as) ) ||

                         ( (pair2ScoreR1.hasXS == false && pair2ScoreR2.hasXS == true) &&
                           (pair2ScoreR1.as + pair2ScoreR2.as > pair2ScoreR1.as + pair2ScoreR2.xs) )
                      )
                   )
            {
               results.push_back(reads[2]);
               results.push_back(reads[3]);
               if(options.debug == true)
               {
                  cout << "Pushing Read" << endl;
               }
               return 0;
            }
            else // didn't match?
            {
               return 1;
            }
         } 
         else if(reads.size() == 2)
         {
            cerr << reads[0].qName << " are symmetrical?" << endl;
            // presumably here they are symetrical
            return 1;
         }
      }
      else if (count == 3)
      {
         cout << "WARN: 3 mapped reads, something odd? Maybe we just need to take the pair?" << endl;
      }
      else
      {
         return 1;
      }

      /*
      The iterator contains all of the records. So what we probably want to do is
      ensure we know 
      //typedef concurrent_hash_map<int, Result> TempResultsStore;
      struct Result
      {
         bool done_CT_CT, done_CT_GA, done_GA_CT, done_GA_GA;
         map<int, vector<BamAlignmentRecord>> alignments;
      };

      */
   }

   return 0;

}

/*
   Return the record(s) that are to be written out to file.

   Choices about how the records are chosen are called from
   here but once these results are returned, these should be
   written out.

*/
int get_best_alignment(TempResultsStore::accessor &a, 
                       ModifyStringOptions options, 
                       vector<BamAlignmentRecord> &results)
{

   // If the alignment has nothing in it, just exit
   if(a->second.alignments.size() < 1)
      return 1;

   // Using our options, we choose which method of returning our results.
   // The get functions should be aware of single end or paired end!
   getBest(options, results, a); // AKA bismark mode
   //   getBestWithRandom(ModifyStringOptions options, results); 
   // AKA bismark with bs-sequel
   //   getAll(ModifyStringOptions options, results); 
   // simply return everything that mapped!

   if(results.size() > 0)
      return 0;
   else
      return 1;
}

int set_tags(BamAlignmentRecord &record, CharString ref, CharString seq)
{
   BamTagsDict tagsDict(record.tags);

   //if single end
   if(!hasFlagAllProper(record))
   {
      if(ref == "GA") setTagValue(tagsDict,"XR","GA");
      if(ref == "CT") setTagValue(tagsDict,"XR","CT");
      if(seq == "GA") setTagValue(tagsDict,"XG","GA");
      if(seq == "CT") setTagValue(tagsDict,"XG","CT");
   }
   else // if PE?
   {
      if(ref == "CT" && hasFlagFirst(record) && !hasFlagRC(record))
         setTagValue(tagsDict,"XR","CT");
      else if(ref == "CT" && hasFlagLast(record) && hasFlagRC(record))
         setTagValue(tagsDict,"XR","GA");
      else if(ref == "CT" && hasFlagFirst(record) && hasFlagRC(record))
         setTagValue(tagsDict,"XR","CT");
      else if(ref == "CT" && hasFlagLast(record) && !hasFlagRC(record))
         setTagValue(tagsDict,"XR","GA");

      // set the reference conversion
      if(seq == "GA") setTagValue(tagsDict,"XG","GA");
      if(seq == "CT") setTagValue(tagsDict,"XG","CT");
   }

   return 0;
}

// a recursive function to clear out residual unmapped reads from RAM
Dna5String get_read(int line, map<int,Dna5String> &reads, 
                    SeqFileIn &seqFileIn, int &step)
{
   // find line in the map std::map<int,Dna5String> original_reads;
   auto it = reads.find(line);
   if(it != reads.end())
   {
      // return something
      Dna5String result = it->second;
      reads.erase(it);
      //cout << "Found read" << endl;
      return result;
   }
   else
   {
      //cout << "Couldn't find so am going to search for you" << endl;
      String<CharString> ids;
      String<Dna5String> seqs;
      int value = 5000;

      readRecords(ids, seqs, seqFileIn, value);

      // before inserting more, let's remove a bunch of old ones 
      // which presumably haven't mapped
      int keep = 3;
      if(step > keep)
      {
         auto it_low = reads.lower_bound((step-keep-1)*value);
         auto it_high = reads.upper_bound(((step-keep-1)*value)+value);

         for(auto meh = it_low; meh != it_high; meh++)
            reads.erase(meh);
      }

      for(int i = 0; i < length(ids); i++)
      {
         reads.insert(pair<int,Dna5String>((step*value)+i, seqs[i]));
      }
      step++;

      // now, try this again recursivly
      return get_read(line, reads, seqFileIn, step);
   }
}

/*
   This function writes each result to file. It also, depending on the nature of 
   the results and options requested, determine which output file it should be 
   written to, e.g. ambig or final etc.
*/
int writeToFile(ModifyStringOptions options, 
                vector<BamAlignmentRecord> &results,
                BamFileOut &newBamFileOut, bool &headerWritten, 
                SeqFileIn &seqFileIn, SeqFileIn &seqFileInR1, 
                SeqFileIn &seqFileInR2)
{
   if(headerWritten == false)
   {
      BamHeader header;
      /*
        Writing out our flags
      */
      typedef typename BamHeaderRecord::TTag TTag;
      BamHeaderRecord firstRecord;
      firstRecord.type = BAM_HEADER_FIRST;
      appendValue(firstRecord.tags, TTag("VN", "1.0"));
      appendValue(firstRecord.tags, TTag("SO", "unsorted"));
      appendValue(header, firstRecord);

      BamHeaderRecord lastRecord;
      lastRecord.type = BAM_HEADER_PROGRAM;
      appendValue(lastRecord.tags, TTag("ID", "MSM"));
      appendValue(lastRecord.tags, TTag("VN", "0.0.4"));
      appendValue(lastRecord.tags, TTag("CL", options.command));
      appendValue(header, lastRecord);

      writeHeader(newBamFileOut, header);
      headerWritten = true;
   }

   for(auto i : results)
   {
      //writeRecord(newBamFileOut, i);
      BamTagsDict tagsDict(i.tags);
      unsigned tagIdx = 0;

      // if the read is ambig
      // if we have the ambig tag AND it's not ignore
//      if(findTagKey(tagIdx, tagsDict, "AM") && options.ambig != "ignore")
//      {
//         writeRecord(newBamFileOut, i);
//      }
      if(!findTagKey(tagIdx, tagsDict, "AM"))
      {
         writeRecord(newBamFileOut, i);
      }
      else if(findTagKey(tagIdx, tagsDict, "AM") && options.ambig != "ignore")
      {
         writeRecord(newBamFileOut, i);
      }
      // else {} // do nothing, this should be AMBIG tag and options.ambig == "ignore"
   }
   return 0;
}

/*
   This function will prep the alignment record(s) ready for writing out.
   It will;
      * Get the correct orientation - DONE
      * Recover the original sequence - DONE
      * Create the bismark methylation string - DONE
*/
int prepRecords(ModifyStringOptions options, 
                vector<BamAlignmentRecord> &results,
                SeqFileIn &seqFileIn, SeqFileIn &seqFileInR1,
                SeqFileIn &seqFileInR2, FaiIndex &faiIndex)
{
   //cout << "Prepping" << endl;
   for(auto & i : results)
   {
      //cout << i.qName << endl;
      StringSet<CharString> split;
      strSplit(split, i.qName, EqualsChar<'_'>());
      long int line = strtol(toCString(split[length(split)-1]), NULL, 0);
      Dna5String orig;
      if(options.inputFileName != NULL)
      {
         //cout << "It's a Single end one!" << endl;

         /*
            This is the part which is causing a segfault
            when running ambig == "all" because get_read removes the 
            original sequence from the cache meaning when a duplicate
            alignment is encounted, it cannot be found.
         */
         orig = get_read(line, original_reads, seqFileIn, stepSE);
      }
      else if(options.inputLeftFileName != NULL &&
           options.inputRightFileName != NULL &&
           options.inputFileName == NULL)
      {
         //cout << "It's a Paired end one!" << endl;
         // need to determine which mate it is
         if(hasFlagFirst(i))
         {
            //cout << "First mate" << endl;
            orig = get_read(line, original_reads_R1, seqFileInR1, stepPE_R1);
         }
         else if(hasFlagLast(i))
         {
            //cout << "Second mate" << endl;
            orig = get_read(line, original_reads_R2, seqFileInR2, stepPE_R2);
         }
         else
         {
            //cerr << "Odd error" << endl;
            return 1; // shouldn't happen right?
         }
      }

      // now make the Methyl tag
      if(makeMethylTag(i, i.qName, orig, faiIndex) != 0)
      {
         cerr << "Record Not written as it overlaps";
         cerr << " a chromosome edge.\t";
         cerr << i.qName << " ";
         cerr << i.beginPos << " ";
         cerr << i.flag << endl;
         return 1;
      }
   }   
   return 0;
}

// New insert_result function. This should work for both SE/PE reads and not
// care about the number of reads at all! e.g. from bowtie2 all mapped reads
int insert_result(BamAlignmentRecord &record, CharString seq, CharString ref,
                  ModifyStringOptions options, bool &headerWritten,
                  BamFileOut &ambigBamFileOut, BamFileOut &newBamFileOut,
                  FaiIndex &faiIndex, int &last, SeqFileIn &seqFileIn,
                  SeqFileIn &seqFileInR1, SeqFileIn &seqFileInR2)
{
   // find the record in the store, if it doesn't exist, add it
   TempResultsStore::accessor a;
   if(!resultsTable.find(a, getLine(record.qName)))
      resultsTable.insert(a, getLine(record.qName));

   if(options.debug == true)
      cout << "Inserting " << record.qName << " " << seq << ":" << ref << endl;

   // Insert record into map
   if(!hasFlagUnmapped(record))
   {
      BamTagsDict tagsDict(record.tags);
      unsigned tagIdx = 0;
      if(!findTagKey(tagIdx, tagsDict, "AS"))
         cerr << "ERROR: Unknown key! There is no AS tag!! " << record.qName << endl;
      int as;
      if(!extractTagValue(as, tagsDict, tagIdx))
         cerr << "ERROR: There was an error extracting AS from tags!" << endl;
      a->second.alignments[as].push_back(record);
   }

   return 0;
}

/*I'm really terrible at naming things*/
int process(CharString seq, CharString ref,
            ModifyStringOptions options, bool &headerWritten,
            BamFileOut &ambigBamFileOut, BamFileOut &newBamFileOut,
            FaiIndex &faiIndex, int &last, SeqFileIn &seqFileIn,
            SeqFileIn &seqFileInR1, SeqFileIn &seqFileInR2)
{
      TempResultsStore::accessor last_val;
      if(!resultsTable.find(last_val, last))
      {
         //cerr << "Could not find " << last << endl;
      }
      else
      {
         if(seq == "GA" && ref == "GA")
         {
            last_val->second.done_GA_GA = true;
         }
         else if(seq == "GA" && ref == "CT")
         {
            last_val->second.done_GA_CT = true;
         }
         else if(seq == "CT" && ref == "GA")
         {
            last_val->second.done_CT_GA = true;
         }
         else if(seq == "CT" && ref == "CT")
         {
            last_val->second.done_CT_CT = true;
         }

         if( ( last_val->second.done_GA_GA == true &&
               last_val->second.done_GA_CT == true &&
               last_val->second.done_CT_GA == true &&
               last_val->second.done_CT_CT == true &&
               options.non_directional == true) ||
             ( last_val->second.done_CT_GA == true &&
               last_val->second.done_CT_CT == true &&
               options.non_directional == false)
           )
         {
            vector<BamAlignmentRecord> results;
            get_best_alignment(last_val, options, results);
            //cout << "we have the alignments" << endl;
            // get original read and process record
            if(prepRecords(options, results, seqFileIn, seqFileInR1,
                           seqFileInR2, faiIndex) == 0)
            {
               //cout << "prepped" << endl;
               // write out the file
               // It's here that I want to be able to deal with the
               // correct output file, e.g. unmapped, mapped, ambig
               writeToFile(options, results, newBamFileOut, headerWritten,
                           seqFileIn, seqFileInR1, seqFileInR2);
            }

            // remove from memory
            resultsTable.erase(last_val);
         }

      }
   return 0;
}

/*This is my new version*/
void extract_alignments(Iterator<CharString, Rooted>::Type &iterator,
                        CharString ref, CharString seq,
                        BamIOContext<StringSet<CharString>> context,
                        ModifyStringOptions options, bool &headerWritten,
                        BamFileOut &ambigBamFileOut, BamFileOut &newBamFileOut,
                        FaiIndex &faiIndex, int &last, 
                        SeqFileIn &seqFileIn, SeqFileIn &seqFileInR1, 
                        SeqFileIn &seqFileInR2, int done)
{

   while (!atEnd(iterator)) // go over stuff in the iterator
   {
      // read the record
      BamAlignmentRecord record;
      readRecord(record, context, iterator, Sam()); // I need a context?

      if(last == -1)
      {
         //cout << "First run of " << ref << " " << seq << " so adding ";
         last = getLine(record.qName);
         cout << last << endl;
      }

      // set correct tags
      set_tags(record, ref, seq);

      // insert record into result
      insert_result(record, ref, seq, options, headerWritten, ambigBamFileOut,
                    newBamFileOut, faiIndex, last, seqFileIn, seqFileInR1, 
                    seqFileInR2);

      if(getLine(record.qName) != last || done == 1)
      {
         process(ref, seq, options, headerWritten, ambigBamFileOut,
                 newBamFileOut, faiIndex, last, seqFileIn, seqFileInR1,
                 seqFileInR2);
      }

      last = getLine(record.qName);
   }

}

// iterate through the buffer(s)
int parse_buffer(FILE *in_GA_GA, FILE *in_CT_GA, FILE *in_GA_CT, FILE *in_CT_CT,
                 ModifyStringOptions options, BamHeader &header_GA_GA,
                 BamHeader &header_GA_CT, BamHeader &header_CT_GA,
                 BamHeader &header_CT_CT,
                 BamIOContext<StringSet<CharString>> &bamIOContext_GA_GA,
                 BamIOContext<StringSet<CharString>> &bamIOContext_CT_GA,
                 BamIOContext<StringSet<CharString>> &bamIOContext_GA_CT,
                 BamIOContext<StringSet<CharString>> &bamIOContext_CT_CT,
                 BamFileOut &ambigBamFileOut, BamFileOut &newBamFileOut,
                 FaiIndex &faiIndex, SeqFileIn &seqFileIn, 
                 SeqFileIn &seqFileInR1, SeqFileIn &seqFileInR2)
{
   bool headerWritten = false;

   char buff_GA_GA[8192], buff_CT_GA[8192];
   char buff_GA_CT[8192], buff_CT_CT[8192];
   bool done_GA_GA = false, done_CT_GA = false;
   bool done_GA_CT = false, done_CT_CT = false;

   int last_GA_GA = -1, last_CT_GA = -1;
   int last_GA_CT = -1, last_CT_CT = -1;

   while( ( options.non_directional == true &&
            !feof(in_GA_GA) && !feof(in_GA_CT) &&
            !feof(in_CT_GA) && !feof(in_CT_CT) ) ||
          ( options.non_directional == false &&
            !feof(in_CT_GA) && !feof(in_CT_CT) )
        )
   {

      if(options.non_directional == true)
      {
         if(fgets(buff_GA_GA, sizeof(buff_GA_GA), in_GA_GA) == NULL)
            cout << "Finished reading GA_GA." << endl;
         if(fgets(buff_GA_CT, sizeof(buff_GA_CT), in_GA_CT) == NULL)
            cout << "Finished reading GA_CT." << endl;
      }

      if(fgets(buff_CT_GA, sizeof(buff_CT_GA), in_CT_GA) == NULL)
         cout << "Finished reading CT_GA." << endl;
      if(fgets(buff_CT_CT, sizeof(buff_CT_CT), in_CT_CT) == NULL)
         cout << "Finished reading CT_CT." << endl;

         if(options.debug == true)
         {
            cout << "BUFFER " << buff_CT_CT << "\nBUFFER " << buff_CT_GA << endl;
            if(options.non_directional == true)
               cout << "BUFFER " << buff_GA_CT << "\nBUFFER " << buff_GA_GA << endl;
         }

         // convert buffer to charstring
         CharString input_GA_GA = buff_GA_GA;
         CharString input_GA_CT = buff_GA_CT;
         CharString input_CT_GA = buff_CT_GA;
         CharString input_CT_CT = buff_CT_CT;

         // create iterator for the current buffer
         Iterator<CharString, Rooted>::Type iter_GA_GA = begin(input_GA_GA);
         Iterator<CharString, Rooted>::Type iter_GA_CT = begin(input_GA_CT);
         Iterator<CharString, Rooted>::Type iter_CT_GA = begin(input_CT_GA);
         Iterator<CharString, Rooted>::Type iter_CT_CT = begin(input_CT_CT);

         if(options.non_directional == true)
         {
            readHeader(header_GA_GA, bamIOContext_GA_GA, iter_GA_GA, Sam());    
            readHeader(header_GA_CT, bamIOContext_GA_CT, iter_GA_CT, Sam());
         }

         readHeader(header_CT_GA, bamIOContext_CT_GA, iter_CT_GA, Sam());
         readHeader(header_CT_CT, bamIOContext_CT_CT, iter_CT_CT, Sam());

         if(options.non_directional == true)
         {
            extract_alignments(iter_GA_GA, "GA", "GA",
                               bamIOContext_GA_GA, options, headerWritten,
                               ambigBamFileOut, newBamFileOut, faiIndex, 
                               last_GA_GA, seqFileIn, seqFileInR1, 
                               seqFileInR2, feof(in_GA_GA));
            extract_alignments(iter_GA_CT, "GA", "CT",
                               bamIOContext_GA_CT, options, headerWritten,
                               ambigBamFileOut, newBamFileOut, faiIndex, 
                               last_GA_CT, seqFileIn, seqFileInR1, 
                               seqFileInR2, feof(in_GA_CT));
         }

         extract_alignments(iter_CT_GA, "CT", "GA",
                            bamIOContext_CT_GA, options, headerWritten,
                            ambigBamFileOut, newBamFileOut, faiIndex, 
                            last_CT_GA, seqFileIn, seqFileInR1, 
                            seqFileInR2, feof(in_CT_GA));
         extract_alignments(iter_CT_CT, "CT", "CT",
                            bamIOContext_CT_CT, options, headerWritten,
                            ambigBamFileOut, newBamFileOut, faiIndex, 
                            last_CT_CT, seqFileIn, seqFileInR1, 
                            seqFileInR2, feof(in_CT_CT));
   }

   cout << "Final size " << resultsTable.size() << endl;

   return 0;
}

#endif
