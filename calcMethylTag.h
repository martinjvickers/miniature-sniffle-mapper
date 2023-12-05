#ifndef CALCMETHYLTAG_H 
#define CALCMETHYLTAG_H

/*
   z unmethylated C in CpG context
   Z methylated C in CpG context
   x unmethylated C in CHG context
   X methylated C in CHG context
   h unmethylated C in CHH context
   H methylated C in CHH context   
*/
int makeMethylTag(BamAlignmentRecord &record, CharString id, Dna5String fragseq,
                  FaiIndex &faiIndex)
{
   BamTagsDict tagsDict(record.tags);
   unsigned tagIdx = 0;
   CharString XR_val, XG_val;

   // get tags
   if(!findTagKey(tagIdx, tagsDict, "XR"))
      cerr << "ERROR: Unknown key! There is no XR! " << record.qName << endl;
   if(!extractTagValue(XR_val, tagsDict, tagIdx))
      cerr << "ERROR: There was an error extracting XR from tags!" << endl;
   if(!findTagKey(tagIdx, tagsDict, "XG"))
      cerr << "ERROR: Unknown key! There is no XG! " << record.qName << endl;
   if(!extractTagValue(XG_val, tagsDict, tagIdx))
      cerr << "ERROR: There was an error extracting XG from tags!" << endl;

   // create blank methylation string
   CharString methyl_tag;
   for(unsigned i = 0; i < length(record.seq); ++i)
      methyl_tag += '.';

   // get the reference with +2 on either side
   Dna5String reference_seq;

   // quick fix, if it's single end, do this
   if(!hasFlagAllProper(record))
   {
      if(XR_val == "GA" && XG_val == "GA")
      {
         record.flag = record.flag ^ BAM_FLAG_RC;
      }
      else if(XR_val == "CT" && XG_val == "CT")
      {
      }
      else if(XR_val == "CT" && XG_val == "GA")
      {
         reverseComplement(fragseq);
         reverseComplement(record.seq);
      }
      else if(XR_val == "GA" && XG_val == "CT")
      {
         reverseComplement(fragseq);
         reverseComplement(record.seq);
         record.flag = record.flag ^ BAM_FLAG_RC;
      }
   }
   else if(hasFlagAllProper(record))
   {
      if(XR_val == "GA" && XG_val == "CT")
      {
         reverseComplement(fragseq);
         reverseComplement(record.seq);
      }
      else if(XR_val == "CT" && XG_val == "GA")
      {
         reverseComplement(fragseq);
         reverseComplement(record.seq);
      }
   }
  

   //quick -- sorry -- hacky hack hack :-(
   //calculate any additional insertions required to the reference
   int inserted_bp = 0;
   if(length(record.cigar) != 1)   // seq matches the reference completely
   {
      for(unsigned j = 0; j < length(record.cigar); j++)
         if(record.cigar[j].operation == 'D')
            inserted_bp = inserted_bp + record.cigar[j].count;
   }

   // for majority of cases, this just works, so why mess with it.
   readRegion(reference_seq, faiIndex, record.rID, record.beginPos - 2,
              record.beginPos + length(record.seq) + 2 + inserted_bp);

   /*
   if we have a reverse compliment flag then we need to look at the -2 -1 
   values of our reference so die
        //if(length(reference_seq) == 0 && hasFlagRC(record))
        //      return 1;

        //if we don't have a RC flag, then we don't. But for the below to 
        work, we need to pad the first 
        //two bases... how about we do it with N's
        //so, let's pad the number we need
        //if(length(reference_seq) == 0 && !hasFlagRC(record))
   */
   if(length(reference_seq) == 0)
   {
      if(record.beginPos - 2 == -1)
      {
         readRegion(reference_seq, faiIndex, record.rID, record.beginPos - 1,
                    record.beginPos + length(record.seq) + 2 + inserted_bp);
         Dna5String pad = "N";
         pad += reference_seq;
         reference_seq = pad;
      }
      else if(record.beginPos - 2 == -2)
      {
         readRegion(reference_seq, faiIndex, record.rID, record.beginPos,
                    record.beginPos + length(record.seq) + 2 + inserted_bp);
         Dna5String pad = "NN";
         pad += reference_seq;
         reference_seq = pad;
      }
      else
      {
         cerr << "ERROR: I've clearly done something wrong here" << endl;
      }
   }

   /*
   So, after all this faffing, if we don't have a record that is 2 bases on 
   either side of the sequence, we've failed
   */
   if(record.beginPos > 0 &&
      !(length(reference_seq) >= length(record.seq) + 4) &&
      !hasFlagRC(record))
      return 1;
   else if(hasFlagRC(record) && (record.beginPos-1 < 0))
      return 1;
   else if(hasFlagRC(record) && (record.beginPos-2 < 0))
      return 1;

   Dna5String corrected_ref = reference_seq;

   // seq matches the reference completely
   if(length(record.cigar) != 1)
   {
      int pos = 0;
      for(unsigned j = 0; j < length(record.cigar); j++)
      {
         if(record.cigar[j].operation == 'M')
         {
            pos = pos + record.cigar[j].count;
         }
         else if(record.cigar[j].operation == 'I')
         {
            for(int i = 0; i < record.cigar[j].count; i++)
               insert(corrected_ref, pos+2, 'N');
            pos = pos + record.cigar[j].count;
         }
         else if(record.cigar[j].operation == 'D')
         {
            for(int i = 0; i < record.cigar[j].count; i++)
               erase(corrected_ref, pos+2);
         }
         else
         {
            cout << "CIGAR ISSUE? Found " << record.cigar[j].operation << endl;
         }
      }
   }

   // now go through the sequence and take a look
   for(unsigned i = 2; i < length(record.seq)+2; ++i) // correct
   {
      if((!hasFlagRC(record) && !hasFlagAllProper(record) ) || 
         (hasFlagRC(record) && hasFlagAllProper(record) && hasFlagLast(record)) ||
         (!hasFlagRC(record) && hasFlagAllProper(record) && hasFlagFirst(record))
        )
      {
         if((fragseq[i-2] == 'T') && (corrected_ref[i] == 'C'))
         {
            if(corrected_ref[i+1] == 'G')
               methyl_tag[i-2] = 'z';   //CG
            else if(corrected_ref[i+1] == 'N')
               methyl_tag[i-2] = 'u';   //Unknown
            else
               if(corrected_ref[i+2] == 'G')
                  methyl_tag[i-2] = 'x';    //CHG
               else if(corrected_ref[i+2] == 'N')
                  methyl_tag[i-2] = 'u';
               else
                  methyl_tag[i-2] = 'h';    //CHH
         }
         else if((fragseq[i-2] == 'C') && (corrected_ref[i] == 'C'))
         {
            if(corrected_ref[i+1] == 'G')
               methyl_tag[i-2] = 'Z';  //CG
            else if(corrected_ref[i+1] == 'N')
               methyl_tag[i-2] = 'U';  //Unknown
            else
               if(corrected_ref[i+2] == 'G')
                  methyl_tag[i-2] = 'X';    //CHG
               else if(corrected_ref[i+2] == 'N')
                  methyl_tag[i-2] = 'U';
               else
                  methyl_tag[i-2] = 'H';    //CHH
         }
      }
      else if((hasFlagRC(record) && !hasFlagAllProper(record)) || 
              (!hasFlagRC(record) && hasFlagAllProper(record) && hasFlagLast(record)) || 
              (hasFlagRC(record) && hasFlagAllProper(record) && hasFlagFirst(record) )
             )
      {
         if((fragseq[i-2] == 'A') && (corrected_ref[i] == 'G'))
         {
            if(corrected_ref[i-1] == 'C')
               methyl_tag[i-2] = 'z';  //CG
            else if(corrected_ref[i-1] == 'N')
               methyl_tag[i-2] = 'u';  //Unknown
            else
               if(corrected_ref[i-2] == 'C')
                  methyl_tag[i-2] = 'x';    //CHG
               else if(corrected_ref[i-2] == 'N')
                  methyl_tag[i-2] = 'u';
               else
                  methyl_tag[i-2] = 'h';    //CHH
         }
         else if((fragseq[i-2] == 'G') && (corrected_ref[i] == 'G'))
         {
            if(corrected_ref[i-1] == 'C')
               methyl_tag[i-2] = 'Z';  //CG
            else if(corrected_ref[i-1] == 'N')
               methyl_tag[i-2] = 'U';  //Unknown
            else
               if(corrected_ref[i-2] == 'C')
                  methyl_tag[i-2] = 'X';    //CHG
               else if(corrected_ref[i-2] == 'N')
                  methyl_tag[i-2] = 'U';
               else
                  methyl_tag[i-2] = 'H';    //CHH
         }
      }
   }

   setTagValue(tagsDict,"XM",methyl_tag);
   record.seq = fragseq;

   return 0;
}

#endif
