// Copyright 2023 Jimmy Eng
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "NovelModeUtils.h"
#include "CometDataInternal.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <unordered_set>

extern comet_fileoffset_t clSizeCometFileOffset;

static string TrimStringLocal(const string& sInput)
{
   size_t iStart = 0;
   while (iStart < sInput.length() && isspace((unsigned char)sInput[iStart]))
      ++iStart;

   size_t iEnd = sInput.length();
   while (iEnd > iStart && isspace((unsigned char)sInput[iEnd - 1]))
      --iEnd;

   return sInput.substr(iStart, iEnd - iStart);
}

bool ParsePeptideIdxEntries(const string& sIdxPath,
                            vector<PeptideMassEntry>& vEntries,
                            bool bLoadProteinIds,
                            string& sErrorMsg)
{
   vEntries.clear();
   FILE* fp = fopen(sIdxPath.c_str(), "rb");
   if (fp == NULL)
   {
      sErrorMsg = " Error - cannot open peptide index file \"" + sIdxPath + "\".\n";
      return false;
   }

   comet_fileoffset_t lEndOfStruct = 0;
   comet_fileoffset_t lProteinsPos = 0;
   if (comet_fseek(fp, -clSizeCometFileOffset * 2, SEEK_END) != 0
         || fread(&lEndOfStruct, clSizeCometFileOffset, 1, fp) != 1
         || fread(&lProteinsPos, clSizeCometFileOffset, 1, fp) != 1)
   {
      fclose(fp);
      sErrorMsg = " Error - invalid peptide index footer in \"" + sIdxPath + "\".\n";
      return false;
   }

   vector<vector<string>> vProteinSets;
   if (bLoadProteinIds)
   {
      if (lProteinsPos <= 0 || comet_fseek(fp, lProteinsPos, SEEK_SET) != 0)
      {
         fclose(fp);
         sErrorMsg = " Error - invalid peptide index protein section in \"" + sIdxPath + "\".\n";
         return false;
      }

      size_t tNumSets = 0;
      if (fread(&tNumSets, clSizeCometFileOffset, 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - failed reading protein-set count from \"" + sIdxPath + "\".\n";
         return false;
      }

      vProteinSets.resize(tNumSets);
      for (size_t iSet = 0; iSet < tNumSets; ++iSet)
      {
         size_t tSetSize = 0;
         if (fread(&tSetSize, clSizeCometFileOffset, 1, fp) != 1)
         {
            fclose(fp);
            sErrorMsg = " Error - failed reading protein-set size from \"" + sIdxPath + "\".\n";
            return false;
         }

         vector<comet_fileoffset_t> vOffsets;
         vOffsets.resize(tSetSize);
         for (size_t i = 0; i < tSetSize; ++i)
         {
            if (fread(&vOffsets[i], clSizeCometFileOffset, 1, fp) != 1)
            {
               fclose(fp);
               sErrorMsg = " Error - failed reading protein offset from \"" + sIdxPath + "\".\n";
               return false;
            }
         }

         unordered_set<string> setSeen;
         for (size_t i = 0; i < vOffsets.size(); ++i)
         {
            if (comet_fseek(fp, vOffsets[i], SEEK_SET) != 0)
            {
               fclose(fp);
               sErrorMsg = " Error - failed seeking protein entry in \"" + sIdxPath + "\".\n";
               return false;
            }

            char szProteinName[512];
            if (fscanf(fp, "%500s", szProteinName) != 1)
            {
               fclose(fp);
               sErrorMsg = " Error - failed reading protein name from \"" + sIdxPath + "\".\n";
               return false;
            }
            szProteinName[500] = '\0';

            string sProteinId = szProteinName;
            if (!sProteinId.empty() && setSeen.insert(sProteinId).second)
               vProteinSets[iSet].push_back(sProteinId);
         }
      }
   }

   if (lEndOfStruct <= 0)
   {
      fclose(fp);
      sErrorMsg = " Error - invalid peptide index structure offset in \"" + sIdxPath + "\".\n";
      return false;
   }

   if (comet_fseek(fp, lEndOfStruct, SEEK_SET) != 0)
   {
      fclose(fp);
      sErrorMsg = " Error - cannot seek peptide index structure in \"" + sIdxPath + "\".\n";
      return false;
   }

   int iMinMass = 0;
   int iMaxMass = 0;
   uint64_t tNumPeptides = 0;
   if (fread(&iMinMass, sizeof(int), 1, fp) != 1
         || fread(&iMaxMass, sizeof(int), 1, fp) != 1
         || fread(&tNumPeptides, sizeof(uint64_t), 1, fp) != 1)
   {
      fclose(fp);
      sErrorMsg = " Error - failed reading peptide index summary for \"" + sIdxPath + "\".\n";
      return false;
   }
   (void)iMinMass;
   (void)tNumPeptides;

   if (iMaxMass <= 0)
   {
      fclose(fp);
      return true;
   }

   int iMaxMass10 = iMaxMass * 10;
   vector<comet_fileoffset_t> vIndex(iMaxMass10, -1);
   if (fread(vIndex.data(), clSizeCometFileOffset, iMaxMass10, fp) != (size_t)iMaxMass10)
   {
      fclose(fp);
      sErrorMsg = " Error - failed reading peptide index offsets from \"" + sIdxPath + "\".\n";
      return false;
   }

   comet_fileoffset_t lPeptideStart = -1;
   for (int i = 0; i < iMaxMass10; ++i)
   {
      if (vIndex[i] >= 0 && (lPeptideStart == -1 || vIndex[i] < lPeptideStart))
         lPeptideStart = vIndex[i];
   }
   if (lPeptideStart < 0)
   {
      fclose(fp);
      return true;
   }

   if (comet_fseek(fp, lPeptideStart, SEEK_SET) != 0)
   {
      fclose(fp);
      sErrorMsg = " Error - failed seeking peptide entries in \"" + sIdxPath + "\".\n";
      return false;
   }

   while (comet_ftell(fp) < lEndOfStruct)
   {
      int iLen = 0;
      if (fread(&iLen, sizeof(int), 1, fp) != 1)
         break;

      if (iLen <= 0 || iLen >= MAX_PEPTIDE_LEN)
      {
         fclose(fp);
         sErrorMsg = " Error - invalid peptide length while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      string sPeptide;
      sPeptide.resize((size_t)iLen);
      if (fread(&sPeptide[0], sizeof(char), (size_t)iLen, fp) != (size_t)iLen)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated peptide sequence while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      char cPrev = '-';
      char cNext = '-';
      if (fread(&cPrev, sizeof(char), 1, fp) != 1
            || fread(&cNext, sizeof(char), 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated peptide flanking residues while parsing \"" + sIdxPath + "\".\n";
         return false;
      }
      (void)cPrev;
      (void)cNext;

      unsigned char cNumMods = 0;
      if (fread(&cNumMods, sizeof(unsigned char), 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated peptide mod count while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      for (unsigned char iMod = 0; iMod < cNumMods; ++iMod)
      {
         unsigned char cPos = 0;
         char cWhichMod = 0;
         if (fread(&cPos, sizeof(unsigned char), 1, fp) != 1
               || fread(&cWhichMod, sizeof(char), 1, fp) != 1)
         {
            fclose(fp);
            sErrorMsg = " Error - truncated peptide mod sites while parsing \"" + sIdxPath + "\".\n";
            return false;
         }
         (void)cPos;
         (void)cWhichMod;
      }

      double dPepMass = 0.0;
      comet_fileoffset_t lProteinSetPos = 0;
      if (fread(&dPepMass, sizeof(double), 1, fp) != 1
            || fread(&lProteinSetPos, clSizeCometFileOffset, 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated peptide mass/protein-set entry while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      PeptideMassEntry entry;
      entry.sPeptide = sPeptide;
      entry.dMass = dPepMass;
      entry.bDecoy = false;
      if (bLoadProteinIds)
      {
         if (lProteinSetPos < 0 || (size_t)lProteinSetPos >= vProteinSets.size())
         {
            fclose(fp);
            sErrorMsg = " Error - invalid protein-set reference while parsing \"" + sIdxPath + "\".\n";
            return false;
         }

         entry.vProteinIds = vProteinSets[(size_t)lProteinSetPos];
      }
      vEntries.push_back(entry);
   }

   fclose(fp);
   return true;
}

bool ParseFragmentIdxEntries(const string& sIdxPath,
                             vector<PeptideMassEntry>& vEntries,
                             string& sErrorMsg)
{
   vEntries.clear();
   FILE* fp = fopen(sIdxPath.c_str(), "rb");
   if (fp == NULL)
   {
      sErrorMsg = " Error - cannot open fragment index file \"" + sIdxPath + "\".\n";
      return false;
   }

   comet_fileoffset_t clPeptidesFilePos = 0;
   comet_fileoffset_t clProteinsFilePos = 0;
   comet_fileoffset_t clPermutationsFilePos = 0;
   if (comet_fseek(fp, -clSizeCometFileOffset * 3, SEEK_END) != 0
         || fread(&clPeptidesFilePos, clSizeCometFileOffset, 1, fp) != 1
         || fread(&clProteinsFilePos, clSizeCometFileOffset, 1, fp) != 1
         || fread(&clPermutationsFilePos, clSizeCometFileOffset, 1, fp) != 1)
   {
      fclose(fp);
      sErrorMsg = " Error - invalid fragment index footer in \"" + sIdxPath + "\".\n";
      return false;
   }
   (void)clProteinsFilePos;
   (void)clPermutationsFilePos;

   if (clPeptidesFilePos <= 0 || comet_fseek(fp, clPeptidesFilePos, SEEK_SET) != 0)
   {
      fclose(fp);
      sErrorMsg = " Error - invalid fragment peptide section in \"" + sIdxPath + "\".\n";
      return false;
   }

   size_t tNumPeptides = 0;
   if (fread(&tNumPeptides, sizeof(size_t), 1, fp) != 1)
   {
      fclose(fp);
      sErrorMsg = " Error - failed reading fragment peptide count from \"" + sIdxPath + "\".\n";
      return false;
   }

   for (size_t i = 0; i < tNumPeptides; ++i)
   {
      int iLen = 0;
      if (fread(&iLen, sizeof(int), 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated fragment peptide entry while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      if (iLen <= 0 || iLen >= MAX_PEPTIDE_LEN)
      {
         fclose(fp);
         sErrorMsg = " Error - invalid fragment peptide length while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      string sPeptide;
      sPeptide.resize((size_t)iLen);
      if (fread(&sPeptide[0], sizeof(char), (size_t)iLen, fp) != (size_t)iLen)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated fragment peptide sequence while parsing \"" + sIdxPath + "\".\n";
         return false;
      }

      char cPrev = '-';
      char cNext = '-';
      double dPepMass = 0.0;
      unsigned short siVarModProteinFilter = 0;
      comet_fileoffset_t lProteinSetPos = 0;

      if (fread(&cPrev, sizeof(char), 1, fp) != 1
            || fread(&cNext, sizeof(char), 1, fp) != 1
            || fread(&dPepMass, sizeof(double), 1, fp) != 1
            || fread(&siVarModProteinFilter, sizeof(unsigned short), 1, fp) != 1
            || fread(&lProteinSetPos, clSizeCometFileOffset, 1, fp) != 1)
      {
         fclose(fp);
         sErrorMsg = " Error - truncated fragment peptide metadata while parsing \"" + sIdxPath + "\".\n";
         return false;
      }
      (void)cPrev;
      (void)cNext;
      (void)siVarModProteinFilter;
      (void)lProteinSetPos;

      PeptideMassEntry entry;
      entry.sPeptide = sPeptide;
      entry.dMass = dPepMass;
      entry.bDecoy = false;
      vEntries.push_back(entry);
   }

   fclose(fp);
   return true;
}

bool WritePeptidesToFasta(const vector<string>& vPeptides,
                          const string& sPrefix,
                          string& sOutPath,
                          string& sErrorMsg)
{
   if (!CreateTempPath(sPrefix, ".fasta", sOutPath, sErrorMsg))
      return false;

   std::ofstream outFile(sOutPath.c_str(), std::ios::out | std::ios::trunc);
   if (!outFile.good())
   {
      sErrorMsg = " Error - cannot create temporary FASTA \"" + sOutPath + "\".\n";
      return false;
   }

   unordered_set<string> setSeen;
   int iWritten = 0;
   for (size_t i = 0; i < vPeptides.size(); ++i)
   {
      string sPep = NormalizePeptideToken(vPeptides.at(i));
      if (sPep.empty())
         continue;
      if (!setSeen.insert(sPep).second)
         continue;

      outFile << ">COMETPLUS_NOVEL_" << (iWritten + 1) << "\n";
      outFile << sPep << "\n";
      iWritten++;
   }

   outFile.flush();
   if (!outFile.good())
   {
      sErrorMsg = " Error - failed writing temporary FASTA \"" + sOutPath + "\".\n";
      return false;
   }

   return true;
}

bool WriteNovelRecordsToFasta(const vector<NovelPeptideRecord>& vRecords,
                              const string& sPrefix,
                              string& sOutPath,
                              string& sErrorMsg)
{
   if (!CreateTempPath(sPrefix, ".fasta", sOutPath, sErrorMsg))
      return false;

   std::ofstream outFile(sOutPath.c_str(), std::ios::out | std::ios::trunc);
   if (!outFile.good())
   {
      sErrorMsg = " Error - cannot create temporary FASTA \"" + sOutPath + "\".\n";
      return false;
   }

   unordered_set<string> setSeenPeptides;
   unordered_set<string> setSeenIds;
   int iAutoId = 1;

   for (size_t i = 0; i < vRecords.size(); ++i)
   {
      string sPep = NormalizePeptideToken(vRecords.at(i).sPeptide);
      if (sPep.empty())
         continue;
      if (!setSeenPeptides.insert(sPep).second)
         continue;

      string sPepId = TrimStringLocal(vRecords.at(i).sPeptideId);
      if (sPepId.empty())
      {
         sPepId = "COMETPLUS_NOVEL_" + std::to_string(iAutoId);
         iAutoId++;
      }

      if (!setSeenIds.insert(sPepId).second)
      {
         sErrorMsg = " Error - duplicate peptide_id encountered while building novel scoring FASTA: \"" + sPepId + "\".\n";
         return false;
      }

      outFile << ">" << sPepId << "\n";
      outFile << sPep << "\n";
   }

   outFile.flush();
   if (!outFile.good())
   {
      sErrorMsg = " Error - failed writing temporary FASTA \"" + sOutPath + "\".\n";
      return false;
   }

   return true;
}
