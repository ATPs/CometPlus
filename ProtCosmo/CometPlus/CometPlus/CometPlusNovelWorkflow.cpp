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

#include "CometPlusNovelWorkflow.h"
#include "CometPlusApp.h"
#include "CometPlusRuntimeUtils.h"
#include "CometPlusMultiDB.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>

using namespace std;

static string JoinProteinIdsField(const vector<string>& vProteinIds)
{
   string sField;
   for (size_t i = 0; i < vProteinIds.size(); ++i)
   {
      if (i > 0)
         sField += ";";
      sField += vProteinIds.at(i);
   }
   return sField;
}

static bool ParseVariableModMasses(CometInterfaces::ICometSearchManager* pSearchMgr,
                                   vector<double>& vMasses,
                                   vector<bool>& vMassKnown)
{
   vMasses.assign(16, 0.0);
   vMassKnown.assign(16, false);
   if (pSearchMgr == NULL)
      return false;

   for (int iMod = 1; iMod <= 15; ++iMod)
   {
      char szParamName[32];
      snprintf(szParamName, sizeof(szParamName), "variable_mod%02d", iMod);

      VarMods vm;
      if (!pSearchMgr->GetParamValue(szParamName, vm))
         continue;

      vMasses.at((size_t)iMod) = vm.dVarModMass;
      vMassKnown.at((size_t)iMod) = true;
   }

   return true;
}

static string BuildVarModToken(int iModIndex,
                               const vector<double>& vModMasses,
                               const vector<bool>& vModMassKnown)
{
   if (iModIndex <= 0)
      return "";

   if ((size_t)iModIndex < vModMassKnown.size()
         && vModMassKnown.at((size_t)iModIndex))
   {
      char szBuf[64];
      snprintf(szBuf, sizeof(szBuf), "[%+.4f]", vModMasses.at((size_t)iModIndex));
      return szBuf;
   }

   return "[mod" + std::to_string(iModIndex) + "]";
}

static string BuildPeptideWithModMassAnnotation(const PeptideMassEntry& entry,
                                                const vector<double>& vModMasses,
                                                const vector<bool>& vModMassKnown)
{
   string sPeptide = NormalizePeptideToken(entry.sPeptide);
   if (sPeptide.empty())
      return "";

   size_t iLen = sPeptide.length();
   if (entry.vVarModSites.size() < iLen + 2)
      return sPeptide;

   bool bAnyMods = false;
   for (size_t i = 0; i < iLen + 2; ++i)
   {
      if (entry.vVarModSites.at(i) != 0)
      {
         bAnyMods = true;
         break;
      }
   }
   if (!bAnyMods)
      return sPeptide;

   string sOut;
   int iNtermMod = (int)entry.vVarModSites.at(iLen);
   if (iNtermMod > 0)
   {
      sOut += "n";
      sOut += BuildVarModToken(iNtermMod, vModMasses, vModMassKnown);
   }

   for (size_t i = 0; i < iLen; ++i)
   {
      sOut.push_back(sPeptide[i]);
      int iMod = (int)entry.vVarModSites.at(i);
      if (iMod > 0)
         sOut += BuildVarModToken(iMod, vModMasses, vModMassKnown);
   }

   int iCtermMod = (int)entry.vVarModSites.at(iLen + 1);
   if (iCtermMod > 0)
   {
      sOut += "c";
      sOut += BuildVarModToken(iCtermMod, vModMasses, vModMassKnown);
   }

   return sOut;
}

static void ComputeToleranceOnlyMzWindow(double dPepMass,
                                         int iCharge,
                                         const NovelMassFilterContext& ctx,
                                         double& dMz,
                                         double& dMzMin,
                                         double& dMzMax)
{
   dMz = (dPepMass + (iCharge - 1) * PROTON_MASS) / iCharge;

   if (ctx.iTolUnits == 2) // ppm
   {
      dMzMin = dMz + dMz * ctx.dTolLower / 1.0E6;
      dMzMax = dMz + dMz * ctx.dTolUpper / 1.0E6;
   }
   else
   {
      double dTolLow = ctx.dTolLower;
      double dTolHigh = ctx.dTolUpper;

      if (ctx.iTolUnits == 1) // mmu -> amu
      {
         dTolLow *= 0.001;
         dTolHigh *= 0.001;
      }

      if (ctx.iTolType == 1)
      {
         dTolLow *= iCharge;
         dTolHigh *= iCharge;
      }

      double dExpMassMin = dPepMass - dTolHigh;
      double dExpMassMax = dPepMass - dTolLow;
      dMzMin = (dExpMassMin + (iCharge - 1) * PROTON_MASS) / iCharge;
      dMzMax = (dExpMassMax + (iCharge - 1) * PROTON_MASS) / iCharge;
   }

   if (dMzMin > dMzMax)
      std::swap(dMzMin, dMzMax);
}

void RunNovelWorkflowIfNeeded(const NovelModeOptions& novelOpts,
                              CometInterfaces::ICometSearchManager* pSearchMgr,
                              const std::vector<string>& vKnownDatabases,
                              bool bKnownAllIdx,
                              int iKnownIdxType,
                              const string& sParamsFile,
                              int iThreadOverride,
                              bool bTreatSameIL,
                              int iDecoySearch,
                              const string& sResolvedOutputInternalNovelPath,
                              const std::chrono::steady_clock::time_point& tProgramStart,
                              std::vector<string>& vTempArtifacts,
                              string& sMergedDatabasePath,
                              std::vector<NovelPeptideRecord>& vNovelRecords,
                              std::vector<double>& vNovelMasses,
                              string& sNovelScoringDatabase)
{
   if (!novelOpts.HasNovelMode())
      return;

   if (vKnownDatabases.empty())
   {
      string strErrorMsg = " Error - known database must be provided through --database or database_name in params when using novel options.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   auto AppendUniqueProtein = [](vector<string>& vProteinIds, const string& sProteinId)
   {
      if (sProteinId.empty())
         return;
      for (size_t i = 0; i < vProteinIds.size(); ++i)
      {
         if (vProteinIds.at(i) == sProteinId)
            return;
      }
      vProteinIds.push_back(sProteinId);
   };

   auto BuildNoVarModOverrides = []() -> map<string, string>
   {
      map<string, string> mOverrides;
      mOverrides["max_variable_mods_in_peptide"] = "0";
      mOverrides["require_variable_mod"] = "0";
      for (int iMod = 1; iMod <= 15; ++iMod)
      {
         char szParamName[32];
         snprintf(szParamName, sizeof(szParamName), "variable_mod%02d", iMod);
         mOverrides[szParamName] = "0.0 X 0 3 -1 0 0 0.0";
      }
      return mOverrides;
   };

   if (novelOpts.HasNovelInputOptions())
   {
      string sTmpNoVarModParams;
      {
         string sErrorMsg;
         map<string, string> mNoVarOverrides = BuildNoVarModOverrides();
         if (!BuildTemporaryParamsFile(sParamsFile, mNoVarOverrides, sTmpNoVarModParams, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }
         vTempArtifacts.push_back(sTmpNoVarModParams);
      }

      set<string> setKnownPeptides;
      {
         auto tStageStart = std::chrono::steady_clock::now();
         vector<PeptideMassEntry> vKnownEntries;
         string sErrorMsg;

         if (bKnownAllIdx)
         {
            for (size_t iDb = 0; iDb < vKnownDatabases.size(); ++iDb)
            {
               const string& sDb = vKnownDatabases.at(iDb);

               int iType = 0;
               string sProbeError;
               if (!CometPlusProbeIdxType(sDb, iType, sProbeError))
               {
                  logerr(sProbeError);
                  exit(1);
               }

               vector<PeptideMassEntry> vEntries;
               bool bOk = false;
               if (iType == 1)
                  bOk = ParseFragmentIdxEntries(sDb, vEntries, sErrorMsg);
               else if (iType == 2)
                  bOk = ParsePeptideIdxEntries(sDb, vEntries, false, sErrorMsg);
               else
               {
                  sErrorMsg = " Error - unknown .idx type encountered while parsing known database.\n";
                  bOk = false;
               }

               if (!bOk)
               {
                  logerr(sErrorMsg);
                  exit(1);
               }

               vKnownEntries.insert(vKnownEntries.end(), vEntries.begin(), vEntries.end());
            }
         }
         else
         {
            string sTmpKnownFasta;
            if (!BuildMergedFasta(vKnownDatabases, sTmpKnownFasta, sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }
            vTempArtifacts.push_back(sTmpKnownFasta);

            if (!RunCometForIndexGeneration(g_sCometPlusExecutablePath,
                                            sTmpNoVarModParams,
                                            sTmpKnownFasta,
                                            false,
                                            iThreadOverride,
                                            sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }

            string sTmpKnownIdx = sTmpKnownFasta + ".idx";
            vTempArtifacts.push_back(sTmpKnownIdx);
            if (!ParsePeptideIdxEntries(sTmpKnownIdx, vKnownEntries, false, sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }
         }

         for (size_t i = 0; i < vKnownEntries.size(); ++i)
         {
            string sKey = NormalizePeptideForCompare(vKnownEntries.at(i).sPeptide, bTreatSameIL);
            if (!sKey.empty())
               setKnownPeptides.insert(sKey);
         }

         LogStageTiming("known peptide extraction", tStageStart, tProgramStart);
      }

      struct NovelCandidateAggregate
      {
         string sRepresentativePeptide;
         vector<string> vProteinIds;
      };

      unordered_map<string, NovelCandidateAggregate> mNovelCandidates;
      {
         auto tStageStart = std::chrono::steady_clock::now();
         if (!novelOpts.sNovelProteinPath.empty())
         {
            string sErrorMsg;
            string sTmpNovelProtein;
            if (!BuildMergedFasta(vector<string>(1, novelOpts.sNovelProteinPath), sTmpNovelProtein, sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }
            vTempArtifacts.push_back(sTmpNovelProtein);

            if (!RunCometForIndexGeneration(g_sCometPlusExecutablePath,
                                            sTmpNoVarModParams,
                                            sTmpNovelProtein,
                                            false,
                                            iThreadOverride,
                                            sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }
            string sTmpNovelProteinIdx = sTmpNovelProtein + ".idx";
            vTempArtifacts.push_back(sTmpNovelProteinIdx);

            vector<PeptideMassEntry> vTmp;
            if (!ParsePeptideIdxEntries(sTmpNovelProteinIdx, vTmp, true, sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }

            for (size_t i = 0; i < vTmp.size(); ++i)
            {
               string sNormPep = NormalizePeptideToken(vTmp.at(i).sPeptide);
               string sKey = NormalizePeptideForCompare(sNormPep, bTreatSameIL);
               if (sKey.empty())
                  continue;

               auto it = mNovelCandidates.find(sKey);
               if (it == mNovelCandidates.end())
               {
                  NovelCandidateAggregate agg;
                  agg.sRepresentativePeptide = sNormPep;
                  mNovelCandidates[sKey] = agg;
                  it = mNovelCandidates.find(sKey);
               }

               for (size_t iProt = 0; iProt < vTmp.at(i).vProteinIds.size(); ++iProt)
                  AppendUniqueProtein(it->second.vProteinIds, vTmp.at(i).vProteinIds.at(iProt));
            }
         }

         if (!novelOpts.sNovelPeptidePath.empty())
         {
            string sErrorMsg;
            vector<string> vInputPeptides;
            if (!ParseNovelPeptideFile(novelOpts.sNovelPeptidePath, vInputPeptides, sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }

            if (vInputPeptides.empty())
            {
               string strErrorMsg = " Error - no peptide entries were parsed from --novel_peptide input.\n";
               logerr(strErrorMsg);
               exit(1);
            }

            for (size_t i = 0; i < vInputPeptides.size(); ++i)
            {
               string sNormPep = NormalizePeptideToken(vInputPeptides.at(i));
               string sKey = NormalizePeptideForCompare(sNormPep, bTreatSameIL);
               if (sKey.empty())
                  continue;

               auto it = mNovelCandidates.find(sKey);
               if (it == mNovelCandidates.end())
               {
                  NovelCandidateAggregate agg;
                  agg.sRepresentativePeptide = sNormPep;
                  mNovelCandidates[sKey] = agg;
                  it = mNovelCandidates.find(sKey);
               }

               AppendUniqueProtein(it->second.vProteinIds, sNormPep);
            }
         }

         LogStageTiming("novel candidate assembly", tStageStart, tProgramStart);
      }

      {
         auto tStageStart = std::chrono::steady_clock::now();
         int iSubtractedCount = 0;
         vector<string> vKeys;
         vKeys.reserve(mNovelCandidates.size());
         for (auto it = mNovelCandidates.begin(); it != mNovelCandidates.end(); ++it)
            vKeys.push_back(it->first);
         std::sort(vKeys.begin(), vKeys.end());

         vNovelRecords.clear();
         for (size_t i = 0; i < vKeys.size(); ++i)
         {
            if (setKnownPeptides.find(vKeys.at(i)) != setKnownPeptides.end())
            {
               iSubtractedCount++;
               continue;
            }

            NovelPeptideRecord rec;
            rec.sPeptide = mNovelCandidates[vKeys.at(i)].sRepresentativePeptide;
            rec.vProteinIds = mNovelCandidates[vKeys.at(i)].vProteinIds;
            vNovelRecords.push_back(rec);
         }

         std::sort(vNovelRecords.begin(),
                   vNovelRecords.end(),
                   [](const NovelPeptideRecord& a, const NovelPeptideRecord& b)
                   {
                      return a.sPeptide < b.sPeptide;
                   });

         for (size_t i = 0; i < vNovelRecords.size(); ++i)
            vNovelRecords.at(i).sPeptideId = "COMETPLUS_NOVEL_" + std::to_string(i + 1);

         char szLogBuf[512];
         snprintf(szLogBuf, sizeof(szLogBuf),
                  " Novel mode: %zu unique novel peptides parsed, %d removed by known-db subtraction, %zu retained.\n",
                  mNovelCandidates.size(),
                  iSubtractedCount,
                  vNovelRecords.size());
         logout(szLogBuf);

         LogStageTiming("novel candidate subtraction", tStageStart, tProgramStart);
      }
   }
   vector<double> vImportedPrecomputedMasses;
   bool bImportedHasDetailedMzColumns = false;
   if (novelOpts.HasInternalNovelInput())
   {
      auto tStageStart = std::chrono::steady_clock::now();
      string sErrorMsg;
      if (!ParseInternalNovelPeptideFile(novelOpts.sInternalNovelPeptidePath,
                                         vNovelRecords,
                                         sErrorMsg,
                                         &vImportedPrecomputedMasses,
                                         &bImportedHasDetailedMzColumns))
      {
         logerr(sErrorMsg);
         exit(1);
      }

      char szLogBuf[512];
      snprintf(szLogBuf, sizeof(szLogBuf),
               " Novel mode: imported %zu internal novel peptide records from \"%s\".\n",
               vNovelRecords.size(),
               novelOpts.sInternalNovelPeptidePath.c_str());
      logout(szLogBuf);

      if (bImportedHasDetailedMzColumns && !vImportedPrecomputedMasses.empty())
      {
         char szMassBuf[512];
         snprintf(szMassBuf,
                  sizeof(szMassBuf),
                  " Novel mode: detailed internal TSV detected; imported %zu precomputed novel masses.\n",
                  vImportedPrecomputedMasses.size());
         logout(szMassBuf);
      }

      LogStageTiming("internal novel peptide import", tStageStart, tProgramStart);
   }

   if (iDecoySearch == 0)
   {
      logout(" Warning - decoy_search=0 in novel mode; novel peptide inputs are treated as target entries unless already encoded as decoy sequences.\n");
   }

   vector<PeptideMassEntry> vNovelMassEntriesForExport;
   if (!vNovelRecords.empty())
   {
      auto tStageStart = std::chrono::steady_clock::now();
      bool bDisableNtermMethionineClipForNovel = (!novelOpts.sNovelPeptidePath.empty()
            || novelOpts.HasInternalNovelInput());
      string sNovelFastaPath;
      string sErrorMsg;
      if (!WriteNovelRecordsToFasta(vNovelRecords, "cometplus_novel_scoring", sNovelFastaPath, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
      vTempArtifacts.push_back(sNovelFastaPath);

      bool bUseImportedDetailedMasses = false;
      if (novelOpts.HasInternalNovelInput()
            && bImportedHasDetailedMzColumns
            && !vImportedPrecomputedMasses.empty())
      {
         vNovelMasses = vImportedPrecomputedMasses;
         bUseImportedDetailedMasses = true;

         char szMassLog[512];
         snprintf(szMassLog,
                  sizeof(szMassLog),
                  " Novel mode: using %zu precomputed novel masses from detailed internal TSV (skip temporary mass-index generation).\n",
                  vNovelMasses.size());
         logout(szMassLog);
      }

      if (!bUseImportedDetailedMasses)
      {
         int iNoCutEnzyme = -1;
         if (!FindNoCutEnzymeNumber(sParamsFile, iNoCutEnzyme, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }

         map<string, string> mNoCutOverrides;
         mNoCutOverrides["search_enzyme_number"] = std::to_string(iNoCutEnzyme);
         mNoCutOverrides["search_enzyme2_number"] = "0";
         mNoCutOverrides["allowed_missed_cleavage"] = "0";
         mNoCutOverrides["num_enzyme_termini"] = "2";
         if (bDisableNtermMethionineClipForNovel)
            mNoCutOverrides["clip_nterm_methionine"] = "0";

         string sTmpNoCutParams;
         if (!BuildTemporaryParamsFile(sParamsFile, mNoCutOverrides, sTmpNoCutParams, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }
         vTempArtifacts.push_back(sTmpNoCutParams);

         if (!RunCometForIndexGeneration(g_sCometPlusExecutablePath,
                                         sTmpNoCutParams,
                                         sNovelFastaPath,
                                         false,
                                         iThreadOverride,
                                         sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }

         string sTmpNovelMassIdx = sNovelFastaPath + ".idx";
         vTempArtifacts.push_back(sTmpNovelMassIdx);

         if (!ParsePeptideIdxEntries(sTmpNovelMassIdx, vNovelMassEntriesForExport, false, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }

         set<double> setNovelMasses;
         for (size_t i = 0; i < vNovelMassEntriesForExport.size(); ++i)
         {
            if (vNovelMassEntriesForExport.at(i).dMass > 0.0)
               setNovelMasses.insert(vNovelMassEntriesForExport.at(i).dMass);
         }
         vNovelMasses.assign(setNovelMasses.begin(), setNovelMasses.end());
      }

      if (!(novelOpts.bStopAfterSavingNovelPeptide && !sResolvedOutputInternalNovelPath.empty()))
      {
         if (bKnownAllIdx)
         {
            bool bFragment = (iKnownIdxType == 1);
            string sScoringParamsPath = sParamsFile;
            if (bDisableNtermMethionineClipForNovel)
            {
               map<string, string> mScoringOverrides;
               mScoringOverrides["clip_nterm_methionine"] = "0";
               if (!BuildTemporaryParamsFile(sParamsFile, mScoringOverrides, sScoringParamsPath, sErrorMsg))
               {
                  logerr(sErrorMsg);
                  exit(1);
               }
               vTempArtifacts.push_back(sScoringParamsPath);
            }

            if (!RunCometForIndexGeneration(g_sCometPlusExecutablePath,
                                            sScoringParamsPath,
                                            sNovelFastaPath,
                                            bFragment,
                                            iThreadOverride,
                                            sErrorMsg))
            {
               logerr(sErrorMsg);
               exit(1);
            }

            sNovelScoringDatabase = sNovelFastaPath + ".idx";
            vTempArtifacts.push_back(sNovelScoringDatabase);
         }
         else
         {
            sNovelScoringDatabase = sNovelFastaPath;
         }
      }

      LogStageTiming("novel mass calculation", tStageStart, tProgramStart);
   }
   else
   {
      logout(" Warning - no novel peptides remain after subtraction against known database(s); spectrum filtering will likely remove all scans.\n");
   }

   if (!sResolvedOutputInternalNovelPath.empty())
   {
      if (pSearchMgr == NULL)
      {
         logerr(" Error - internal novel peptide export requires initialized search manager.\n");
         exit(1);
      }

      auto tStageStart = std::chrono::steady_clock::now();
      string sOutputDir = GetDirNameLocal(sResolvedOutputInternalNovelPath);
      if (!sOutputDir.empty())
      {
         string sErrorMsg;
         if (!EnsureDirectoryExistsRecursive(sOutputDir, sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }
      }

      std::ofstream outFile(sResolvedOutputInternalNovelPath.c_str(),
                            std::ios::out | std::ios::trunc);
      if (!outFile.good())
      {
         string strErrorMsg = " Error - cannot create internal novel peptide file \"" + sResolvedOutputInternalNovelPath + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }

      NovelMassFilterContext exportCtx = {};
      string sErrorMsg;
      if (!BuildNovelMassFilterContext(pSearchMgr, exportCtx, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }

      int iMinCharge = exportCtx.iMinPrecursorCharge;
      int iMaxCharge = exportCtx.iMaxPrecursorCharge;
      if (iMinCharge < 1)
         iMinCharge = 1;
      if (iMaxCharge < iMinCharge)
         iMaxCharge = iMinCharge;

      vector<double> vVarModMasses;
      vector<bool> vVarModMassKnown;
      if (!ParseVariableModMasses(pSearchMgr, vVarModMasses, vVarModMassKnown))
      {
         logerr(" Error - failed to parse variable mod masses for internal novel peptide export.\n");
         exit(1);
      }

      unordered_map<string, size_t> mPeptideToRecord;
      unordered_map<string, size_t> mMethionineClippedToRecord;
      set<string> setAmbiguousMethionineClipped;
      for (size_t i = 0; i < vNovelRecords.size(); ++i)
      {
         string sPeptide = NormalizePeptideToken(vNovelRecords.at(i).sPeptide);
         if (!sPeptide.empty())
         {
            mPeptideToRecord[sPeptide] = i;

            if (sPeptide.length() > 1 && sPeptide[0] == 'M')
            {
               string sClipped = sPeptide.substr(1);
               auto itClipped = mMethionineClippedToRecord.find(sClipped);
               if (itClipped == mMethionineClippedToRecord.end())
               {
                  if (setAmbiguousMethionineClipped.find(sClipped) == setAmbiguousMethionineClipped.end())
                     mMethionineClippedToRecord[sClipped] = i;
               }
               else if (itClipped->second != i)
               {
                  mMethionineClippedToRecord.erase(itClipped);
                  setAmbiguousMethionineClipped.insert(sClipped);
               }
            }
         }
      }

      outFile << "peptide\tpeptide_id\tprotein_id\tpeptide_with_mod\tcharge\tmz\tmz_window_min\tmz_window_max\n";

      size_t tRowCount = 0;
      size_t tMappedByMethionineClip = 0;
      size_t tUnmappedEntries = 0;
      unordered_map<string, string> mSyntheticPeptideIds;
      size_t tNextSyntheticId = 1;
      for (size_t i = 0; i < vNovelMassEntriesForExport.size(); ++i)
      {
         const PeptideMassEntry& entry = vNovelMassEntriesForExport.at(i);
         if (entry.dMass <= 0.0)
            continue;

         string sPeptide = NormalizePeptideToken(entry.sPeptide);
         if (sPeptide.empty())
            continue;

         const NovelPeptideRecord* pRecord = NULL;
         auto itRecord = mPeptideToRecord.find(sPeptide);
         if (itRecord != mPeptideToRecord.end())
         {
            pRecord = &vNovelRecords.at(itRecord->second);
         }
         else
         {
            auto itClipped = mMethionineClippedToRecord.find(sPeptide);
            if (itClipped != mMethionineClippedToRecord.end()
                  && setAmbiguousMethionineClipped.find(sPeptide) == setAmbiguousMethionineClipped.end())
            {
               pRecord = &vNovelRecords.at(itClipped->second);
               ++tMappedByMethionineClip;
            }
         }

         string sPeptideWithMod = BuildPeptideWithModMassAnnotation(entry, vVarModMasses, vVarModMassKnown);
         if (sPeptideWithMod.empty())
            sPeptideWithMod = sPeptide;

         string sPeptideIdField;
         string sProteinIdField;
         if (pRecord != NULL)
         {
            sPeptideIdField = pRecord->sPeptideId;
            sProteinIdField = JoinProteinIdsField(pRecord->vProteinIds);
         }
         else
         {
            ++tUnmappedEntries;

            auto itSynthetic = mSyntheticPeptideIds.find(sPeptide);
            if (itSynthetic == mSyntheticPeptideIds.end())
            {
               string sSyntheticId = "COMETPLUS_NOVEL_UNMAPPED_" + std::to_string(tNextSyntheticId);
               tNextSyntheticId++;
               mSyntheticPeptideIds[sPeptide] = sSyntheticId;
               sPeptideIdField = sSyntheticId;
            }
            else
            {
               sPeptideIdField = itSynthetic->second;
            }

            sProteinIdField = sPeptide;
         }

         for (int z = iMinCharge; z <= iMaxCharge; ++z)
         {
            double dMz = 0.0;
            double dMzMin = 0.0;
            double dMzMax = 0.0;
            ComputeToleranceOnlyMzWindow(entry.dMass, z, exportCtx, dMz, dMzMin, dMzMax);

            char szMz[64];
            char szMzMin[64];
            char szMzMax[64];
            snprintf(szMz, sizeof(szMz), "%.17g", dMz);
            snprintf(szMzMin, sizeof(szMzMin), "%.17g", dMzMin);
            snprintf(szMzMax, sizeof(szMzMax), "%.17g", dMzMax);

            outFile << sPeptide
                    << "\t" << sPeptideIdField
                    << "\t" << sProteinIdField
                    << "\t" << sPeptideWithMod
                    << "\t" << z
                    << "\t" << szMz
                    << "\t" << szMzMin
                    << "\t" << szMzMax
                    << "\n";
            ++tRowCount;
         }
      }

      outFile.flush();
      if (!outFile.good())
      {
         string strErrorMsg = " Error - failed writing internal novel peptide file \"" + sResolvedOutputInternalNovelPath + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }

      char szLogBuf[1024];
      snprintf(szLogBuf,
               sizeof(szLogBuf),
               " Saved internal novel peptide file: %s (%zu rows)\n",
               sResolvedOutputInternalNovelPath.c_str(),
               tRowCount);
      logout(szLogBuf);

      if (tMappedByMethionineClip > 0)
      {
         char szMapBuf[512];
         snprintf(szMapBuf,
                  sizeof(szMapBuf),
                  " Note - internal novel peptide export mapped %zu mass entries via N-term methionine clipping.\n",
                  tMappedByMethionineClip);
         logout(szMapBuf);
      }

      if (tUnmappedEntries > 0)
      {
         char szWarnBuf[512];
         snprintf(szWarnBuf,
                  sizeof(szWarnBuf),
                  " Warning - internal novel peptide export could not map %zu mass entries to source peptides; synthetic peptide_id/protein_id were used.\n",
                  tUnmappedEntries);
         logout(szWarnBuf);
      }
      LogStageTiming("internal novel peptide export", tStageStart, tProgramStart);

      if (novelOpts.bStopAfterSavingNovelPeptide)
      {
         char szStopBuf[1024];
         snprintf(szStopBuf,
                  sizeof(szStopBuf),
                  " [%s] stop-after-saving-novel-peptide enabled; exiting without spectrum prefilter/search.\n",
                  GetLocalTimestampString().c_str());
         logout(szStopBuf);

         if (g_bCometPlusKeepTempFiles)
         {
            LogRetainedTempArtifacts(vTempArtifacts, sMergedDatabasePath);
         }
         else
         {
            for (size_t i = 0; i < vTempArtifacts.size(); ++i)
            {
               if (!vTempArtifacts.at(i).empty())
                  remove(vTempArtifacts.at(i).c_str());
            }
            if (!sMergedDatabasePath.empty())
               remove(sMergedDatabasePath.c_str());
         }
         exit(0);
      }
   }
}
