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
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <unordered_map>

using namespace std;

void RunNovelWorkflowIfNeeded(const NovelModeOptions& novelOpts,
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
   else if (novelOpts.HasInternalNovelInput())
   {
      auto tStageStart = std::chrono::steady_clock::now();
      string sErrorMsg;
      if (!ParseInternalNovelPeptideFile(novelOpts.sInternalNovelPeptidePath, vNovelRecords, sErrorMsg))
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
      LogStageTiming("internal novel peptide import", tStageStart, tProgramStart);
   }

   if (!sResolvedOutputInternalNovelPath.empty())
   {
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

      outFile << "peptide\tpeptide_id\tprotein_id\n";
      int iAutoId = 1;
      for (size_t i = 0; i < vNovelRecords.size(); ++i)
      {
         string sPeptide = NormalizePeptideToken(vNovelRecords.at(i).sPeptide);
         if (sPeptide.empty())
            continue;

         string sPeptideId = vNovelRecords.at(i).sPeptideId;
         if (sPeptideId.empty())
         {
            sPeptideId = "COMETPLUS_NOVEL_" + std::to_string(iAutoId);
            iAutoId++;
         }

         string sProteinIdField;
         for (size_t j = 0; j < vNovelRecords.at(i).vProteinIds.size(); ++j)
         {
            if (j > 0)
               sProteinIdField += ";";
            sProteinIdField += vNovelRecords.at(i).vProteinIds.at(j);
         }

         outFile << sPeptide << "\t" << sPeptideId << "\t" << sProteinIdField << "\n";
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
               vNovelRecords.size());
      logout(szLogBuf);
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

   if (iDecoySearch == 0)
   {
      logout(" Warning - decoy_search=0 in novel mode; novel peptide inputs are treated as target entries unless already encoded as decoy sequences.\n");
   }

   if (!vNovelRecords.empty())
   {
      auto tStageStart = std::chrono::steady_clock::now();
      string sNovelFastaPath;
      string sErrorMsg;
      if (!WriteNovelRecordsToFasta(vNovelRecords, "cometplus_novel_scoring", sNovelFastaPath, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
      vTempArtifacts.push_back(sNovelFastaPath);

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

      vector<PeptideMassEntry> vNovelMassEntries;
      if (!ParsePeptideIdxEntries(sTmpNovelMassIdx, vNovelMassEntries, false, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }

      set<double> setNovelMasses;
      for (size_t i = 0; i < vNovelMassEntries.size(); ++i)
      {
         if (vNovelMassEntries.at(i).dMass > 0.0)
            setNovelMasses.insert(vNovelMassEntries.at(i).dMass);
      }
      vNovelMasses.assign(setNovelMasses.begin(), setNovelMasses.end());

      if (bKnownAllIdx)
      {
         bool bFragment = (iKnownIdxType == 1);
         if (!RunCometForIndexGeneration(g_sCometPlusExecutablePath,
                                         sParamsFile,
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

      LogStageTiming("novel mass calculation", tStageStart, tProgramStart);
   }
   else
   {
      logout(" Warning - no novel peptides remain after subtraction against known database(s); spectrum filtering will likely remove all scans.\n");
   }
}
