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

#include "CometPlusPrefilterWorkflow.h"
#include "CometPlusApp.h"
#include "CometPlusRuntimeUtils.h"
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <set>
#include <string>
#include <thread>

using namespace std;

void PrepareInputFilesWithPrefilter(vector<InputFileInfo*>& vParsedInputs,
                                    std::vector<InputFileInfo*>& pvInputFiles,
                                    CometInterfaces::ICometSearchManager* pSearchMgr,
                                    const NovelModeOptions& novelOpts,
                                    bool bSetThreadOverride,
                                    int iThreadOverride,
                                    bool bNovelMode,
                                    const std::vector<double>& vNovelMasses,
                                    bool bUseNovelMergedSearch,
                                    const string& sNovelMergedOutputBase,
                                    bool bCreateFragmentIndex,
                                    bool bCreatePeptideIndex,
                                    const std::chrono::steady_clock::time_point& tProgramStart,
                                    std::vector<string>& vTempArtifacts)
{
   bool bNeedScanPrefilter = bNovelMode || novelOpts.HasExplicitScanFilter();
   if (!bNeedScanPrefilter)
   {
      pvInputFiles = vParsedInputs;

      if (bCreateFragmentIndex || bCreatePeptideIndex)
         return;

      char szSummary[1024];
      snprintf(szSummary,
               sizeof(szSummary),
               " [%s] search setup: spectra_inputs=%zu, novel_mode=%s, explicit_scan_filter=%s, novel_mass_count=%zu\n",
               GetLocalTimestampString().c_str(),
               pvInputFiles.size(),
               BoolToOnOffString(bNovelMode),
               BoolToOnOffString(novelOpts.HasExplicitScanFilter()),
               vNovelMasses.size());
      logout(szSummary);
      return;
   }

   set<int> setExplicitScans;
   string sErrorMsg;
   if (!novelOpts.sScanNumbers.empty())
   {
      if (!ParseScanIntegersFromString(novelOpts.sScanNumbers, setExplicitScans, sErrorMsg, "--scan_numbers"))
      {
         logerr(sErrorMsg);
         exit(1);
      }
   }
   if (!novelOpts.sScanFilePath.empty())
   {
      if (!ParseScanIntegersFromFile(novelOpts.sScanFilePath, setExplicitScans, sErrorMsg, "--scan"))
      {
         logerr(sErrorMsg);
         exit(1);
      }
   }
   bool bUseExplicitScans = novelOpts.HasExplicitScanFilter();

   NovelMassFilterContext novelMassCtx = {};
   if (!pSearchMgr->GetParamValue("ms_level", novelMassCtx.iMSLevel))
      novelMassCtx.iMSLevel = 2;
   if (bNovelMode)
   {
      if (!BuildNovelMassFilterContext(pSearchMgr, novelMassCtx, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
   }

   struct PrefilterResult
   {
      string sOriginalPath;
      string sOutputBaseName;
      string sTempMgfPath;
      string sWorkerJobPath;
      string sWorkerResultPath;
      int iNumScansKept = 0;
      double dElapsedSeconds = 0.0;
      bool bSuccess = false;
      string sErrorMsg;
   };

   int iPrefilterThreads = 1;
   if (bSetThreadOverride && iThreadOverride > 0)
   {
      iPrefilterThreads = iThreadOverride;
   }
   else
   {
      int iNumThreadsParam = 0;
      if (pSearchMgr->GetParamValue("num_threads", iNumThreadsParam) && iNumThreadsParam > 0)
         iPrefilterThreads = iNumThreadsParam;
   }

   if (iPrefilterThreads < 1)
      iPrefilterThreads = 1;
   if (iPrefilterThreads > (int)vParsedInputs.size())
      iPrefilterThreads = (int)vParsedInputs.size();
   if (iPrefilterThreads < 1)
      iPrefilterThreads = 1;

   bool bUseProcessPrefilter = (iPrefilterThreads > 1 && HasMzMLbInputs(vParsedInputs));
   string sPrefilterWorkerExePath;
   string sExplicitScansFileForWorker;
   string sNovelMassesFileForWorker;
   string sMassOffsetsFileForWorker;

   if (bUseProcessPrefilter)
   {
      if (!ResolvePrefilterWorkerExecutablePath(g_sCometPlusExecutablePath,
                                                sPrefilterWorkerExePath,
                                                sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }

      if (bUseExplicitScans)
      {
         if (!WriteIntegerSetToTempFile(setExplicitScans,
                                        "cometplus_prefilter_scan_ids",
                                        sExplicitScansFileForWorker,
                                        sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }
      }

      if (bNovelMode)
      {
         if (!WriteDoubleVectorToTempFile(vNovelMasses,
                                          "cometplus_prefilter_novel_masses",
                                          sNovelMassesFileForWorker,
                                          sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }
      }

      if (!novelMassCtx.vMassOffsets.empty())
      {
         if (!WriteDoubleVectorToTempFile(novelMassCtx.vMassOffsets,
                                          "cometplus_prefilter_mass_offsets",
                                          sMassOffsetsFileForWorker,
                                          sErrorMsg))
         {
            logerr(sErrorMsg);
            exit(1);
         }
      }

      char szWarnBuf[1024];
      snprintf(szWarnBuf,
               sizeof(szWarnBuf),
               " [%s] mzMLb input detected: using process-parallel prefilter workers=%d (worker=\"%s\").\n",
               GetLocalTimestampString().c_str(),
               iPrefilterThreads,
               sPrefilterWorkerExePath.c_str());
      logout(szWarnBuf);
   }

   char szPrefilterStart[1024];
   snprintf(szPrefilterStart,
            sizeof(szPrefilterStart),
            " [%s] scan prefilter (%s): %zu input files, worker_threads=%d\n",
            GetLocalTimestampString().c_str(),
            bUseProcessPrefilter ? "process-parallel" : "parallel",
            vParsedInputs.size(),
            iPrefilterThreads);
   logout(szPrefilterStart);

   auto tPrefilterStart = std::chrono::steady_clock::now();
   vector<PrefilterResult> vPrefilterResults(vParsedInputs.size());
   std::atomic<size_t> iNextInput(0);
   std::atomic<bool> bPrefilterError(false);
   std::mutex mError;
   string sFirstPrefilterError;

   vector<std::thread> vWorkers;
   vWorkers.reserve((size_t)iPrefilterThreads);
   for (int iWorker = 0; iWorker < iPrefilterThreads; ++iWorker)
   {
      vWorkers.emplace_back([&]()
      {
         while (true)
         {
            size_t iInput = iNextInput.fetch_add(1);
            if (iInput >= vParsedInputs.size())
               break;

            if (bPrefilterError.load())
               break;

            InputFileInfo* pInputFileInfo = vParsedInputs.at(iInput);
            PrefilterResult result;
            result.sOriginalPath = pInputFileInfo->szFileName;
            result.sOutputBaseName = pInputFileInfo->szBaseName;

            auto tInputPrefilterStart = std::chrono::steady_clock::now();
            string sTempMgfPath;
            int iNumScansKept = 0;
            string sLocalError;
            bool bOk = false;

            if (bUseProcessPrefilter)
            {
               if (!CreateTempPath("cometplus_prefilter_result", ".tsv", result.sWorkerResultPath, sLocalError))
               {
                  bOk = false;
               }
               else if (!WritePrefilterWorkerJobFile(*pInputFileInfo,
                                                     bUseExplicitScans,
                                                     sExplicitScansFileForWorker,
                                                     bNovelMode,
                                                     sNovelMassesFileForWorker,
                                                     novelOpts.sInternalNovelPeptidePath,
                                                     novelMassCtx,
                                                     sMassOffsetsFileForWorker,
                                                     novelOpts.sOutputFolder,
                                                     result.sWorkerResultPath,
                                                     result.sWorkerJobPath,
                                                     sLocalError))
               {
                  bOk = false;
               }
               else if (!RunPrefilterWorkerJob(sPrefilterWorkerExePath, result.sWorkerJobPath, sLocalError))
               {
                  bOk = false;
               }
               else
               {
                  PrefilterWorkerOutput workerOutput;
                  if (!ParsePrefilterWorkerResultFile(result.sWorkerResultPath, workerOutput, sLocalError))
                  {
                     bOk = false;
                  }
                  else if (!workerOutput.bSuccess)
                  {
                     sLocalError = workerOutput.sErrorMsg;
                     if (sLocalError.empty())
                        sLocalError = " Error - prefilter worker reported failure.\n";
                     bOk = false;
                  }
                  else
                  {
                     sTempMgfPath = workerOutput.sTempMgfPath;
                     iNumScansKept = workerOutput.iNumScansKept;
                     bOk = true;
                  }
               }
            }
            else
            {
               bOk = FilterInputFileToTempMgf(*pInputFileInfo,
                                              setExplicitScans,
                                              bUseExplicitScans,
                                              vNovelMasses,
                                              bNovelMode,
                                              novelMassCtx,
                                              sTempMgfPath,
                                              iNumScansKept,
                                              sLocalError);
            }
            result.dElapsedSeconds = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - tInputPrefilterStart).count() / 1000.0;

            if (!bOk)
            {
               result.bSuccess = false;
               result.sErrorMsg = sLocalError;
               bPrefilterError.store(true);
               std::lock_guard<std::mutex> lk(mError);
               if (sFirstPrefilterError.empty())
                  sFirstPrefilterError = sLocalError;
            }
            else
            {
               result.sTempMgfPath = sTempMgfPath;
               result.iNumScansKept = iNumScansKept;
               result.bSuccess = true;
            }

            vPrefilterResults.at(iInput) = result;
         }
      });
   }

   for (size_t i = 0; i < vWorkers.size(); ++i)
      vWorkers.at(i).join();

   if (!sExplicitScansFileForWorker.empty())
      vTempArtifacts.push_back(sExplicitScansFileForWorker);
   if (!sNovelMassesFileForWorker.empty())
      vTempArtifacts.push_back(sNovelMassesFileForWorker);
   if (!sMassOffsetsFileForWorker.empty())
      vTempArtifacts.push_back(sMassOffsetsFileForWorker);
   for (size_t i = 0; i < vPrefilterResults.size(); ++i)
   {
      if (!vPrefilterResults.at(i).sWorkerJobPath.empty())
         vTempArtifacts.push_back(vPrefilterResults.at(i).sWorkerJobPath);
      if (!vPrefilterResults.at(i).sWorkerResultPath.empty())
         vTempArtifacts.push_back(vPrefilterResults.at(i).sWorkerResultPath);
   }

   if (bPrefilterError.load())
   {
      if (sFirstPrefilterError.empty())
      {
         for (size_t i = 0; i < vPrefilterResults.size(); ++i)
         {
            if (!vPrefilterResults.at(i).bSuccess && !vPrefilterResults.at(i).sErrorMsg.empty())
            {
               sFirstPrefilterError = vPrefilterResults.at(i).sErrorMsg;
               break;
            }
         }
      }

      if (sFirstPrefilterError.empty())
         sFirstPrefilterError = " Error - scan prefilter failed.\n";
      logerr(sFirstPrefilterError);
      exit(1);
   }

   for (size_t i = 0; i < vPrefilterResults.size(); ++i)
   {
      vTempArtifacts.push_back(vPrefilterResults.at(i).sTempMgfPath);

      char szLogBuf[1024];
      snprintf(szLogBuf, sizeof(szLogBuf),
               " [%s] scan prefilter: \"%s\" -> %d scans retained (%.3f sec)\n",
               GetLocalTimestampString().c_str(),
               vPrefilterResults.at(i).sOriginalPath.c_str(),
               vPrefilterResults.at(i).iNumScansKept,
               vPrefilterResults.at(i).dElapsedSeconds);
      logout(szLogBuf);
   }

   LogStageTiming("scan prefilter (parallel)", tPrefilterStart, tProgramStart);

   if (bUseNovelMergedSearch)
   {
      auto tMergeStart = std::chrono::steady_clock::now();

      vector<std::pair<string, string>> vMergeInputs;
      vMergeInputs.reserve(vPrefilterResults.size());
      for (size_t i = 0; i < vPrefilterResults.size(); ++i)
      {
         // Keep original TITLE values from filtered MGF shards when basenames are unique.
         vMergeInputs.push_back(std::make_pair(vPrefilterResults.at(i).sTempMgfPath, ""));
      }

      string sMergedMgfPath;
      if (!MergeFilteredMgfFilesWithSourceTag(vMergeInputs, sMergedMgfPath, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
      vTempArtifacts.push_back(sMergedMgfPath);

      InputFileInfo* pMergedInputFile = new InputFileInfo();
      pMergedInputFile->iAnalysisType = AnalysisType_EntireFile;
      pMergedInputFile->iFirstScan = 0;
      pMergedInputFile->iLastScan = 0;
      strncpy(pMergedInputFile->szFileName, sMergedMgfPath.c_str(), SIZE_FILE - 1);
      pMergedInputFile->szFileName[SIZE_FILE - 1] = '\0';
      strncpy(pMergedInputFile->szBaseName, sNovelMergedOutputBase.c_str(), SIZE_FILE - 1);
      pMergedInputFile->szBaseName[SIZE_FILE - 1] = '\0';

      pvInputFiles.push_back(pMergedInputFile);

      for (size_t i = 0; i < vParsedInputs.size(); ++i)
      {
         delete vParsedInputs.at(i);
         vParsedInputs.at(i) = NULL;
      }
      vParsedInputs.clear();

      char szMergeLog[2048];
      snprintf(szMergeLog,
               sizeof(szMergeLog),
               " [%s] filtered MGF merge: %zu files -> \"%s\"\n",
               GetLocalTimestampString().c_str(),
               vMergeInputs.size(),
               sMergedMgfPath.c_str());
      logout(szMergeLog);

      LogStageTiming("filtered MGF merge", tMergeStart, tProgramStart);
   }
   else
   {
      for (size_t i = 0; i < vParsedInputs.size(); ++i)
      {
         InputFileInfo* pInputFileInfo = vParsedInputs.at(i);
         strncpy(pInputFileInfo->szFileName, vPrefilterResults.at(i).sTempMgfPath.c_str(), SIZE_FILE - 1);
         pInputFileInfo->szFileName[SIZE_FILE - 1] = '\0';
         strncpy(pInputFileInfo->szBaseName, vPrefilterResults.at(i).sOutputBaseName.c_str(), SIZE_FILE - 1);
         pInputFileInfo->szBaseName[SIZE_FILE - 1] = '\0';
         pvInputFiles.push_back(pInputFileInfo);
      }
   }

   char szSummary[1024];
   if (bUseNovelMergedSearch)
   {
      snprintf(szSummary,
               sizeof(szSummary),
               " [%s] search setup: spectra_inputs=%zu (merged from %zu source files), novel_mode=%s, explicit_scan_filter=%s, novel_mass_count=%zu\n",
               GetLocalTimestampString().c_str(),
               pvInputFiles.size(),
               vPrefilterResults.size(),
               BoolToOnOffString(bNovelMode),
               BoolToOnOffString(novelOpts.HasExplicitScanFilter()),
               vNovelMasses.size());
   }
   else
   {
      snprintf(szSummary,
               sizeof(szSummary),
               " [%s] search setup: spectra_inputs=%zu, novel_mode=%s, explicit_scan_filter=%s, novel_mass_count=%zu\n",
               GetLocalTimestampString().c_str(),
               pvInputFiles.size(),
               BoolToOnOffString(bNovelMode),
               BoolToOnOffString(novelOpts.HasExplicitScanFilter()),
               vNovelMasses.size());
   }
   logout(szSummary);
}
