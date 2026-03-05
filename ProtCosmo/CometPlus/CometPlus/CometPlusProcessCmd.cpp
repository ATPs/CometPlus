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

#include "CometPlusApp.h"
#include "CometPlusMultiDB.h"
#include "CometPlusParams.h"
#include "CometPlusNovelWorkflow.h"
#include "CometPlusInputFiles.h"
#include "CometPlusPrefilterWorkflow.h"
#include "CometPlusRuntimeUtils.h"
#include "NovelModeUtils.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <mutex>
#include <set>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace CometInterfaces;
using namespace std;

void ProcessCmdLine(int argc,
                    char *argv[],
                    char *szParamsFile,
                    vector<InputFileInfo*> &pvInputFiles,
                    ICometSearchManager *pSearchMgr,
                    string &sMergedDatabasePath,
                    vector<string> &vTempArtifacts,
                    bool &bSearchAlreadyRun,
                    bool &bSearchSucceeded)
{
   const auto tProgramStart = std::chrono::steady_clock::now();
   bSearchAlreadyRun = false;
   bSearchSucceeded = false;

   if (argc == 1)
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + " Error - no input files specified so nothing to do.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   string sParamsFile = "comet.params";
   int iPrintParams = 0;
   bool bCreateFragmentIndex = false;
   bool bCreatePeptideIndex = false;
   bool bSetOutputName = false;
   bool bSetFirstScan = false;
   bool bSetLastScan = false;
   bool bSetBatchSize = false;
   bool bSetThreadOverride = false;
   bool bHelpFull = false;
   bool bOutputFolderSpecified = false;
   bool bRunCometEachRequested = false;
   bool bForceNovelOutputOnly = false;
   bool bInputFilesOptionSeen = false;
   int iFirstScan = 0;
   int iLastScan = 0;
   int iBatchSize = 0;
   int iThreadOverride = 0;
   string sOutputName;
   string sInputFilesPath;
   string sResolvedOutputInternalNovelPath;
   vector<string> vDatabases;
   vector<string> vInputArgs;
   vector<CmdParamOverride> vCliParamOverrides;
   NovelModeOptions novelOpts;
   g_bCometPlusKeepTempFiles = false;

   for (int iArg = 1; iArg < argc; ++iArg)
   {
      string sArg = argv[iArg];

      if (sArg == "--help" || sArg == "-h")
      {
         Usage(argv[0]);
      }

      if (sArg.rfind("--", 0) == 0)
      {
         string sName;
         string sValue;
         size_t iPos = sArg.find('=');
         if (iPos == string::npos)
         {
            sName = sArg.substr(2);
         }
         else
         {
            sName = sArg.substr(2, iPos - 2);
            sValue = sArg.substr(iPos + 1);
         }

         if (sName.empty())
         {
            string strErrorMsg = " Error - unknown option \"" + sArg + "\".\n";
            logerr(strErrorMsg);
            exit(1);
         }

         auto RequireValue = [&](const string& sOptName) -> string
         {
            if (!sValue.empty())
               return sValue;

            if (iArg + 1 >= argc)
            {
               string strErrorMsg = " Error - missing value for option --" + sOptName + ".\n";
               logerr(strErrorMsg);
               exit(1);
            }
            ++iArg;
            return argv[iArg];
         };

         if (sName == "database")
         {
            vDatabases.push_back(RequireValue(sName));
         }
         else if (sName == "help-full")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --help-full does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            bHelpFull = true;
         }
         else if (sName == "thread")
         {
            iThreadOverride = atoi(RequireValue(sName).c_str());
            bSetThreadOverride = true;
         }
         else if (sName == "params")
         {
            sParamsFile = RequireValue(sName);
         }
         else if (sName == "input_files")
         {
            if (bInputFilesOptionSeen)
            {
               string strErrorMsg = " Error - option --input_files cannot be provided more than once.\n";
               logerr(strErrorMsg);
               exit(1);
            }

            sInputFilesPath = RequireValue(sName);
            bInputFilesOptionSeen = true;
         }
         else if (sName == "name")
         {
            sOutputName = RequireValue(sName);
            bSetOutputName = true;
         }
         else if (sName == "output-folder")
         {
            novelOpts.sOutputFolder = RequireValue(sName);
            bOutputFolderSpecified = true;
         }
         else if (sName == "first-scan")
         {
            iFirstScan = atoi(RequireValue(sName).c_str());
            bSetFirstScan = true;
         }
         else if (sName == "last-scan")
         {
            iLastScan = atoi(RequireValue(sName).c_str());
            bSetLastScan = true;
         }
         else if (sName == "novel_protein")
         {
            novelOpts.sNovelProteinPath = RequireValue(sName);
         }
         else if (sName == "novel_peptide")
         {
            novelOpts.sNovelPeptidePath = RequireValue(sName);
         }
         else if (sName == "output_internal_novel_peptide")
         {
            novelOpts.sOutputInternalNovelPeptidePath = RequireValue(sName);
         }
         else if (sName == "internal_novel_peptide")
         {
            novelOpts.sInternalNovelPeptidePath = RequireValue(sName);
         }
         else if (sName == "stop-after-saving-novel-peptide")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --stop-after-saving-novel-peptide does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            novelOpts.bStopAfterSavingNovelPeptide = true;
         }
         else if (sName == "keep-tmp")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --keep-tmp does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            g_bCometPlusKeepTempFiles = true;
         }
         else if (sName == "scan")
         {
            novelOpts.sScanFilePath = RequireValue(sName);
         }
         else if (sName == "scan_numbers")
         {
            novelOpts.sScanNumbers = RequireValue(sName);
         }
         else if (sName == "run-comet-each")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --run-comet-each does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            bRunCometEachRequested = true;
         }
         else if (sName == "cometplus-novel-output-only")
         {
            if (!sValue.empty())
            {
               string strErrorMsg = " Error - option --cometplus-novel-output-only does not take a value.\n";
               logerr(strErrorMsg);
               exit(1);
            }
            bForceNovelOutputOnly = true;
         }
         else
         {
            CmdParamOverride cliOverride;
            cliOverride.sName = sName;
            cliOverride.sValue = RequireValue(sName);
            vCliParamOverrides.push_back(cliOverride);
         }

         continue;
      }

      if (sArg.length() > 1 && sArg[0] == '-')
      {
         auto RequireShortValue = [&](char cOpt) -> string
         {
            if (sArg.length() > 2)
               return sArg.substr(2);

            if (iArg + 1 >= argc)
            {
               string strErrorMsg = " Error - missing value for option -" + string(1, cOpt) + ".\n";
               logerr(strErrorMsg);
               exit(1);
            }
            ++iArg;
            return argv[iArg];
         };

         switch (sArg[1])
         {
            case 'p':
               iPrintParams = 1;
               break;
            case 'q':
               iPrintParams = 2;
               break;
            case 'P':
               sParamsFile = RequireShortValue('P');
               break;
            case 'D':
               vDatabases.push_back(RequireShortValue('D'));
               break;
            case 'N':
               sOutputName = RequireShortValue('N');
               bSetOutputName = true;
               break;
            case 'F':
               iFirstScan = atoi(RequireShortValue('F').c_str());
               bSetFirstScan = true;
               break;
            case 'L':
               iLastScan = atoi(RequireShortValue('L').c_str());
               bSetLastScan = true;
               break;
            case 'B':
               iBatchSize = atoi(RequireShortValue('B').c_str());
               bSetBatchSize = true;
               break;
            case 'i':
               bCreateFragmentIndex = true;
               bCreatePeptideIndex = false;
               break;
            case 'j':
               bCreateFragmentIndex = false;
               bCreatePeptideIndex = true;
               break;
            default:
            {
               string strErrorMsg = " Error - unknown option \"" + sArg + "\".\n";
               logerr(strErrorMsg);
               exit(1);
            }
         }

         continue;
      }

      vInputArgs.push_back(sArg);
   }

   if (bHelpFull)
   {
      Usage(argv[0], true, sParamsFile.c_str());
   }

   if (iPrintParams)
   {
      PrintParams(iPrintParams);
      exit(0);
   }

   string sInputFilesError;
   if (!ResolveInputArgsFromInputFilesOption(bInputFilesOptionSeen,
                                             sInputFilesPath,
                                             vInputArgs,
                                             sInputFilesError))
   {
      logerr(sInputFilesError);
      exit(1);
   }

   strncpy(szParamsFile, sParamsFile.c_str(), SIZE_FILE - 1);
   szParamsFile[SIZE_FILE - 1] = '\0';

   set<string> setParamKeys;
   string sKeyErrMsg;
   if (!CollectParamsFileKeys(szParamsFile, setParamKeys, sKeyErrMsg))
   {
      string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
         + sKeyErrMsg;
      logerr(strErrorMsg);
      exit(1);
   }

   auto DedicatedHint = [](const string &sParamKey) -> string
   {
      if (sParamKey == "database_name")
         return "--database";
      if (sParamKey == "num_threads")
         return "--thread";
      if (sParamKey == "scan_range")
         return "--first-scan/--last-scan";
      if (sParamKey == "spectrum_batch_size")
         return "-B";
      return "";
   };

   const set<string> &setDedicatedKeys = GetDedicatedOverrideKeys();
   for (size_t i = 0; i < vCliParamOverrides.size(); ++i)
   {
      const string &sOverrideName = vCliParamOverrides.at(i).sName;

      if (setDedicatedKeys.find(sOverrideName) != setDedicatedKeys.end())
      {
         string strErrorMsg = " Error - --" + sOverrideName
            + " is not supported as a generic parameter override; use "
            + DedicatedHint(sOverrideName) + " instead.\n";
         logerr(strErrorMsg);
         exit(1);
      }

         if (setParamKeys.find(sOverrideName) == setParamKeys.end())
         {
            string strErrorMsg = " Error - unknown option \"--" + sOverrideName
               + "\"; key is not present in params file \"" + sParamsFile + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }
   }

   LoadParameters(szParamsFile, pSearchMgr, vCliParamOverrides);

   int iNumThreadsFromParams = 0;
   if (!pSearchMgr->GetParamValue("num_threads", iNumThreadsFromParams))
      iNumThreadsFromParams = 0;

   if (bSetFirstScan || bSetLastScan)
   {
      IntRange iScanRange;
      pSearchMgr->GetParamValue("scan_range", iScanRange);
      if (bSetFirstScan)
         iScanRange.iStart = iFirstScan;
      if (bSetLastScan)
         iScanRange.iEnd = iLastScan;

      char szParamStringVal[512];
      snprintf(szParamStringVal, sizeof(szParamStringVal), "%d %d", iScanRange.iStart, iScanRange.iEnd);
      pSearchMgr->SetParam("scan_range", szParamStringVal, iScanRange);
   }

   if (bSetBatchSize)
   {
      char szParamStringVal[512];
      snprintf(szParamStringVal, sizeof(szParamStringVal), "%d", iBatchSize);
      pSearchMgr->SetParam("spectrum_batch_size", szParamStringVal, iBatchSize);
   }

   if (bSetThreadOverride)
   {
      char szParamStringVal[512];
      snprintf(szParamStringVal, sizeof(szParamStringVal), "%d", iThreadOverride);
      pSearchMgr->SetParam("num_threads", szParamStringVal, iThreadOverride);
   }

   if (bCreateFragmentIndex)
   {
      pSearchMgr->SetParam("create_fragment_index", "1", 1);
      pSearchMgr->SetParam("create_peptide_index", "0", 0);
   }
   else if (bCreatePeptideIndex)
   {
      pSearchMgr->SetParam("create_fragment_index", "0", 0);
      pSearchMgr->SetParam("create_peptide_index", "1", 1);
   }

   if (!novelOpts.sOutputInternalNovelPeptidePath.empty() && !novelOpts.HasNovelInputOptions())
   {
      string strErrorMsg = " Error - --output_internal_novel_peptide requires --novel_protein and/or --novel_peptide.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (novelOpts.HasInternalNovelInput() && novelOpts.HasNovelInputOptions())
   {
      string strErrorMsg = " Error - --internal_novel_peptide cannot be used together with --novel_protein or --novel_peptide.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (novelOpts.bStopAfterSavingNovelPeptide && novelOpts.sOutputInternalNovelPeptidePath.empty())
   {
      string strErrorMsg = " Error - --stop-after-saving-novel-peptide requires --output_internal_novel_peptide.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (bSetOutputName && bOutputFolderSpecified && PathHasSeparatorForName(sOutputName))
   {
      string strErrorMsg = " Error - when using --name with --output-folder, --name must be a base name without path separators.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if ((bCreateFragmentIndex || bCreatePeptideIndex)
         && (novelOpts.HasNovelMode() || novelOpts.HasExplicitScanFilter()
            || !novelOpts.sOutputInternalNovelPeptidePath.empty()))
   {
      string strErrorMsg = " Error - novel/scan subset options cannot be used with -i or -j.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (novelOpts.sOutputFolder.empty())
      novelOpts.sOutputFolder = ".";

   {
      string sErrorMsg;
      if (!EnsureDirectoryExistsRecursive(novelOpts.sOutputFolder, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
   }
   SetCometPlusTempDirectory(novelOpts.sOutputFolder);

   if (!novelOpts.sOutputInternalNovelPeptidePath.empty())
      sResolvedOutputInternalNovelPath = ResolveInternalOutputPath(novelOpts.sOutputInternalNovelPeptidePath,
                                                                   novelOpts.sOutputFolder);

   bool bNeedInputSpectra = !novelOpts.bStopAfterSavingNovelPeptide;
   if ((novelOpts.HasNovelMode() || novelOpts.HasExplicitScanFilter()) && vInputArgs.empty() && bNeedInputSpectra)
   {
      string strErrorMsg = " Error - at least one input spectrum file is required when using novel or scan subset options.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   vector<string> vKnownDatabases = vDatabases;
   if (vKnownDatabases.empty())
   {
      string sDatabaseFromParams;
      if (pSearchMgr->GetParamValue("database_name", sDatabaseFromParams)
            && !sDatabaseFromParams.empty())
      {
         vKnownDatabases.push_back(sDatabaseFromParams);
      }
   }

   bool bKnownAllIdx = true;
   bool bKnownAllFasta = true;
   int iKnownIdxType = 0;
   for (size_t i = 0; i < vKnownDatabases.size(); ++i)
   {
      bool bIdx = IsIdxDatabasePath(vKnownDatabases.at(i));
      if (bIdx)
      {
         bKnownAllFasta = false;
      }
      else
      {
         bKnownAllIdx = false;
      }
   }
   if (!vKnownDatabases.empty() && !(bKnownAllIdx || bKnownAllFasta))
   {
      string strErrorMsg = " Error - mixed database types are not allowed; provide either all FASTA files or all .idx files.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   if (bKnownAllIdx)
   {
      string sProbeError;
      for (size_t i = 0; i < vKnownDatabases.size(); ++i)
      {
         int iObservedType = 0;
         if (!CometPlusProbeIdxType(vKnownDatabases.at(i), iObservedType, sProbeError))
         {
            logerr(sProbeError);
            exit(1);
         }

         if (i == 0)
            iKnownIdxType = iObservedType;
         else if (iObservedType != iKnownIdxType)
         {
            string strErrorMsg = " Error - mixed .idx types are not supported in one search invocation.\n";
            logerr(strErrorMsg);
            exit(1);
         }
      }
   }

   vector<InputFileInfo*> vParsedInputs;
   for (size_t i = 0; i < vInputArgs.size(); ++i)
   {
      InputFileInfo* pInputFileInfo = new InputFileInfo();
      std::vector<char> vInput(vInputArgs.at(i).begin(), vInputArgs.at(i).end());
      vInput.push_back('\0');

      if (!ParseCmdLine(vInput.data(), pInputFileInfo, pSearchMgr))
      {
         string strErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n"
            + " Error - input file \"" + std::string(pInputFileInfo->szFileName) + "\" not found.\n";
         logerr(strErrorMsg);
         pvInputFiles.clear();
         exit(1);
      }

      vParsedInputs.push_back(pInputFileInfo);
   }

   const bool bNovelMode = novelOpts.HasNovelMode();
   const bool bUseNovelMergedSearch = (bNovelMode && vParsedInputs.size() > 1);
   bool bRunCometEach = false;
   int iRunCometEachTotalThreads = 0;

   if (bRunCometEachRequested)
   {
      vector<string> vDisableReasons;
      if (!bNovelMode)
         vDisableReasons.push_back("novel mode is off");
      if (vParsedInputs.size() <= 1)
         vDisableReasons.push_back("input_count<=1");
      if (!bKnownAllIdx)
         vDisableReasons.push_back("known database is not all .idx");
      if (iKnownIdxType != 2)
         vDisableReasons.push_back("known .idx type is not peptide index (-j)");

      if (!vDisableReasons.empty())
      {
         string sReasonText;
         for (size_t i = 0; i < vDisableReasons.size(); ++i)
         {
            if (i > 0)
               sReasonText += "; ";
            sReasonText += vDisableReasons.at(i);
         }

         char szWarnBuf[2048];
         snprintf(szWarnBuf,
                  sizeof(szWarnBuf),
                  " Warning - --run-comet-each requested but disabled (%s). Falling back to merged-MGF single-search workflow.\n",
                  sReasonText.c_str());
         logout(szWarnBuf);
      }
      else
      {
         bRunCometEach = true;
         if (bSetThreadOverride && iThreadOverride > 0)
            iRunCometEachTotalThreads = iThreadOverride;
         else if (iNumThreadsFromParams > 0)
            iRunCometEachTotalThreads = iNumThreadsFromParams;

         if (iRunCometEachTotalThreads <= 0)
         {
            string strErrorMsg = " Error - --run-comet-each requires a positive total thread count: set --thread > 0 or num_threads > 0 in params.\n";
            logerr(strErrorMsg);
            exit(1);
         }
      }
   }

   if (bSetOutputName && vParsedInputs.size() != 1 && !bUseNovelMergedSearch)
   {
      string strErrorMsg = " Error - --name/-N can be used only when one input spectrum file is provided.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   string sNovelMergedOutputBase;
   if (bUseNovelMergedSearch)
   {
      string sMergedBaseStem = bSetOutputName ? sOutputName : "cometplus_novel_merged";
      if (sMergedBaseStem.empty())
      {
         string strErrorMsg = " Error - cannot derive merged output basename for novel multi-input search.\n";
         logerr(strErrorMsg);
         exit(1);
      }

      sNovelMergedOutputBase = JoinPathLocal(novelOpts.sOutputFolder, sMergedBaseStem);
      if (sNovelMergedOutputBase.length() >= SIZE_FILE)
      {
         string strErrorMsg = " Error - output basename path is too long: \"" + sNovelMergedOutputBase + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }
   }

   for (size_t i = 0; i < vParsedInputs.size(); ++i)
   {
      string sInputBaseStem;
      if (bUseNovelMergedSearch)
         sInputBaseStem = ComputeInputBaseName(vParsedInputs.at(i)->szFileName);
      else if (bSetOutputName)
         sInputBaseStem = sOutputName;
      else
         sInputBaseStem = ComputeInputBaseName(vParsedInputs.at(i)->szFileName);

      if (sInputBaseStem.empty())
      {
         string strErrorMsg = " Error - cannot derive output basename for input file \"" + string(vParsedInputs.at(i)->szFileName) + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }

      string sOutputBase = JoinPathLocal(novelOpts.sOutputFolder, sInputBaseStem);
      if (sOutputBase.length() >= SIZE_FILE)
      {
         string strErrorMsg = " Error - output basename path is too long: \"" + sOutputBase + "\".\n";
         logerr(strErrorMsg);
         exit(1);
      }

      strncpy(vParsedInputs.at(i)->szBaseName, sOutputBase.c_str(), SIZE_FILE - 1);
      vParsedInputs.at(i)->szBaseName[SIZE_FILE - 1] = '\0';
   }

   int iDecoySearch = 0;
   pSearchMgr->GetParamValue("decoy_search", iDecoySearch);

   int iOutputInt = 0;
   bool bOutputSqtFile = false;
   bool bOutputTxtFile = false;
   bool bOutputPepXmlFile = false;
   int iOutputMzidFile = 0;
   bool bOutputPercolatorFile = false;

   if (pSearchMgr->GetParamValue("output_sqtfile", iOutputInt))
      bOutputSqtFile = (iOutputInt != 0);
   if (pSearchMgr->GetParamValue("output_txtfile", iOutputInt))
      bOutputTxtFile = (iOutputInt != 0);
   if (pSearchMgr->GetParamValue("output_pepxmlfile", iOutputInt))
      bOutputPepXmlFile = (iOutputInt != 0);
   if (pSearchMgr->GetParamValue("output_mzidentmlfile", iOutputInt))
      iOutputMzidFile = iOutputInt;
   if (pSearchMgr->GetParamValue("output_percolatorfile", iOutputInt))
      bOutputPercolatorFile = (iOutputInt != 0);

   string sOutputSuffix;
   pSearchMgr->GetParamValue("output_suffix", sOutputSuffix);
   string sTxtExt = "txt";
   if (pSearchMgr->GetParamValue("text_file_extension", sTxtExt))
   {
      if (sTxtExt.empty())
         sTxtExt = "txt";
   }

   if (!(bCreateFragmentIndex || bCreatePeptideIndex))
   {
      bool bPlanOutputSqtFile = bOutputSqtFile;
      bool bPlanOutputTxtFile = bOutputTxtFile;
      bool bPlanOutputPepXmlFile = bOutputPepXmlFile;
      int iPlanOutputMzidFile = iOutputMzidFile;
      bool bPlanOutputPercolatorFile = bOutputPercolatorFile;

      if (bRunCometEach)
      {
         // run-comet-each mode guarantees merged pin output only.
         bPlanOutputSqtFile = false;
         bPlanOutputTxtFile = false;
         bPlanOutputPepXmlFile = false;
         iPlanOutputMzidFile = 0;
         bPlanOutputPercolatorFile = true;
      }

      struct PlannedOutputInput
      {
         InputFileInfo inputInfo;
         string sBaseName;
      };

      vector<PlannedOutputInput> vPlannedOutputInputs;
      if (bUseNovelMergedSearch)
      {
         PlannedOutputInput mergedOutputInput;
         mergedOutputInput.inputInfo.iAnalysisType = AnalysisType_EntireFile;
         mergedOutputInput.sBaseName = sNovelMergedOutputBase;
         vPlannedOutputInputs.push_back(mergedOutputInput);
      }
      else
      {
         for (size_t i = 0; i < vParsedInputs.size(); ++i)
         {
            PlannedOutputInput outputInput;
            outputInput.inputInfo = *vParsedInputs.at(i);
            outputInput.sBaseName = vParsedInputs.at(i)->szBaseName;
            vPlannedOutputInputs.push_back(outputInput);
         }
      }

      vector<string> vPlannedOutputs;
      for (size_t i = 0; i < vPlannedOutputInputs.size(); ++i)
      {
         AppendPlannedOutputFilesForInput(vPlannedOutputInputs.at(i).inputInfo,
                                          vPlannedOutputInputs.at(i).sBaseName,
                                          sOutputSuffix,
                                          sTxtExt,
                                          bPlanOutputSqtFile,
                                          bPlanOutputTxtFile,
                                          bPlanOutputPepXmlFile,
                                          iPlanOutputMzidFile,
                                          bPlanOutputPercolatorFile,
                                          iDecoySearch,
                                          vPlannedOutputs);
      }

      if (!sResolvedOutputInternalNovelPath.empty())
         vPlannedOutputs.push_back(sResolvedOutputInternalNovelPath);

      unordered_map<string, vector<string>> mPathToSources;
      unordered_set<string> setExistingSeen;
      vector<string> vExisting;
      for (size_t i = 0; i < vPlannedOutputs.size(); ++i)
      {
         string sPathKey = NormalizePathKeyLocal(vPlannedOutputs.at(i));
         mPathToSources[sPathKey].push_back(vPlannedOutputs.at(i));
         if (FileExistsLocal(vPlannedOutputs.at(i)) && setExistingSeen.insert(sPathKey).second)
            vExisting.push_back(vPlannedOutputs.at(i));
      }

      vector<string> vInternalConflicts;
      for (auto it = mPathToSources.begin(); it != mPathToSources.end(); ++it)
      {
         if (it->second.size() > 1)
         {
            string sMsg = it->second.at(0) + " (planned " + std::to_string(it->second.size()) + " times)";
            vInternalConflicts.push_back(sMsg);
         }
      }

      if (!vInternalConflicts.empty() || !vExisting.empty())
      {
         string strErrorMsg = " Error - output file conflict(s) detected before search.\n";
         if (!vInternalConflicts.empty())
         {
            strErrorMsg += "  Internal planned conflicts:\n";
            for (size_t i = 0; i < vInternalConflicts.size(); ++i)
               strErrorMsg += "   - " + vInternalConflicts.at(i) + "\n";
         }
         if (!vExisting.empty())
         {
            strErrorMsg += "  Existing files on disk:\n";
            for (size_t i = 0; i < vExisting.size(); ++i)
               strErrorMsg += "   - " + vExisting.at(i) + "\n";
         }
         logerr(strErrorMsg);
         exit(1);
      }
   }

   int iEqualIL = 1;
   pSearchMgr->GetParamValue("equal_I_and_L", iEqualIL);
   bool bTreatSameIL = (iEqualIL != 0);

   vector<double> vNovelMasses;
   string sNovelScoringDatabase;
   vector<NovelPeptideRecord> vNovelRecords;

   RunNovelWorkflowIfNeeded(novelOpts,
                            pSearchMgr,
                            vKnownDatabases,
                            bKnownAllIdx,
                            iKnownIdxType,
                            sParamsFile,
                            iThreadOverride,
                            bTreatSameIL,
                            iDecoySearch,
                            sResolvedOutputInternalNovelPath,
                            tProgramStart,
                            vTempArtifacts,
                            sMergedDatabasePath,
                            vNovelRecords,
                            vNovelMasses,
                            sNovelScoringDatabase);

   vector<string> vFinalDatabases = vKnownDatabases;
   if (!sNovelScoringDatabase.empty())
      vFinalDatabases.push_back(sNovelScoringDatabase);

   if (!vFinalDatabases.empty())
   {
      string sErrorMsg;
      if (!ConfigureDatabaseInputs(vFinalDatabases, pSearchMgr, sMergedDatabasePath, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
      if (!sMergedDatabasePath.empty())
         vTempArtifacts.push_back(sMergedDatabasePath);
   }

   const bool bEnableNovelOutputOnly = (bNovelMode || bForceNovelOutputOnly);
   char szNovelOutputOnlyParam[32];
   snprintf(szNovelOutputOnlyParam, sizeof(szNovelOutputOnlyParam), "%d", bEnableNovelOutputOnly ? 1 : 0);
   pSearchMgr->SetParam("cometplus_novel_output_only", szNovelOutputOnlyParam, bEnableNovelOutputOnly ? 1 : 0);
   CometPlusSetNovelOutputOnly(bEnableNovelOutputOnly);

   PrepareInputFilesWithPrefilter(vParsedInputs,
                                  pvInputFiles,
                                  pSearchMgr,
                                  novelOpts,
                                  bSetThreadOverride,
                                  iThreadOverride,
                                  bNovelMode,
                                  vNovelMasses,
                                  bUseNovelMergedSearch,
                                  bRunCometEach,
                                  sNovelMergedOutputBase,
                                  bCreateFragmentIndex,
                                  bCreatePeptideIndex,
                                  tProgramStart,
                                  vTempArtifacts);

   if (!bRunCometEach)
      return;

   if (pvInputFiles.empty())
   {
      string strErrorMsg = " Error - --run-comet-each did not receive any filtered MGF shards.\n";
      logerr(strErrorMsg);
      exit(1);
   }

   struct RunEachSourceShard
   {
      string sInputMgfPath;
      unsigned long long ullSizeBytes = 0;
   };

   struct RunEachTask
   {
      string sInputMgfPath;
      int iTaskThreads = 1;
      string sShardBaseName;
      string sShardPinPath;
      bool bSuccess = false;
      string sErrorMsg;
   };

   auto ReadFileSizeBytes = [&](const string& sPath, unsigned long long& ullSizeBytes, string& sErrorMsg) -> bool
   {
      ullSizeBytes = 0;
      sErrorMsg.clear();
      std::ifstream inFile(sPath.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
      if (!inFile.good())
      {
         sErrorMsg = " Error - cannot read filtered MGF shard \"" + sPath + "\" for --run-comet-each grouping.\n";
         return false;
      }

      std::ifstream::pos_type iEndPos = inFile.tellg();
      if (iEndPos < 0)
      {
         sErrorMsg = " Error - cannot determine file size for filtered MGF shard \"" + sPath + "\".\n";
         return false;
      }

      ullSizeBytes = static_cast<unsigned long long>(iEndPos);
      return true;
   };

   vector<RunEachSourceShard> vSourceShards;
   vSourceShards.reserve(pvInputFiles.size());
   for (size_t i = 0; i < pvInputFiles.size(); ++i)
   {
      RunEachSourceShard shard;
      shard.sInputMgfPath = pvInputFiles.at(i)->szFileName;
      string sFileSizeError;
      if (!ReadFileSizeBytes(shard.sInputMgfPath, shard.ullSizeBytes, sFileSizeError))
      {
         logerr(sFileSizeError);
         exit(1);
      }
      vSourceShards.push_back(shard);
   }

   const size_t tSourceShardCount = vSourceShards.size();
   int iTargetTaskCount = iRunCometEachTotalThreads / 2;
   if (iTargetTaskCount < 1)
      iTargetTaskCount = 1;
   if ((size_t)iTargetTaskCount > tSourceShardCount)
      iTargetTaskCount = (int)tSourceShardCount;

   vector<vector<RunEachSourceShard>> vTaskGroups((size_t)iTargetTaskCount);
   vector<unsigned long long> vTaskGroupBytes((size_t)iTargetTaskCount, 0);

   if ((size_t)iTargetTaskCount == tSourceShardCount)
   {
      for (size_t i = 0; i < vSourceShards.size(); ++i)
      {
         vTaskGroups.at(i).push_back(vSourceShards.at(i));
         vTaskGroupBytes.at(i) = vSourceShards.at(i).ullSizeBytes;
      }
   }
   else
   {
      vector<RunEachSourceShard> vSortedShards = vSourceShards;
      std::sort(vSortedShards.begin(),
                vSortedShards.end(),
                [](const RunEachSourceShard& a, const RunEachSourceShard& b)
                {
                   if (a.ullSizeBytes != b.ullSizeBytes)
                      return a.ullSizeBytes > b.ullSizeBytes;
                   return a.sInputMgfPath < b.sInputMgfPath;
                });

      for (size_t i = 0; i < vSortedShards.size(); ++i)
      {
         size_t iBestGroup = 0;
         for (size_t iGroup = 1; iGroup < vTaskGroups.size(); ++iGroup)
         {
            if (vTaskGroupBytes.at(iGroup) < vTaskGroupBytes.at(iBestGroup))
               iBestGroup = iGroup;
            else if (vTaskGroupBytes.at(iGroup) == vTaskGroupBytes.at(iBestGroup)
                  && vTaskGroups.at(iGroup).size() < vTaskGroups.at(iBestGroup).size())
               iBestGroup = iGroup;
         }

         vTaskGroups.at(iBestGroup).push_back(vSortedShards.at(i));
         vTaskGroupBytes.at(iBestGroup) += vSortedShards.at(i).ullSizeBytes;
      }
   }

   vector<string> vRunEachInputMgfPaths;
   vector<unsigned long long> vRunEachInputMgfBytes;
   vRunEachInputMgfPaths.reserve(vTaskGroups.size());
   vRunEachInputMgfBytes.reserve(vTaskGroups.size());
   for (size_t iGroup = 0; iGroup < vTaskGroups.size(); ++iGroup)
   {
      const vector<RunEachSourceShard>& vGroup = vTaskGroups.at(iGroup);
      if (vGroup.empty())
      {
         string strErrorMsg = " Error - internal run-comet-each grouping produced an empty MGF task.\n";
         logerr(strErrorMsg);
         exit(1);
      }

      if (vGroup.size() == 1)
      {
         vRunEachInputMgfPaths.push_back(vGroup.at(0).sInputMgfPath);
         vRunEachInputMgfBytes.push_back(vTaskGroupBytes.at(iGroup));
      }
      else
      {
         vector<std::pair<string, string>> vMergeInputs;
         vMergeInputs.reserve(vGroup.size());
         for (size_t iShard = 0; iShard < vGroup.size(); ++iShard)
            vMergeInputs.push_back(std::make_pair(vGroup.at(iShard).sInputMgfPath, ""));

         string sMergedTaskMgfPath;
         string sMergeMgfError;
         if (!MergeFilteredMgfFilesWithSourceTag(vMergeInputs, sMergedTaskMgfPath, sMergeMgfError))
         {
            logerr(sMergeMgfError);
            exit(1);
         }
         vTempArtifacts.push_back(sMergedTaskMgfPath);
         vRunEachInputMgfPaths.push_back(sMergedTaskMgfPath);
         vRunEachInputMgfBytes.push_back(vTaskGroupBytes.at(iGroup));

         char szGroupMergeLog[2048];
         snprintf(szGroupMergeLog,
                  sizeof(szGroupMergeLog),
                  " [%s] run-comet-each grouped merge: group=%zu, source_shards=%zu, bytes=%llu -> \"%s\"\n",
                  GetLocalTimestampString().c_str(),
                  iGroup + 1,
                  vGroup.size(),
                  vTaskGroupBytes.at(iGroup),
                  sMergedTaskMgfPath.c_str());
         logout(szGroupMergeLog);
      }
   }

   const size_t tNumTasks = vRunEachInputMgfPaths.size();
   int iMaxParallelJobs = (int)tNumTasks;
   if (iMaxParallelJobs < 1)
      iMaxParallelJobs = 1;

   vector<int> vTaskThreads(tNumTasks, 1);
   int iBaseTaskThreads = iRunCometEachTotalThreads / (int)tNumTasks;
   int iExtraThreadTasks = iRunCometEachTotalThreads % (int)tNumTasks;
   if (iBaseTaskThreads < 1)
   {
      iBaseTaskThreads = 1;
      iExtraThreadTasks = 0;
   }

   vector<size_t> vTaskOrder;
   vTaskOrder.reserve(tNumTasks);
   for (size_t i = 0; i < tNumTasks; ++i)
      vTaskOrder.push_back(i);

   std::sort(vTaskOrder.begin(),
             vTaskOrder.end(),
             [&](size_t a, size_t b)
             {
                if (vRunEachInputMgfBytes.at(a) != vRunEachInputMgfBytes.at(b))
                   return vRunEachInputMgfBytes.at(a) > vRunEachInputMgfBytes.at(b);
                return a < b;
             });

   int iThreadAssignedSum = 0;
   int iMinTaskThreads = INT_MAX;
   int iMaxTaskThreads = 0;
   for (size_t iOrder = 0; iOrder < vTaskOrder.size(); ++iOrder)
   {
      size_t iTask = vTaskOrder.at(iOrder);
      int iTaskThreads = iBaseTaskThreads;
      if ((int)iOrder < iExtraThreadTasks)
         ++iTaskThreads;
      if (iTaskThreads < 1)
         iTaskThreads = 1;

      vTaskThreads.at(iTask) = iTaskThreads;
      iThreadAssignedSum += iTaskThreads;
      if (iTaskThreads < iMinTaskThreads)
         iMinTaskThreads = iTaskThreads;
      if (iTaskThreads > iMaxTaskThreads)
         iMaxTaskThreads = iTaskThreads;
   }

   vector<RunEachTask> vTasks(tNumTasks);
   for (size_t i = 0; i < tNumTasks; ++i)
   {
      vTasks.at(i).sInputMgfPath = vRunEachInputMgfPaths.at(i);
      vTasks.at(i).iTaskThreads = vTaskThreads.at(i);

      string sTempMarkerPath;
      string sErrorMsg;
      if (!CreateTempPath("cometplus_run_each_shard", "", sTempMarkerPath, sErrorMsg))
      {
         logerr(sErrorMsg);
         exit(1);
      }
      vTempArtifacts.push_back(sTempMarkerPath);

      string sShardBaseName = ComputeInputBaseName(sTempMarkerPath);
      if (sShardBaseName.empty())
      {
         string strErrorMsg = " Error - failed to derive shard output basename for --run-comet-each.\n";
         logerr(strErrorMsg);
         exit(1);
      }
      vTasks.at(i).sShardBaseName = sShardBaseName;
      vTasks.at(i).sShardPinPath = JoinPathLocal(novelOpts.sOutputFolder,
                                                 sShardBaseName + sOutputSuffix + ".pin");
      vTempArtifacts.push_back(vTasks.at(i).sShardPinPath);
   }

   char szRunEachStart[1024];
   snprintf(szRunEachStart,
            sizeof(szRunEachStart),
            " [%s] run-comet-each: source_shards=%zu, grouped_tasks=%zu, total_threads=%d, assigned_threads=%d, task_threads=[%d,%d], max_parallel_jobs=%d\n",
            GetLocalTimestampString().c_str(),
            tSourceShardCount,
            tNumTasks,
            iRunCometEachTotalThreads,
            iThreadAssignedSum,
            iMinTaskThreads == INT_MAX ? 0 : iMinTaskThreads,
            iMaxTaskThreads,
            iMaxParallelJobs);
   logout(szRunEachStart);

   auto tRunEachStart = std::chrono::steady_clock::now();
   std::atomic<size_t> iNextTask(0);
   std::atomic<bool> bRunEachError(false);
   std::mutex mRunEachError;
   string sFirstRunEachError;
   vector<std::thread> vWorkers;
   vWorkers.reserve((size_t)iMaxParallelJobs);
   for (int iWorker = 0; iWorker < iMaxParallelJobs; ++iWorker)
   {
      vWorkers.emplace_back([&]()
      {
         while (true)
         {
            size_t iTask = iNextTask.fetch_add(1);
            if (iTask >= tNumTasks)
               break;

            if (bRunEachError.load())
               break;

            RunEachTask& task = vTasks.at(iTask);
            vector<string> vCmdArgs;
            vCmdArgs.push_back(g_sCometPlusExecutablePath);
            vCmdArgs.push_back("--params");
            vCmdArgs.push_back(sParamsFile);
            for (size_t iOverride = 0; iOverride < vCliParamOverrides.size(); ++iOverride)
            {
               vCmdArgs.push_back("--" + vCliParamOverrides.at(iOverride).sName);
               vCmdArgs.push_back(vCliParamOverrides.at(iOverride).sValue);
            }
            vCmdArgs.push_back("--thread");
            vCmdArgs.push_back(std::to_string(task.iTaskThreads));
            vCmdArgs.push_back("--output-folder");
            vCmdArgs.push_back(novelOpts.sOutputFolder);
            vCmdArgs.push_back("--name");
            vCmdArgs.push_back(task.sShardBaseName);
            vCmdArgs.push_back("--output_percolatorfile");
            vCmdArgs.push_back("1");
            vCmdArgs.push_back("--output_sqtfile");
            vCmdArgs.push_back("0");
            vCmdArgs.push_back("--output_txtfile");
            vCmdArgs.push_back("0");
            vCmdArgs.push_back("--output_pepxmlfile");
            vCmdArgs.push_back("0");
            vCmdArgs.push_back("--output_mzidentmlfile");
            vCmdArgs.push_back("0");
            vCmdArgs.push_back("--cometplus-novel-output-only");

            for (size_t iDb = 0; iDb < vFinalDatabases.size(); ++iDb)
            {
               vCmdArgs.push_back("--database");
               vCmdArgs.push_back(vFinalDatabases.at(iDb));
            }

            vCmdArgs.push_back(task.sInputMgfPath);

            string sCmdError;
            if (!RunExternalCommand(vCmdArgs, sCmdError))
            {
               task.bSuccess = false;
               task.sErrorMsg = " Error - --run-comet-each task search failed for input \"" + task.sInputMgfPath + "\".\n" + sCmdError;
               bRunEachError.store(true);
               std::lock_guard<std::mutex> lk(mRunEachError);
               if (sFirstRunEachError.empty())
                  sFirstRunEachError = task.sErrorMsg;
               continue;
            }

            task.bSuccess = true;
         }
      });
   }

   for (size_t i = 0; i < vWorkers.size(); ++i)
      vWorkers.at(i).join();

   if (bRunEachError.load())
   {
      if (sFirstRunEachError.empty())
         sFirstRunEachError = " Error - --run-comet-each task search failed.\n";
      logerr(sFirstRunEachError);
      exit(1);
   }

   vector<string> vShardPinFiles;
   vShardPinFiles.reserve(vTasks.size());
   for (size_t i = 0; i < vTasks.size(); ++i)
      vShardPinFiles.push_back(vTasks.at(i).sShardPinPath);

   string sMergedPinPath = sNovelMergedOutputBase + sOutputSuffix + ".pin";
   string sMergeError;
   if (!MergePercolatorPinFiles(vShardPinFiles, sMergedPinPath, sMergeError))
   {
      logerr(sMergeError);
      exit(1);
   }

   char szRunEachDone[2048];
   snprintf(szRunEachDone,
            sizeof(szRunEachDone),
            " [%s] run-comet-each pin merge: %zu task pin files -> \"%s\"\n",
            GetLocalTimestampString().c_str(),
            vShardPinFiles.size(),
            sMergedPinPath.c_str());
   logout(szRunEachDone);

   LogStageTiming("run-comet-each", tRunEachStart, tProgramStart);

   for (size_t i = 0; i < pvInputFiles.size(); ++i)
   {
      delete pvInputFiles.at(i);
      pvInputFiles.at(i) = NULL;
   }
   pvInputFiles.clear();

   bSearchAlreadyRun = true;
   bSearchSucceeded = true;
} // ProcessCmdLine
